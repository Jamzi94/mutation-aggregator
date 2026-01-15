#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import sqlite3
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass, asdict
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import requests
import yaml
from requests.adapters import HTTPAdapter
from requests.exceptions import RequestException
from urllib3.util.retry import Retry
from tqdm import tqdm

DEFAULT_CONFIG_PATH = os.environ.get("MUT_AGG_CONFIG", "config.yaml")
CACHE_DB = os.environ.get("MUT_AGG_CACHE", "mutation_cache.sqlite3")
LOG_FILE = os.environ.get("MUT_AGG_LOG", "mutation_aggregator.log")
MAX_WORKERS = 6

logger = logging.getLogger("mutation_aggregator")
logger.setLevel(logging.DEBUG)
_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
_stream_h = logging.StreamHandler()
_stream_h.setFormatter(_formatter)
_stream_h.setLevel(logging.INFO)
_file_h = logging.FileHandler(LOG_FILE)
_file_h.setFormatter(_formatter)
_file_h.setLevel(logging.DEBUG)
logger.addHandler(_stream_h)
logger.addHandler(_file_h)


@dataclass
class MutationRecord:
    position_aa: Optional[int]
    ref_residue: Optional[str]
    alt_residue: Optional[str]
    modification_type: Optional[str]
    impact_category: Optional[str]
    impact_score: Optional[float]
    source_databases: List[str]
    validation_status: Optional[str]
    clinvar_id: Optional[str]
    uniprot_variant_id: Optional[str]
    pubmed_references: List[Dict[str, Any]]
    structural_context: Optional[str]
    raw_sources: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["source_databases"] = list(self.source_databases)
        return d


def load_config(path: str = DEFAULT_CONFIG_PATH) -> Dict[str, Any]:
    if not os.path.exists(path):
        logger.warning("Config file %s not found — using defaults", path)
        return {}
    try:
        with open(path) as fh:
            cfg = yaml.safe_load(fh) or {}
            logger.info("Loaded configuration from %s", path)
            return cfg
    except Exception as e:
        logger.exception("Failed to read config %s: %s", path, e)
        return {}


class SQLiteCache:
    def __init__(self, path: str = CACHE_DB):
        self.path = path
        self._lock = threading.RLock()
        self._ensure_db()

    def _ensure_db(self) -> None:
        with self._get_conn() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS cache (
                    key TEXT PRIMARY KEY,
                    url TEXT,
                    params TEXT,
                    response TEXT,
                    timestamp REAL
                )
                """
            )
            conn.commit()

    @contextmanager
    def _get_conn(self):
        conn = sqlite3.connect(self.path, check_same_thread=False)
        try:
            yield conn
        finally:
            conn.close()

    def _make_key(self, url: str, params: Optional[Dict[str, Any]]) -> str:
        obj = {"url": url, "params": params or {}}
        j = json.dumps(obj, sort_keys=True, default=str)
        return hashlib.sha256(j.encode()).hexdigest()

    def get(self, url: str, params: Optional[Dict[str, Any]] = None) -> Optional[Any]:
        key = self._make_key(url, params)
        with self._lock, self._get_conn() as conn:
            cur = conn.execute("SELECT response, timestamp FROM cache WHERE key = ?", (key,))
            row = cur.fetchone()
            if row:
                resp_text, ts = row
                logger.debug("Cache hit for %s (age=%.1fs)", url, time.time() - ts)
                try:
                    return json.loads(resp_text)
                except json.JSONDecodeError:
                    return resp_text
            logger.debug("Cache miss for %s", url)
            return None

    def set(self, url: str, params: Optional[Dict[str, Any]], response: Any) -> None:
        key = self._make_key(url, params)
        resp_text = json.dumps(response, default=str)
        ts = time.time()
        with self._lock, self._get_conn() as conn:
            conn.execute(
                "REPLACE INTO cache (key, url, params, response, timestamp) VALUES (?, ?, ?, ?, ?)",
                (key, url, json.dumps(params, default=str), resp_text, ts),
            )
            conn.commit()
            logger.debug("Cached response for %s", url)


class HTTPClient:
    def __init__(self, config: Dict[str, Any]):
        self.session = requests.Session()
        retries = config.get("retries", 3)
        backoff_factor = config.get("backoff_factor", 0.5)
        status_forcelist = config.get("status_forcelist", [429, 500, 502, 503, 504])
        retry = Retry(
            total=retries,
            read=retries,
            connect=retries,
            backoff_factor=backoff_factor,
            status_forcelist=status_forcelist,
            allowed_methods=frozenset(["GET", "POST"]),
        )
        adapter = HTTPAdapter(max_retries=retry)
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)
        self.rate_limits = config.get("rate_limits", {})
        self._host_locks: Dict[str, threading.Lock] = {}
        self._host_last: Dict[str, float] = {}

    def _throttle(self, host: str) -> None:
        limit = self.rate_limits.get(host)
        if not limit:
            return
        lock = self._host_locks.setdefault(host, threading.Lock())
        with lock:
            last = self._host_last.get(host, 0.0)
            min_interval = 1.0 / float(limit) if limit > 0 else 0
            wait = min_interval - (time.time() - last)
            if wait > 0:
                logger.debug("Throttling %s: sleeping %.3fs", host, wait)
                time.sleep(wait)
            self._host_last[host] = time.time()

    def request(self, method: str, url: str, params: Optional[Dict[str, Any]] = None, **kwargs) -> Any:
        host = requests.utils.urlparse(url).netloc
        self._throttle(host)
        try:
            resp = self.session.request(method, url, params=params, timeout=20, **kwargs)
            resp.raise_for_status()
            try:
                return resp.json()
            except ValueError:
                return resp.text
        except RequestException as e:
            logger.warning("HTTP request failed for %s: %s", url, e)
            raise


class SourceAdapters:
    def __init__(self, client: HTTPClient, cache: SQLiteCache, config: Dict[str, Any]):
        self.client = client
        self.cache = cache
        self.config = config

    def fetch_uniprot(self, uniprot_id: str) -> Dict[str, Any]:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        cached = self.cache.get(url)
        if cached is not None:
            return cached
        try:
            data = self.client.request("GET", url)
            self.cache.set(url, None, data)
            return data
        except Exception as e:
            logger.error("Failed fetching UniProt %s: %s", uniprot_id, e)
            return {}

    def fetch_clinical(self, uniprot_id: str) -> List[Dict[str, Any]]:
        ncbi_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        term = uniprot_id
        search_url = f"{ncbi_base}/esearch.fcgi"
        params = {"db": "clinvar", "term": term, "retmode": "json"}
        cached = self.cache.get(search_url, params)
        results: List[Dict[str, Any]] = []
        try:
            data = cached if cached is not None else self.client.request("GET", search_url, params=params)
            if cached is None:
                self.cache.set(search_url, params, data)
            uid_list = []
            if isinstance(data, dict):
                uid_list = data.get("esearchresult", {}).get("idlist", [])
            for uid in uid_list[:50]:
                fetch_url = f"{ncbi_base}/efetch.fcgi"
                fetch_params = {"db": "clinvar", "id": uid, "retmode": "xml"}
                cached_f = self.cache.get(fetch_url, fetch_params)
                if cached_f is not None:
                    xml_text = cached_f
                else:
                    xml_text = self.client.request("GET", fetch_url, params=fetch_params)
                    self.cache.set(fetch_url, fetch_params, xml_text)
                results.append({"uid": uid, "xml": xml_text})
            return results
        except Exception as e:
            logger.warning("ClinVar fetch failed: %s", e)
            return []

    def fetch_references(self, uniprot_id: str) -> List[Dict[str, Any]]:
        ncbi_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        term = f"{uniprot_id}[All Fields]"
        search_url = f"{ncbi_base}/esearch.fcgi"
        params = {"db": "pubmed", "term": term, "retmode": "json", "retmax": 50}
        cached = self.cache.get(search_url, params)
        results: List[Dict[str, Any]] = []
        try:
            data = cached if cached is not None else self.client.request("GET", search_url, params=params)
            if cached is None:
                self.cache.set(search_url, params, data)
            uid_list = []
            if isinstance(data, dict):
                uid_list = data.get("esearchresult", {}).get("idlist", [])
            for uid in uid_list[:50]:
                fetch_url = f"{ncbi_base}/esummary.fcgi"
                fetch_params = {"db": "pubmed", "id": uid, "retmode": "json"}
                cached_f = self.cache.get(fetch_url, fetch_params)
                if cached_f is not None:
                    summary = cached_f
                else:
                    summary = self.client.request("GET", fetch_url, params=fetch_params)
                    self.cache.set(fetch_url, fetch_params, summary)
                try:
                    doc = summary.get("result", {}).get(str(uid), {})
                    doi = doc.get("elocationid")
                    title = doc.get("title")
                    journal = doc.get("fulljournalname")
                    year = doc.get("pubdate")
                    results.append({"uid": uid, "doi": doi, "title": title, "journal": journal, "year": year})
                except Exception:
                    results.append({"uid": uid, "raw": summary})
            return results
        except Exception as e:
            logger.warning("PubMed fetch failed: %s", e)
            return []


AA_REGEX = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|=)")
AA_1TO3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "E": "Glu", "Q": "Gln",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}
THREE_TO_ONE = {v.capitalize(): k for k, v in AA_1TO3.items()}


def parse_uniprot_variant(variant_feature: Dict[str, Any]) -> Optional[MutationRecord]:
    try:
        description = variant_feature.get("description") or ""
        pos = None
        ref = None
        alt = None
        
        # Extract position from location.start.value (UniProt format)
        location = variant_feature.get("location", {})
        start_obj = location.get("start", {})
        if isinstance(start_obj, dict):
            pos = start_obj.get("value")
        elif isinstance(start_obj, int):
            pos = start_obj
        
        # Extract ref/alt from alternativeSequence
        alt_seq = variant_feature.get("alternativeSequence", {})
        if alt_seq:
            ref = alt_seq.get("originalSequence")
            alt_list = alt_seq.get("alternativeSequences", [])
            alt = alt_list[0] if alt_list else None
        
        # Fallback to old format
        if not ref:
            ref = variant_feature.get("original")
        if not alt:
            variation = variant_feature.get("variation")
            if isinstance(variation, str):
                alt = variation
        
        # Fallback: parse from description
        if pos is None:
            m = AA_REGEX.search(description)
            if m:
                ref3, pos_str, alt3 = m.groups()
                pos = int(pos_str)
                if not ref:
                    ref = THREE_TO_ONE.get(ref3.capitalize(), ref3)
                if not alt:
                    alt = THREE_TO_ONE.get(alt3.capitalize(), alt3) if alt3 != "=" else None
        
        # Extract PubMed IDs from evidences
        pubmed_ids = []
        for ev in variant_feature.get("evidences", []):
            if ev.get("source") == "PubMed" and ev.get("id"):
                pubmed_ids.append(ev.get("id"))
        
        mod_type = variant_feature.get("type") or "variant"
        validation = "experimental" if variant_feature.get("evidences") else None
        
        return MutationRecord(
            position_aa=pos,
            ref_residue=ref,
            alt_residue=alt,
            modification_type=mod_type,
            impact_category=None,
            impact_score=None,
            source_databases=["UniProt"],
            validation_status=validation,
            clinvar_id=None,
            uniprot_variant_id=variant_feature.get("id"),
            pubmed_references=[{"pmid": pmid} for pmid in pubmed_ids],
            structural_context=description,
            raw_sources={"uniprot_feature": variant_feature},
        )
    except Exception as e:
        logger.debug("Failed to parse UniProt variant: %s", e)
        return None


def make_mutation_key(pos: Optional[int], ref: Optional[str], alt: Optional[str]) -> str:
    if pos is None:
        return f"unknown:{ref}:{alt}"
    return f"p.{ref or '?'}{pos}{alt or '?'}"


def extract_pos_ref_alt_from_xml(xml_text: Any) -> Tuple[Optional[int], Optional[str], Optional[str]]:
    try:
        if not xml_text:
            return None, None, None
        if isinstance(xml_text, dict):
            xml_text = json.dumps(xml_text)
        s = str(xml_text)
        m = AA_REGEX.search(s)
        if m:
            ref3, pos_str, alt3 = m.groups()
            pos = int(pos_str)
            ref = THREE_TO_ONE.get(ref3.capitalize(), ref3)
            alt = THREE_TO_ONE.get(alt3.capitalize(), alt3) if alt3 != "=" else None
            return pos, ref, alt
    except Exception:
        pass
    return None, None, None


def normalize_records_from_sources(uniprot_json: Dict[str, Any], clinical: List[Dict[str, Any]], references: List[Dict[str, Any]]) -> Dict[str, MutationRecord]:
    records: Dict[str, MutationRecord] = {}

    features = []
    try:
        features = uniprot_json.get("features", []) if isinstance(uniprot_json, dict) else []
    except Exception:
        features = []
    for feat in features:
        if feat.get("type") and feat.get("type").lower() in ("variant", "mutagenesis"):
            mr = parse_uniprot_variant(feat)
            if mr:
                key = make_mutation_key(mr.position_aa, mr.ref_residue, mr.alt_residue)
                if key in records:
                    rec = records[key]
                    if "UniProt" not in rec.source_databases:
                        rec.source_databases.append("UniProt")
                    rec.raw_sources.setdefault("uniprot_features", []).append(feat)
                else:
                    records[key] = mr

    for item in clinical:
        uid = item.get("uid") or item.get("id")
        xml = item.get("xml")
        pos, ref, alt = extract_pos_ref_alt_from_xml(xml)
        key = make_mutation_key(pos, ref, alt)
        if key not in records:
            records[key] = MutationRecord(
                position_aa=pos,
                ref_residue=ref,
                alt_residue=alt,
                modification_type="clinical",
                impact_category=None,
                impact_score=None,
                source_databases=["ClinVar"],
                validation_status="clinical",
                clinvar_id=uid,
                uniprot_variant_id=None,
                pubmed_references=[],
                structural_context=None,
                raw_sources={"clinvar": [item]},
            )
        else:
            rec = records[key]
            if "ClinVar" not in rec.source_databases:
                rec.source_databases.append("ClinVar")
            if "clinvar" not in rec.raw_sources:
                rec.raw_sources["clinvar"] = []
            rec.raw_sources["clinvar"].append(item)

    for ref in references:
        for rec in records.values():
            rec.pubmed_references.append(ref)

    for _, rec in records.items():
        if rec.validation_status == "experimental":
            rec.impact_category = rec.impact_category or "high"
        elif "ClinVar" in rec.source_databases:
            rec.impact_category = rec.impact_category or "medium"
        else:
            rec.impact_category = rec.impact_category or "low"
    return records


class MutationAggregator:
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        self.config = config or load_config()
        client_cfg = self.config.get("http_client", {})
        self.client = HTTPClient(client_cfg)
        cache_path = self.config.get("cache_db", CACHE_DB)
        self.cache = SQLiteCache(cache_path)
        self.adapters = SourceAdapters(self.client, self.cache, self.config)
        self.output_dir = self.config.get("output_dir", "outputs")
        os.makedirs(self.output_dir, exist_ok=True)

    def analyze_uniprot_mutations(self, uniprot_id: str, write_outputs: bool = True) -> pd.DataFrame:
        logger.info("Starting analysis for UniProt ID %s", uniprot_id)
        sources = {
            "uniprot": lambda: self.adapters.fetch_uniprot(uniprot_id),
            "clinical": lambda: self.adapters.fetch_clinical(uniprot_id),
            "references": lambda: self.adapters.fetch_references(uniprot_id),
        }

        fetched: Dict[str, Any] = {}
        with ThreadPoolExecutor(max_workers=self.config.get("max_workers", MAX_WORKERS)) as ex:
            futures = {ex.submit(func): name for name, func in sources.items()}
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Fetching sources"):
                name = futures[fut]
                try:
                    fetched[name] = fut.result()
                    logger.info("Fetched %s (items=%s)", name, (len(fetched[name]) if hasattr(fetched[name], '__len__') else '1'))
                except Exception as e:
                    logger.exception("Source %s failed: %s", name, e)
                    fetched[name] = [] if name != "uniprot" else {}

        records = normalize_records_from_sources(
            uniprot_json=fetched.get("uniprot", {}),
            clinical=fetched.get("clinical", []),
            references=fetched.get("references", []),
        )

        df_rows = []
        for key, rec in records.items():
            # Format mutation as REF→ALT
            mutation_str = ""
            if rec.ref_residue and rec.alt_residue:
                mutation_str = f"{rec.ref_residue}→{rec.alt_residue}"
            elif rec.ref_residue:
                mutation_str = f"{rec.ref_residue}→?"
            elif rec.alt_residue:
                mutation_str = f"?→{rec.alt_residue}"
            
            # Extract PubMed IDs
            pubmed_ids = []
            for ref in (rec.pubmed_references or []):
                if isinstance(ref, dict):
                    pmid = ref.get("pmid") or ref.get("uid")
                    if pmid:
                        pubmed_ids.append(str(pmid))
            
            df_rows.append(
                {
                    "uniprot_id": uniprot_id,
                    "position": rec.position_aa,
                    "mutation": mutation_str,
                    "type": rec.modification_type,
                    "impact": rec.impact_category,
                    "sources": ",".join(sorted(set(rec.source_databases))),
                    "validation": rec.validation_status,
                    "description": rec.structural_context,
                    "pubmed_ids": ",".join(sorted(set(pubmed_ids))) if pubmed_ids else None,
                    "clinvar_id": rec.clinvar_id,
                    "uniprot_variant_id": rec.uniprot_variant_id,
                }
            )

        df = pd.DataFrame(df_rows)
        metadata = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "uniprot_id": uniprot_id,
            "sources_requested": list(sources.keys()),
            "config": {k: v for k, v in self.config.items() if k != "api_keys"},
        }

        if write_outputs:
            out_prefix = os.path.join(self.output_dir, f"{uniprot_id}_mutations")
            json_path = out_prefix + ".json"
            xlsx_path = out_prefix + ".xlsx"
            try:
                export_obj = {"metadata": metadata, "records": [r for r in df.to_dict(orient="records")]}
                with open(json_path, "w") as fh:
                    json.dump(export_obj, fh, default=str, indent=2)
                df.to_excel(xlsx_path, index=False)
                logger.info("Wrote outputs: %s, %s", json_path, xlsx_path)
            except Exception as e:
                logger.exception("Failed to write output files: %s", e)

        return df


def main():
    parser = argparse.ArgumentParser(description="Aggregate mutation data for UniProt IDs")
    parser.add_argument("uniprot_id", nargs="?", help="UniProt accession (e.g., P38398)")
    parser.add_argument("--batch", help="CSV file with 'uniprot_id' column for batch processing")
    parser.add_argument("--config", default=DEFAULT_CONFIG_PATH, help="YAML config path")
    parser.add_argument("--no-write", action="store_true", help="Do not write outputs to disk")
    args = parser.parse_args()

    if not args.uniprot_id and not args.batch:
        parser.error("Either provide a uniprot_id or use --batch with a CSV file")

    cfg = load_config(args.config)
    agg = MutationAggregator(cfg)

    if args.batch:
        # Batch mode: read CSV and process each UniProt ID
        logger.info("Batch mode: reading UniProt IDs from %s", args.batch)
        try:
            batch_df = pd.read_csv(args.batch)
            if "uniprot_id" not in batch_df.columns:
                logger.error("CSV must have a 'uniprot_id' column")
                return
            uniprot_ids = batch_df["uniprot_id"].dropna().unique().tolist()
            logger.info("Found %d UniProt IDs to process", len(uniprot_ids))
            
            all_results = []
            for uid in tqdm(uniprot_ids, desc="Processing UniProt IDs"):
                try:
                    df = agg.analyze_uniprot_mutations(uid, write_outputs=not args.no_write)
                    df["batch_uniprot_id"] = uid
                    all_results.append(df)
                    logger.info("Completed %s: %d mutations", uid, len(df))
                except Exception as e:
                    logger.exception("Failed to process %s: %s", uid, e)
            
            # Create combined summary
            if all_results and not args.no_write:
                combined_df = pd.concat(all_results, ignore_index=True)
                summary_path = os.path.join(agg.output_dir, "batch_summary.xlsx")
                combined_df.to_excel(summary_path, index=False)
                logger.info("Wrote batch summary: %s", summary_path)
                print(f"Batch complete: {len(uniprot_ids)} IDs processed, {len(combined_df)} total mutations")
        except Exception as e:
            logger.exception("Batch processing failed: %s", e)
    else:
        # Single ID mode
        df = agg.analyze_uniprot_mutations(args.uniprot_id, write_outputs=not args.no_write)
        print(df.to_string(index=False))


if __name__ == "__main__":
    main()
