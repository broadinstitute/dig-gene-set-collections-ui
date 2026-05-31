from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

import duckdb
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


APP_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = APP_DIR / "data"


def normalize_gene_list(raw: str | Iterable[str]) -> list[str]:
    if isinstance(raw, str):
        tokens = raw.replace(",", "\n").replace("\t", "\n").splitlines()
    else:
        tokens = list(raw)
    out: list[str] = []
    seen: set[str] = set()
    for token in tokens:
        gene = str(token).strip().upper()
        if not gene:
            continue
        if gene in seen:
            continue
        seen.add(gene)
        out.append(gene)
    return out


class RevealRetrievalIndex:
    def __init__(self, data_dir: Path | None = None) -> None:
        self.data_dir = data_dir or DATA_DIR
        self.db = duckdb.connect(str(self.data_dir / "cards.duckdb"), read_only=True)
        self.cards = self.db.execute("select * from cards").df()
        self.signature_genes = self.db.execute("select * from signature_genes").df()

        model = np.load(self.data_dir / "latent_model.npz")
        self.idf = model["idf"]
        self.components = model["components"]
        self.card_embeddings = model["card_embeddings"]
        self.card_embedding_norms = model["card_embedding_norms"]

        metadata = json.loads((self.data_dir / "retrieval_metadata.json").read_text(encoding="utf-8"))
        self.card_order = metadata["cards"]
        self.genes = metadata["genes"]
        self.gene_to_idx = {gene: idx for idx, gene in enumerate(self.genes)}
        self.hypergeom_universe_size = int(metadata["hypergeom_universe_size"])

    def _filter_cards(self, filters: dict[str, object] | None) -> pd.DataFrame:
        if not filters:
            return self.cards.copy()
        df = self.cards.copy()
        for key, value in filters.items():
            if value in (None, "", "All"):
                continue
            df = df[df[key] == value]
        return df

    def _query_embedding(self, genes: list[str]) -> np.ndarray:
        vec = np.zeros(len(self.genes), dtype=float)
        for gene in genes:
            idx = self.gene_to_idx.get(gene)
            if idx is not None:
                vec[idx] = 1.0
        query_tfidf = vec * self.idf
        return query_tfidf @ self.components.T

    def score_query(self, genes: str | Iterable[str], filters: dict[str, object] | None = None, top_k: int = 20) -> pd.DataFrame:
        query_genes = normalize_gene_list(genes)
        if not query_genes:
            return pd.DataFrame()

        candidate_cards = self._filter_cards(filters)
        if candidate_cards.empty:
            return pd.DataFrame()

        query_set = set(query_genes)
        sig = self.signature_genes[self.signature_genes["card_id"].isin(candidate_cards["card_id"])].copy()
        sig["gene_symbol_norm"] = sig["gene_symbol_norm"].astype(str).str.upper()
        sig["is_overlap"] = sig["gene_symbol_norm"].isin(query_set)

        overlap = (
            sig[sig["is_overlap"]]
            .groupby("card_id")
            .agg(
                overlap_score=("weight", "sum"),
                overlap_count=("gene_symbol_norm", "nunique"),
                overlapping_genes=("gene_symbol_norm", lambda x: ",".join(sorted(set(x))[:20])),
            )
            .reset_index()
        )
        result = candidate_cards.merge(overlap, on="card_id", how="left").fillna(
            {"overlap_score": 0.0, "overlap_count": 0, "overlapping_genes": ""}
        )

        card_sizes = sig.groupby("card_id")["gene_symbol_norm"].nunique().rename("card_gene_count").reset_index()
        result = result.merge(card_sizes, on="card_id", how="left")
        q = len(query_set)
        M = self.hypergeom_universe_size
        enrich_scores = []
        for _, row in result.iterrows():
            k = int(row["card_gene_count"])
            x = int(row["overlap_count"])
            pval = float(hypergeom.sf(x - 1, M, k, q)) if x > 0 else 1.0
            enrich_scores.append(-math.log10(max(pval, 1e-300)))
        result["enrichment_score"] = enrich_scores

        query_embedding = self._query_embedding(query_genes)
        query_norm = np.linalg.norm(query_embedding)
        latent_scores = []
        card_idx = {card_id: idx for idx, card_id in enumerate(self.card_order)}
        for _, row in result.iterrows():
            idx = card_idx[row["card_id"]]
            denom = float(self.card_embedding_norms[idx] * query_norm)
            latent_scores.append(float(np.dot(self.card_embeddings[idx], query_embedding) / denom) if denom else 0.0)
        result["latent_score"] = latent_scores

        for column in ("overlap_score", "enrichment_score", "latent_score"):
            values = result[column].to_numpy(dtype=float)
            std = values.std()
            if std == 0:
                result[f"{column}_z"] = 0.0
            else:
                result[f"{column}_z"] = (values - values.mean()) / std

        result["final_score"] = (
            0.45 * result["overlap_score_z"]
            + 0.30 * result["enrichment_score_z"]
            + 0.25 * result["latent_score_z"]
        )
        result = result.sort_values(["final_score", "overlap_score", "enrichment_score"], ascending=False)
        return result.head(top_k).reset_index(drop=True)
