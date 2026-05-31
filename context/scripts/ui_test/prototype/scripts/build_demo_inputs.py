#!/usr/bin/env python3
from __future__ import annotations

import io
import json
import math
from pathlib import Path
from typing import Any
import subprocess
from urllib.parse import urlencode

import pandas as pd
import requests

from common import PROTOTYPE_DIR, ensure_dir, load_catalog, sha256_file, write_json


GTEX_API = "https://gtexportal.org/api/v2"
MYGENE_QUERY = "https://mygene.info/v3/query"


def _url_with_query(base_url: str, params: list[tuple[str, Any]]) -> str:
    return f"{base_url}?{urlencode(params, doseq=True)}"


def _request_json(url: str, params: list[tuple[str, Any]] | None = None, method: str = "GET", data: dict[str, Any] | None = None) -> Any:
    if method == "POST":
        response = requests.post(url, data=data, timeout=120)
    else:
        response = requests.get(url, params=params, timeout=120)
    response.raise_for_status()
    return response.json()


def _paginate_json(base_url: str, params: list[tuple[str, Any]], items_per_page: int) -> list[dict[str, Any]]:
    page = 0
    all_rows: list[dict[str, Any]] = []
    while True:
        page_params = list(params) + [("itemsPerPage", items_per_page), ("page", page)]
        payload = _request_json(base_url, params=page_params)
        all_rows.extend(payload["data"])
        paging = payload["paging_info"]
        if page + 1 >= int(paging["numberOfPages"]):
            break
        page += 1
    return all_rows


def _fetch_top_n_rows(base_url: str, params: list[tuple[str, Any]], items_per_page: int, top_n: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    page = 0
    while len(rows) < top_n:
        page_params = list(params) + [("itemsPerPage", items_per_page), ("page", page)]
        payload = _request_json(base_url, params=page_params)
        batch = payload["data"]
        if not batch:
            break
        rows.extend(batch)
        paging = payload["paging_info"]
        if page + 1 >= int(paging["numberOfPages"]):
            break
        page += 1
    return rows[:top_n]


def build_gtex_inputs(source_cfg: dict[str, Any]) -> None:
    raw_dir = ensure_dir(PROTOTYPE_DIR / "data" / "raw" / "gtex")
    out_root = ensure_dir(PROTOTYPE_DIR / "data" / "intermediate" / "gtex")
    selected = source_cfg["selected_demo_units"]
    top_n = 400
    per_page = 200

    union_rows: list[dict[str, Any]] = []
    for unit in selected:
        tissue_id = unit["tissue_id"]
        print(f"[gtex] fetching top genes for {tissue_id}", flush=True)
        rows = _fetch_top_n_rows(
            f"{GTEX_API}/expression/topExpressedGene",
            [
                ("tissueSiteDetailId", tissue_id),
                ("datasetId", source_cfg["dataset_id"]),
                ("filterMtGene", "true"),
            ],
            per_page,
            top_n,
        )
        cache_path = raw_dir / f"{tissue_id}.top_genes.json"
        write_json(cache_path, {"tissue_id": tissue_id, "rows": rows})
        union_rows.extend(rows)

    union_gene_ids = sorted({row["gencodeId"] for row in union_rows})
    print(f"[gtex] fetching median expression for {len(union_gene_ids)} union genes", flush=True)
    median_rows: list[dict[str, Any]] = []
    chunk_size = 150
    for start in range(0, len(union_gene_ids), chunk_size):
        chunk = union_gene_ids[start : start + chunk_size]
        median_rows.extend(
            _paginate_json(
                f"{GTEX_API}/expression/medianGeneExpression",
                [("datasetId", source_cfg["dataset_id"])] + [("gencodeId", gene_id) for gene_id in chunk],
                5000,
            )
        )
    write_json(raw_dir / "selected_union_median_expression.json", {"n_genes": len(union_gene_ids), "rows": median_rows})

    median_df = pd.DataFrame(median_rows)
    pivot = median_df.pivot_table(
        index=["gencodeId", "geneSymbol"],
        columns="tissueSiteDetailId",
        values="median",
        aggfunc="first",
    ).reset_index()
    tissue_columns = [col for col in pivot.columns if col not in {"gencodeId", "geneSymbol"}]

    for unit in selected:
        card_id = unit["card_id"]
        tissue_id = unit["tissue_id"]
        print(f"[gtex] writing {card_id}", flush=True)
        out_dir = ensure_dir(out_root / card_id)
        df = pivot.copy()
        df["target_median"] = df[tissue_id].fillna(0.0)
        other_cols = [col for col in tissue_columns if col != tissue_id]
        df["other_tissues_mean"] = df[other_cols].fillna(0.0).mean(axis=1)
        df["stat"] = (df["target_median"] + 1.0).map(math.log2) - (df["other_tissues_mean"] + 1.0).map(math.log2)
        df = df.rename(columns={"gencodeId": "gene_id", "geneSymbol": "gene_symbol"})
        deg = (
            df[["gene_id", "gene_symbol", "stat", "target_median", "other_tissues_mean"]]
            .query("stat > 0")
            .sort_values("stat", ascending=False)
            .reset_index(drop=True)
        )
        deg_path = out_dir / "deg.tsv"
        deg.to_csv(deg_path, sep="\t", index=False)

        source_meta = {
            "card_id": card_id,
            "resource_name": "GTEx",
            "dataset_unit_title": source_cfg["candidate_units"][next(i for i, row in enumerate(source_cfg["candidate_units"]) if row["tissue_id"] == tissue_id)]["title"],
            "dataset_unit_type": "tissue_specific_signature",
            "organism": "human",
            "comparison_space_organism": "human",
            "modality": "transcriptomics",
            "tissue_or_system": tissue_id,
            "contrast_label": "tissue-specific expression versus other GTEx tissues",
            "landing_page": f"https://gtexportal.org/home/tissue/{tissue_id}",
            "primary_access_url": _url_with_query(
                f"{GTEX_API}/expression/topExpressedGene",
                [
                    ("tissueSiteDetailId", tissue_id),
                    ("datasetId", source_cfg["dataset_id"]),
                    ("filterMtGene", "true"),
                ],
            ),
            "access_route": "GTEx Portal V2 API topExpressedGene + medianGeneExpression",
            "publication_ids": [],
            "focus_node": tissue_id,
            "source_files": [
                {
                    "path": str(raw_dir / f"{tissue_id}.top_genes.json"),
                    "sha256": sha256_file(raw_dir / f"{tissue_id}.top_genes.json"),
                    "role": "api_cache_top_expressed_gene",
                    "access_route": "GTEx Portal V2 API tissue-specific topExpressedGene query",
                    "landing_page_url": f"https://gtexportal.org/home/tissue/{tissue_id}",
                },
                {
                    "path": str(raw_dir / "selected_union_median_expression.json"),
                    "sha256": sha256_file(raw_dir / "selected_union_median_expression.json"),
                    "role": "api_cache_median_gene_expression",
                    "access_route": "GTEx Portal V2 API medianGeneExpression endpoint (assembled from paginated gene batches)",
                    "landing_page_url": f"https://gtexportal.org/home/tissue/{tissue_id}",
                },
                {
                    "path": str(deg_path),
                    "sha256": sha256_file(deg_path),
                    "role": "comparison_space_deg_tsv",
                },
            ],
        }
        write_json(out_dir / "source_meta.json", source_meta)


def _download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    ensure_dir(path.parent)
    response = requests.get(url, timeout=300)
    response.raise_for_status()
    path.write_bytes(response.content)


def _dump_rda_to_tsv(rda_path: Path, tsv_path: Path) -> None:
    ensure_dir(tsv_path.parent)
    expr = (
        "e <- new.env(); "
        f"nm <- load('{rda_path}', envir=e); "
        "obj <- e[[nm[1]]]; "
        f"utils::write.table(obj, file='{tsv_path}', sep='\\t', row.names=FALSE, quote=FALSE)"
    )
    subprocess.run(
        ["bash", "-lc", f"R_LIBS_USER=ui_test/.Rlib /opt/homebrew/bin/Rscript --vanilla -e \"{expr}\""],
        cwd=str(PROTOTYPE_DIR.parent.parent),
        check=True,
    )


def _query_mygene(ids: list[str], species: str, fields: str) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    chunk_size = 250
    for start in range(0, len(ids), chunk_size):
        chunk = ids[start : start + chunk_size]
        payload = _request_json(
            MYGENE_QUERY,
            method="POST",
            data={
                "q": ",".join(chunk),
                "species": species,
                "fields": fields,
                "size": "1",
            },
        )
        if isinstance(payload, list):
            out.extend(payload)
        elif isinstance(payload, dict) and "hits" in payload:
            out.extend(payload["hits"])
    return out


def build_motrpac_inputs(source_cfg: dict[str, Any]) -> None:
    raw_dir = ensure_dir(PROTOTYPE_DIR / "data" / "raw" / "motrpac")
    out_root = ensure_dir(PROTOTYPE_DIR / "data" / "intermediate" / "motrpac")
    selected = source_cfg["selected_demo_units"]
    unique_objects = sorted({row["object_name"] for row in selected})

    object_tables: dict[str, pd.DataFrame] = {}
    for object_name in unique_objects:
        print(f"[motrpac] preparing {object_name}", flush=True)
        rda_path = raw_dir / object_name
        _download(source_cfg["access_method"]["raw_template"].format(object_name=object_name), rda_path)
        raw_tsv_path = raw_dir / object_name.replace(".rda", ".tsv")
        if not raw_tsv_path.exists():
            _dump_rda_to_tsv(rda_path, raw_tsv_path)
        object_tables[object_name] = pd.read_csv(raw_tsv_path, sep="\t")

    candidate_frames: dict[str, pd.DataFrame] = {}
    all_native_gene_ids: set[str] = set()
    for unit in selected:
        object_name = unit["object_name"]
        card_id = unit["card_id"]
        df = object_tables[object_name].copy()
        df = df[(df["sex"] == unit["sex"]) & (df["comparison_group"] == unit["comparison_group"])].copy()
        df = df.dropna(subset=["feature_ID", "zscore"]).copy()
        df["abs_zscore"] = df["zscore"].abs()
        df = df.sort_values("abs_zscore", ascending=False).head(3000).copy()
        candidate_frames[card_id] = df
        all_native_gene_ids.update(df["feature_ID"].astype(str).tolist())

    rat_hits = _query_mygene(sorted(all_native_gene_ids), "rat", "symbol,ensembl.gene,homologene.genes")
    rat_rows = []
    human_entrez: set[str] = set()
    for hit in rat_hits:
        query = str(hit.get("query", ""))
        homologene = hit.get("homologene", {})
        genes = homologene.get("genes", []) if isinstance(homologene, dict) else []
        human_entrez_id = None
        for row in genes:
            if isinstance(row, list) and len(row) >= 2 and int(row[0]) == 9606:
                human_entrez_id = str(row[1])
                break
        if human_entrez_id:
            human_entrez.add(human_entrez_id)
        rat_rows.append(
            {
                "native_gene_id": query,
                "native_gene_symbol": hit.get("symbol"),
                "human_entrez_id": human_entrez_id,
            }
        )
    rat_map = pd.DataFrame(rat_rows)

    human_hits = _query_mygene(sorted(human_entrez), "human", "symbol,ensembl.gene,entrezgene")
    human_rows = []
    for hit in human_hits:
        human_rows.append(
            {
                "human_entrez_id": str(hit.get("_id") or hit.get("entrezgene") or hit.get("query")),
                "human_symbol": hit.get("symbol"),
                "human_ensembl_gene": (
                    hit.get("ensembl", {}).get("gene")
                    if isinstance(hit.get("ensembl"), dict)
                    else None
                ),
            }
        )
    human_map = pd.DataFrame(human_rows).drop_duplicates(subset=["human_entrez_id"])

    ortholog_map = rat_map.merge(human_map, on="human_entrez_id", how="left")
    ortholog_map.to_csv(raw_dir / "rat_to_human_ortholog_map.tsv", sep="\t", index=False)

    title_by_object = {row["object_name"]: row["title"] for row in source_cfg["candidate_units"]}
    for unit in selected:
        card_id = unit["card_id"]
        print(f"[motrpac] writing {card_id}", flush=True)
        out_dir = ensure_dir(out_root / card_id)
        native_df = candidate_frames[card_id].copy()
        merged = native_df.merge(ortholog_map, left_on="feature_ID", right_on="native_gene_id", how="left")
        mapped = merged.dropna(subset=["human_symbol"]).copy()
        mapped["gene_id"] = mapped["human_ensembl_gene"].fillna(mapped["human_symbol"])
        mapped["gene_symbol"] = mapped["human_symbol"]
        mapped["stat"] = mapped["zscore"]
        mapped["logfc"] = mapped["shrunk_logFC"].fillna(mapped["logFC"])
        mapped["pvalue"] = mapped["p_value"]
        mapped["padj"] = mapped["adj_p_value"]
        deg = (
            mapped[["gene_id", "gene_symbol", "stat", "logfc", "pvalue", "padj"]]
            .dropna(subset=["gene_id", "gene_symbol", "stat"])
            .sort_values("stat", ascending=False)
            .reset_index(drop=True)
        )
        deg_path = out_dir / "deg.tsv"
        deg.to_csv(deg_path, sep="\t", index=False)
        merged.to_csv(out_dir / "native_to_human_map.tsv", sep="\t", index=False)

        source_meta = {
            "card_id": card_id,
            "resource_name": "MoTrPAC",
            "dataset_unit_title": title_by_object[unit["object_name"]],
            "dataset_unit_type": "tissue_assay_contrast_signature",
            "organism": "rat",
            "comparison_space_organism": "human",
            "modality": "transcriptomics",
            "tissue_or_system": unit["tissue_code"],
            "contrast_label": f"{unit['sex']} {unit['comparison_group']} trained versus control",
            "landing_page": "https://motrpac.github.io/MotrpacRatTraining6moData/",
            "primary_access_url": source_cfg["access_method"]["raw_template"].format(object_name=unit["object_name"]),
            "access_route": f"MoTrPAC public data package object {unit['object_name']}",
            "publication_ids": [],
            "focus_node": unit["tissue_code"],
            "ortholog_mapping_method": "mygene.info rat Ensembl to human ortholog via homologene",
            "source_files": [
                {
                    "path": str(raw_dir / unit["object_name"]),
                    "sha256": sha256_file(raw_dir / unit["object_name"]),
                    "role": "public_rda",
                    "obtain_from_url": source_cfg["access_method"]["raw_template"].format(object_name=unit["object_name"]),
                    "access_route": f"Direct raw download for {unit['object_name']}",
                    "landing_page_url": "https://motrpac.github.io/MotrpacRatTraining6moData/",
                },
                {
                    "path": str(raw_dir / unit["object_name"].replace(".rda", ".tsv")),
                    "sha256": sha256_file(raw_dir / unit["object_name"].replace(".rda", ".tsv")),
                    "role": "expanded_public_tsv",
                },
                {
                    "path": str(raw_dir / "rat_to_human_ortholog_map.tsv"),
                    "sha256": sha256_file(raw_dir / "rat_to_human_ortholog_map.tsv"),
                    "role": "ortholog_map",
                    "obtain_from_url": "https://mygene.info/",
                    "access_route": "mygene.info ortholog mapping service queries",
                    "landing_page_url": "https://mygene.info/",
                },
                {
                    "path": str(deg_path),
                    "sha256": sha256_file(deg_path),
                    "role": "comparison_space_deg_tsv",
                },
            ],
            "mapping_stats": {
                "n_native_rows_considered": int(len(native_df)),
                "n_rows_mapped_to_human": int(len(mapped)),
                "n_unique_native_genes": int(native_df["feature_ID"].nunique()),
                "n_unique_human_symbols": int(mapped["human_symbol"].nunique()),
            },
        }
        write_json(out_dir / "source_meta.json", source_meta)


def build_geo_inputs(source_cfg: dict[str, Any]) -> None:
    out_root = ensure_dir(PROTOTYPE_DIR / "data" / "intermediate" / "geo")
    print("[geo] preparing selected Kang IFN-beta comparisons", flush=True)
    source_long = PROTOTYPE_DIR.parent.parent / source_cfg["access_method"]["local_input_long_tsv"]
    df = pd.read_csv(source_long, sep="\t")
    selected_ids = [row["comparison_id"] for row in source_cfg["selected_demo_units"]]
    df = df[df["comparison_id"].isin(selected_ids)].copy()
    df["gene_symbol"] = df["gene_id"]
    out_path = out_root / "GSE96583_ifnb_selected_long.tsv"
    df.to_csv(out_path, sep="\t", index=False)

    for unit in source_cfg["selected_demo_units"]:
        out_dir = ensure_dir(out_root / unit["card_id"])
        comparison_id = unit["comparison_id"]
        comp_df = df[df["comparison_id"] == comparison_id].copy()
        comp_path = out_dir / "deg.tsv"
        comp_df.to_csv(comp_path, sep="\t", index=False)
        source_meta = {
            "card_id": unit["card_id"],
            "resource_name": "GEO",
            "dataset_unit_title": f"GSE96583 {comparison_id.replace('_stim_vs_ctrl', '').replace('_', ' ')} IFN-beta stimulation",
            "dataset_unit_type": "study_contrast_signature",
            "organism": "human",
            "comparison_space_organism": "human",
            "modality": "transcriptomics",
            "tissue_or_system": "PBMC",
            "contrast_label": comparison_id.replace("_", " "),
            "landing_page": source_cfg["access_method"]["series_page"],
            "primary_access_url": source_cfg["access_method"]["supplementary_dir"],
            "access_route": "NCBI GEO supplemental files plus local Kang IFN-beta standardization workflow",
            "publication_ids": [],
            "focus_node": comparison_id,
            "source_files": [
                {
                    "path": str(source_long),
                    "sha256": sha256_file(source_long),
                    "role": "local_long_deg_tsv",
                    "obtain_from_url": source_cfg["access_method"]["supplementary_dir"],
                    "access_route": "NCBI GEO supplementary files for GSE96583, standardized locally into the long DEG table",
                    "landing_page_url": source_cfg["access_method"]["series_page"],
                },
                {
                    "path": str(comp_path),
                    "sha256": sha256_file(comp_path),
                    "role": "comparison_deg_tsv",
                },
            ],
        }
        write_json(out_dir / "source_meta.json", source_meta)


def main() -> None:
    catalog = load_catalog()
    source_map = {row["source_id"]: row for row in catalog["sources"]}
    build_gtex_inputs(source_map["gtex"])
    build_motrpac_inputs(source_map["motrpac"])
    build_geo_inputs(source_map["geo"])


if __name__ == "__main__":
    main()
