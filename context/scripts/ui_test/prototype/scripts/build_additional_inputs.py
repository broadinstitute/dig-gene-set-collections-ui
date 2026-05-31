#!/usr/bin/env python3
from __future__ import annotations

import io
import json
import os
import shutil
import subprocess
import zipfile
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import requests
from scipy import stats
from scipy import sparse

from common import (
    PROTOTYPE_DIR,
    PROTOTYPE_INPUTS_DIR,
    REPO_ROOT,
    ensure_dir,
    sha256_file,
    write_json,
)


MYGENE_QUERY = "https://mygene.info/v3/query"
MOTRPAC_WAT_REPO = "https://github.com/MoTrPAC/MotrpacRatTraining6moWATData"
MOTRPAC_WAT_RAW = "https://raw.githubusercontent.com/MoTrPAC/MotrpacRatTraining6moWATData/master"
FOUR_DN_FILE_JSON = "https://data.4dnucleome.org/files-processed/4DNFIUALWN8X/?format=json"
FOUR_DN_FILE_PAGE = "https://data.4dnucleome.org/files-processed/4DNFIUALWN8X/"
FOUR_DN_OPEN_DATA = "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/8a879fa5-de06-4f08-9c8d-6e33e7bfde18/4DNFIUALWN8X.zip"
FOUR_DN_CHIP_FILE_JSON = "https://data.4dnucleome.org/files-processed/4DNFIGINV1VI/?format=json"
FOUR_DN_CHIP_FILE_PAGE = "https://data.4dnucleome.org/files-processed/4DNFIGINV1VI/"
FOUR_DN_CHIP_DOWNLOAD = "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/f9e48636-43f5-42cb-a6c1-709660d61491/4DNFIGINV1VI.bb"
LOCAL_HG38_GTF = REPO_ROOT / "atac_seq_gene_extractor" / "ref" / "gencode.v47.annotation.gtf.gz"
GENCODE_V47_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz"
UCSC_BIGBEDTOBED_URL = "https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.arm64/bigBedToBed"
JUMP_TARGET_REPO = "https://github.com/jump-cellpainting/JUMP-Target"
JUMP_RESULTS_PAGE = "https://jump-cellpainting.broadinstitute.org/results"
CALR_REPOSITORY_URL = "https://github.com/CalRTeam/cal-repository"
CALR_LANDING_PAGE = "https://github.com/CalRTeam/cal-repository"
GSE42752_PAGE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42752"
HUBMAP_LIVER_PRODUCT_UUID = "6327fd11-974c-48d1-910d-9b981a84c28e"
HUBMAP_LIVER_PAGE = f"https://data-products.hubmapconsortium.org/data_products/{HUBMAP_LIVER_PRODUCT_UUID}/"
HUBMAP_LIVER_PROCESSED_URL = (
    f"https://assets.hubmapconsortium.org/hubmap-data-products/{HUBMAP_LIVER_PRODUCT_UUID}/LV_processed.h5ad"
)


def _download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    ensure_dir(path.parent)
    response = requests.get(url, timeout=300)
    response.raise_for_status()
    path.write_bytes(response.content)


def _request_json(url: str, *, method: str = "GET", data: dict[str, Any] | None = None) -> Any:
    if method == "POST":
        response = requests.post(url, data=data, timeout=120)
    else:
        response = requests.get(url, timeout=120)
    response.raise_for_status()
    return response.json()


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
        cwd=str(REPO_ROOT),
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


def _rat_to_human_map(entrez_ids: list[str]) -> pd.DataFrame:
    rat_hits = _query_mygene(sorted(set(entrez_ids)), "rat", "symbol,ensembl.gene,homologene.genes")
    rat_rows: list[dict[str, Any]] = []
    human_entrez: set[str] = set()
    for hit in rat_hits:
        query = str(hit.get("query", ""))
        genes = hit.get("homologene", {}).get("genes", []) if isinstance(hit.get("homologene"), dict) else []
        human_entrez_id = None
        for row in genes:
            if isinstance(row, list) and len(row) >= 2 and int(row[0]) == 9606:
                human_entrez_id = str(row[1])
                break
        if human_entrez_id:
            human_entrez.add(human_entrez_id)
        rat_rows.append(
            {
                "native_entrez_gene": query,
                "native_gene_symbol": hit.get("symbol"),
                "human_entrez_id": human_entrez_id,
            }
        )

    human_hits = _query_mygene(sorted(human_entrez), "human", "symbol,ensembl.gene,entrezgene")
    human_rows: list[dict[str, Any]] = []
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

    return pd.DataFrame(rat_rows).merge(pd.DataFrame(human_rows).drop_duplicates(subset=["human_entrez_id"]), on="human_entrez_id", how="left")


def _write_readme(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def _write_overlay(input_dir: Path, source_meta: dict[str, Any]) -> Path:
    overlay = {
        "card_id": source_meta["card_id"],
        "source_resource": source_meta["resource_name"],
        "source_dataset_unit": source_meta["dataset_unit_title"],
        "landing_page": source_meta["landing_page"],
        "access_route": source_meta["access_route"],
        "publication_ids": source_meta.get("publication_ids", []),
        "organism": source_meta["organism"],
        "comparison_space_organism": source_meta["comparison_space_organism"],
        "extractor_notes": source_meta.get("extractor_notes", ""),
        "resource_name": source_meta["resource_name"],
        "dataset_unit_title": source_meta["dataset_unit_title"],
        "dataset_unit_type": source_meta["dataset_unit_type"],
        "modality": source_meta["modality"],
        "tissue_or_system": source_meta["tissue_or_system"],
        "contrast_label": source_meta["contrast_label"],
        "focus_node": source_meta.get("focus_node"),
    }
    overlay_path = input_dir / "provenance_overlay.json"
    write_json(overlay_path, overlay)
    return overlay_path


def _preferred_gene_ids(var_names: list[str], hugo_symbols: list[str]) -> np.ndarray:
    gene_ids: list[str] = []
    for var_name, symbol in zip(var_names, hugo_symbols):
        cleaned_symbol = str(symbol).strip()
        if cleaned_symbol and cleaned_symbol.lower() != "nan":
            gene_ids.append(cleaned_symbol.upper())
        else:
            gene_ids.append(str(var_name).split(".")[0].upper())
    return np.asarray(gene_ids, dtype=object)


def _write_sparse_long_counts(matrix: sparse.spmatrix, gene_ids: np.ndarray, barcodes: np.ndarray, out_path: Path) -> None:
    ensure_dir(out_path.parent)
    coo = matrix.tocoo(copy=False)
    chunk_size = 250_000
    for start in range(0, coo.nnz, chunk_size):
        end = min(start + chunk_size, coo.nnz)
        frame = pd.DataFrame(
            {
                "gene_id": gene_ids[coo.col[start:end]],
                "barcode": barcodes[coo.row[start:end]],
                "count": np.asarray(coo.data[start:end], dtype=np.float32),
            }
        )
        frame.to_csv(
            out_path,
            sep="\t",
            index=False,
            header=start == 0,
            mode="w" if start == 0 else "a",
        )


def _stage_link(src: Path, dst: Path) -> Path:
    ensure_dir(dst.parent)
    if dst.exists() or dst.is_symlink():
        return dst
    target = src.resolve()
    try:
        dst.symlink_to(target)
    except OSError:
        shutil.copy2(target, dst)
    return dst


def _ensure_executable(url: str, dst: Path) -> Path:
    if not dst.exists() or dst.stat().st_size == 0:
        _download(url, dst)
        dst.chmod(0o755)
    elif not os.access(dst, os.X_OK):
        dst.chmod(0o755)
    return dst


def _bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    if pvalues.size == 0:
        return pvalues
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    n = float(len(ranked))
    adjusted = np.empty_like(ranked)
    prev = 1.0
    for idx in range(len(ranked) - 1, -1, -1):
        rank = idx + 1.0
        value = min(prev, ranked[idx] * (n / rank))
        adjusted[idx] = value
        prev = value
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return out


def _summarize_gene_by_sample_matrix(path: Path) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, sep="\t", compression="infer", chunksize=512):
        values = chunk.iloc[:, 2:].to_numpy(dtype=np.float64, copy=False)
        log_values = np.log2(values + 1.0)
        n = log_values.shape[1]
        mean = log_values.mean(axis=1)
        if n > 1:
            var = log_values.var(axis=1, ddof=1)
            m2 = var * (n - 1.0)
        else:
            var = np.zeros(log_values.shape[0], dtype=np.float64)
            m2 = np.zeros(log_values.shape[0], dtype=np.float64)
        rows.append(
            pd.DataFrame(
                {
                    "gene_id": chunk.iloc[:, 0].astype(str).to_numpy(),
                    "gene_symbol": chunk.iloc[:, 1].fillna("").astype(str).to_numpy(),
                    "n": np.full(log_values.shape[0], float(n)),
                    "mean": mean,
                    "m2": m2,
                    "var": var,
                }
            )
        )
    summary = pd.concat(rows, ignore_index=True)
    return summary


def _combine_summaries(summary_frames: list[pd.DataFrame]) -> pd.DataFrame:
    if not summary_frames:
        raise ValueError("No summary frames provided")
    combined = summary_frames[0][["gene_id", "gene_symbol", "n", "mean", "m2"]].copy()
    for frame in summary_frames[1:]:
        other = frame[["gene_id", "gene_symbol", "n", "mean", "m2"]].copy()
        if not combined["gene_id"].equals(other["gene_id"]):
            combined = combined.merge(
                other,
                on="gene_id",
                suffixes=("", "_other"),
                how="inner",
            )
            combined["gene_symbol"] = combined["gene_symbol"].where(combined["gene_symbol"] != "", combined["gene_symbol_other"])
        else:
            combined["gene_symbol_other"] = other["gene_symbol"].to_numpy()
            combined["n_other"] = other["n"].to_numpy()
            combined["mean_other"] = other["mean"].to_numpy()
            combined["m2_other"] = other["m2"].to_numpy()
        n1 = combined["n"].to_numpy(dtype=np.float64)
        mean1 = combined["mean"].to_numpy(dtype=np.float64)
        m21 = combined["m2"].to_numpy(dtype=np.float64)
        n2 = combined["n_other"].to_numpy(dtype=np.float64)
        mean2 = combined["mean_other"].to_numpy(dtype=np.float64)
        m22 = combined["m2_other"].to_numpy(dtype=np.float64)
        total_n = n1 + n2
        delta = mean2 - mean1
        combined["mean"] = mean1 + (delta * (n2 / total_n))
        combined["m2"] = m21 + m22 + ((delta * delta) * (n1 * n2 / total_n))
        combined["n"] = total_n
        combined = combined[["gene_id", "gene_symbol", "n", "mean", "m2"]]
    n = combined["n"].to_numpy(dtype=np.float64)
    combined["var"] = np.where(n > 1.0, combined["m2"] / (n - 1.0), 0.0)
    return combined


def build_motrpac_wat_proteomics() -> None:
    slug = "motrpac_wat_proteomics_diff"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")
    staged_dir = ensure_dir(input_dir / "staged")

    rda_path = raw_dir / "PROT_DA.rda"
    flat_tsv = raw_dir / "PROT_DA.tsv"
    _download(f"{MOTRPAC_WAT_RAW}/data/PROT_DA.rda", rda_path)
    if not flat_tsv.exists():
        _dump_rda_to_tsv(rda_path, flat_tsv)

    df = pd.read_csv(flat_tsv, sep="\t")
    prefix = "trained_vs_SED"
    rename_map = {
        f"{prefix}.feature_ID": "feature_id",
        f"{prefix}.gene_symbol": "native_gene_symbol",
        f"{prefix}.entrez_gene": "native_entrez_gene",
        f"{prefix}.logFC": "log2fc",
        f"{prefix}.t": "stat",
        f"{prefix}.P.Value": "pvalue",
        f"{prefix}.adj.P.Val": "padj",
        f"{prefix}.contrast": "contrast",
    }
    selected = df[list(rename_map)].rename(columns=rename_map).dropna(subset=["native_entrez_gene", "log2fc"]).copy()
    selected["native_entrez_gene"] = selected["native_entrez_gene"].astype(str)
    standardized = selected.dropna(subset=["native_gene_symbol"]).copy()
    standardized["gene_id"] = standardized["native_gene_symbol"].astype(str).str.upper()
    standardized["gene_symbol"] = standardized["native_gene_symbol"].astype(str).str.upper()
    standardized_path = staged_dir / "proteomics_diff.tsv"
    standardized[["gene_id", "gene_symbol", "log2fc", "stat", "pvalue", "padj"]].to_csv(standardized_path, sep="\t", index=False)

    contrast_label = str(standardized["contrast"].iloc[0]).replace("_", " ")
    card_id = "MoTrPAC_WAT__Proteomics__F_1W_vs_F_SED"
    source_meta = {
        "card_id": card_id,
        "resource_name": "MoTrPAC",
        "source_resource": "MoTrPAC",
        "source_dataset_unit": "MoTrPAC rat training 6 month WAT proteomics differential analysis",
        "dataset_unit_title": "MoTrPAC WAT proteomics differential abundance: female 1 week trained vs sedentary",
        "dataset_unit_type": "tissue_assay_contrast_signature",
        "organism": "rat",
        "comparison_space_organism": "human",
        "modality": "proteomics",
        "tissue_or_system": "subcutaneous white adipose tissue",
        "contrast_label": contrast_label,
        "landing_page": MOTRPAC_WAT_REPO,
        "primary_access_url": f"{MOTRPAC_WAT_RAW}/data/PROT_DA.rda",
        "access_route": "Public GitHub R package object PROT_DA.rda",
        "publication_ids": ["doi:10.1038/s42255-023-00959-9"],
        "focus_node": "proteomics_diff.tsv",
        "ortholog_mapping_method": "Prototype approximation only: native rat gene symbols were uppercased into a human-comparable symbol space; no explicit ortholog table was applied for this proteomics slice.",
        "extractor_notes": "Prepared from the MoTrPAC WAT public PROT_DA object by selecting the trained_vs_SED contrast block and using normalized shared gene symbols for cross-dataset retrieval.",
        "source_files": [
            {
                "path": str(rda_path),
                "sha256": sha256_file(rda_path),
                "role": "public_rda",
                "obtain_from_url": f"{MOTRPAC_WAT_RAW}/data/PROT_DA.rda",
                "access_route": "Direct raw download of the public WAT proteomics differential object",
                "landing_page_url": MOTRPAC_WAT_REPO,
            },
            {
                "path": str(flat_tsv),
                "sha256": sha256_file(flat_tsv),
                "role": "expanded_public_tsv",
            },
            {
                "path": str(standardized_path),
                "sha256": sha256_file(standardized_path),
                "role": "standardized_proteomics_tsv",
            },
        ],
        "preparation_steps": [
            {
                "source_roles": ["public_rda"],
                "target_role": "expanded_public_tsv",
                "label": "expand public R object",
                "description": "Expand the MoTrPAC WAT PROT_DA R object into a flat TSV for contrast-specific selection.",
                "command": "R_LIBS_USER=ui_test/.Rlib /opt/homebrew/bin/Rscript --vanilla -e \"load(...); write.table(...)\"",
                "script_path": str((Path(__file__).resolve().parent / "build_additional_inputs.py")),
                "working_directory": str(REPO_ROOT),
            },
            {
                "source_roles": ["expanded_public_tsv"],
                "target_role": "standardized_proteomics_tsv",
                "label": "select contrast and normalize gene symbols",
                "description": "Select the female 1 week versus sedentary contrast and normalize rat gene symbols into a retrieval-friendly uppercase symbol space for the proteomics extractor.",
                "command": f"{REPO_ROOT.parent / '.venv' / 'bin' / 'python'} ui_test/prototype/scripts/build_additional_inputs.py",
                "script_path": str((Path(__file__).resolve().parent / "build_additional_inputs.py")),
                "working_directory": str(REPO_ROOT),
            },
        ],
        "mapping_stats": {
            "n_rows_total": int(len(selected)),
            "n_rows_retained": int(len(standardized)),
            "n_unique_native_genes": int(selected["native_gene_symbol"].nunique()),
            "n_unique_normalized_symbols": int(standardized["gene_symbol"].nunique()),
        },
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Downloaded:
- `data/PROT_DA.rda` from `{MOTRPAC_WAT_RAW}/data/PROT_DA.rda`

Prepared:
- expanded the R object to `raw/PROT_DA.tsv`
- selected the `trained_vs_SED` contrast block corresponding to `{contrast_label}`
- normalized native rat gene symbols to uppercase extractor-ready identifiers
- wrote extractor-ready input to `staged/proteomics_diff.tsv`
""",
    )


def build_motrpac_wat_ptm() -> None:
    slug = "motrpac_wat_ptm"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")
    staged_dir = ensure_dir(input_dir / "staged")

    rda_path = raw_dir / "PHOSPHO_DA.rda"
    flat_tsv = raw_dir / "PHOSPHO_DA.tsv"
    _download(f"{MOTRPAC_WAT_RAW}/data/PHOSPHO_DA.rda", rda_path)
    if not flat_tsv.exists():
        _dump_rda_to_tsv(rda_path, flat_tsv)

    df = pd.read_csv(flat_tsv, sep="\t")
    prefix = "trained_vs_SED"
    rename_map = {
        f"{prefix}.feature_ID": "native_site_id",
        f"{prefix}.gene_symbol": "native_gene_symbol",
        f"{prefix}.entrez_gene": "native_entrez_gene",
        f"{prefix}.site": "native_site",
        f"{prefix}.human_feature_ID": "site_id",
        f"{prefix}.human_uniprot": "protein_accession",
        f"{prefix}.human_site": "human_site",
        f"{prefix}.logFC": "log2fc",
        f"{prefix}.t": "stat",
        f"{prefix}.P.Value": "pvalue",
        f"{prefix}.adj.P.Val": "padj",
        f"{prefix}.contrast": "contrast",
    }
    selected = df[list(rename_map)].rename(columns=rename_map).dropna(subset=["native_entrez_gene", "log2fc", "site_id"]).copy()
    selected["native_entrez_gene"] = selected["native_entrez_gene"].astype(str)
    standardized = selected.dropna(subset=["native_gene_symbol", "protein_accession"]).copy()
    standardized["gene_id"] = standardized["native_gene_symbol"].astype(str).str.upper()
    standardized["gene_symbol"] = standardized["native_gene_symbol"].astype(str).str.upper()
    standardized_path = staged_dir / "ptm_site_diff.tsv"
    standardized[["site_id", "gene_id", "gene_symbol", "protein_accession", "log2fc", "stat", "pvalue", "padj"]].to_csv(
        standardized_path, sep="\t", index=False
    )

    contrast_label = str(standardized["contrast"].iloc[0]).replace("_", " ")
    card_id = "MoTrPAC_WAT__Phosphoproteomics__F_1W_vs_F_SED"
    source_meta = {
        "card_id": card_id,
        "resource_name": "MoTrPAC",
        "source_resource": "MoTrPAC",
        "source_dataset_unit": "MoTrPAC rat training 6 month WAT phosphoproteomics differential analysis",
        "dataset_unit_title": "MoTrPAC WAT phosphoproteomics differential abundance: female 1 week trained vs sedentary",
        "dataset_unit_type": "tissue_assay_contrast_signature",
        "organism": "rat",
        "comparison_space_organism": "human",
        "modality": "phosphoproteomics",
        "tissue_or_system": "subcutaneous white adipose tissue",
        "contrast_label": contrast_label,
        "landing_page": MOTRPAC_WAT_REPO,
        "primary_access_url": f"{MOTRPAC_WAT_RAW}/data/PHOSPHO_DA.rda",
        "access_route": "Public GitHub R package object PHOSPHO_DA.rda",
        "publication_ids": ["doi:10.1038/s42255-023-00959-9"],
        "focus_node": "ptm_site_diff.tsv",
        "ortholog_mapping_method": "Prototype approximation only: human phosphosite identifiers from the public MoTrPAC object were retained while native rat gene symbols were normalized to uppercase for gene-level retrieval.",
        "extractor_notes": "Prepared from the MoTrPAC WAT PHOSPHO_DA object by selecting the trained_vs_SED contrast block, retaining rows with human phosphosite identifiers, and normalizing gene symbols for retrieval.",
        "source_files": [
            {
                "path": str(rda_path),
                "sha256": sha256_file(rda_path),
                "role": "public_rda",
                "obtain_from_url": f"{MOTRPAC_WAT_RAW}/data/PHOSPHO_DA.rda",
                "access_route": "Direct raw download of the public WAT phosphoproteomics differential object",
                "landing_page_url": MOTRPAC_WAT_REPO,
            },
            {
                "path": str(flat_tsv),
                "sha256": sha256_file(flat_tsv),
                "role": "expanded_public_tsv",
            },
            {
                "path": str(standardized_path),
                "sha256": sha256_file(standardized_path),
                "role": "standardized_ptm_tsv",
            },
        ],
        "preparation_steps": [
            {
                "source_roles": ["public_rda"],
                "target_role": "expanded_public_tsv",
                "label": "expand public R object",
                "description": "Expand the MoTrPAC WAT PHOSPHO_DA R object into a flat TSV for contrast-specific site selection.",
                "command": "R_LIBS_USER=ui_test/.Rlib /opt/homebrew/bin/Rscript --vanilla -e \"load(...); write.table(...)\"",
                "script_path": str((Path(__file__).resolve().parent / "build_additional_inputs.py")),
                "working_directory": str(REPO_ROOT),
            },
            {
                "source_roles": ["expanded_public_tsv"],
                "target_role": "standardized_ptm_tsv",
                "label": "select contrast and retain human site mappings",
                "description": "Select the female 1 week versus sedentary phosphosite contrast, keep rows with human site identifiers, and normalize gene symbols for retrieval.",
                "command": f"{REPO_ROOT.parent / '.venv' / 'bin' / 'python'} ui_test/prototype/scripts/build_additional_inputs.py",
                "script_path": str((Path(__file__).resolve().parent / "build_additional_inputs.py")),
                "working_directory": str(REPO_ROOT),
            },
        ],
        "mapping_stats": {
            "n_rows_total": int(len(selected)),
            "n_rows_with_human_site_and_gene": int(len(standardized)),
            "n_unique_sites": int(standardized["site_id"].nunique()),
            "n_unique_normalized_symbols": int(standardized["gene_symbol"].nunique()),
        },
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Downloaded:
- `data/PHOSPHO_DA.rda` from `{MOTRPAC_WAT_RAW}/data/PHOSPHO_DA.rda`

Prepared:
- expanded the R object to `raw/PHOSPHO_DA.tsv`
- selected the `trained_vs_SED` contrast block corresponding to `{contrast_label}`
- kept rows with human phosphosite identifiers from the public MoTrPAC object
- normalized native rat gene symbols to uppercase extractor-ready identifiers
- wrote extractor-ready input to `staged/ptm_site_diff.tsv`
""",
    )


def build_4dn_atac() -> None:
    slug = "4dn_atac_bulk_4DNFIUALWN8X_pvalb"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")
    staged_dir = ensure_dir(input_dir / "staged")

    metadata_path = raw_dir / "4DNFIUALWN8X.metadata.json"
    zip_path = raw_dir / "4DNFIUALWN8X.zip"
    if not metadata_path.exists():
        write_json(metadata_path, _request_json(FOUR_DN_FILE_JSON))
    _download(FOUR_DN_OPEN_DATA, zip_path)

    selected_name = "human_individual_cell-type_peaks/PVALB_human_ATAC_peaks_SPM.bed.gz"
    peaks_path = staged_dir / "PVALB_human_ATAC_peaks_SPM.bed.gz"
    if not peaks_path.exists():
        with zipfile.ZipFile(zip_path) as zf:
            with zf.open(selected_name) as src, peaks_path.open("wb") as dst:
                dst.write(src.read())

    card_id = "4DN_4DNFIUALWN8X__ATAC__PVALB"
    source_meta = {
        "card_id": card_id,
        "resource_name": "4DN",
        "source_resource": "4DN",
        "source_dataset_unit": "4DN released processed ATAC peak archive 4DNFIUALWN8X",
        "dataset_unit_title": "4DN released human ATAC peaks 4DNFIUALWN8X: PVALB cell-type peak set",
        "dataset_unit_type": "released_processed_peak_file",
        "organism": "human",
        "comparison_space_organism": "human",
        "modality": "chromatin_accessibility",
        "tissue_or_system": "human individual cell types",
        "contrast_label": "PVALB cell-type accessibility peaks",
        "landing_page": FOUR_DN_FILE_PAGE,
        "primary_access_url": FOUR_DN_OPEN_DATA,
        "access_route": "4DN released processed file open-data download",
        "publication_ids": [],
        "focus_node": "PVALB_human_ATAC_peaks_SPM.bed.gz",
        "extractor_notes": "Prepared by extracting one released cell-type peak BED.GZ from the 4DN processed archive and mapping peaks to genes with the shared ATAC bulk extractor.",
        "source_files": [
            {
                "path": str(metadata_path),
                "sha256": sha256_file(metadata_path),
                "role": "dataset_metadata_json",
                "obtain_from_url": FOUR_DN_FILE_JSON,
                "access_route": "4DN JSON metadata endpoint for the released processed file",
                "landing_page_url": FOUR_DN_FILE_PAGE,
            },
            {
                "path": str(zip_path),
                "sha256": sha256_file(zip_path),
                "role": "public_zip",
                "obtain_from_url": FOUR_DN_OPEN_DATA,
                "access_route": "4DN open-data released archive download",
                "landing_page_url": FOUR_DN_FILE_PAGE,
            },
            {
                "path": str(peaks_path),
                "sha256": sha256_file(peaks_path),
                "role": "extracted_peak_bed",
            },
            {
                "path": str(LOCAL_HG38_GTF),
                "sha256": sha256_file(LOCAL_HG38_GTF),
                "role": "reference_gtf",
                "obtain_from_url": GENCODE_V47_URL,
                "access_route": "GENCODE release 47 human annotation",
                "landing_page_url": "https://www.gencodegenes.org/human/release_47.html",
            },
        ],
        "preparation_steps": [
            {
                "source_roles": ["public_zip"],
                "target_role": "extracted_peak_bed",
                "label": "extract selected peak BED",
                "description": "Extract the PVALB peak BED.GZ from the released 4DN archive.",
                "command": f"{REPO_ROOT.parent / '.venv' / 'bin' / 'python'} ui_test/prototype/scripts/build_additional_inputs.py",
                "script_path": str((Path(__file__).resolve().parent / "build_additional_inputs.py")),
                "working_directory": str(REPO_ROOT),
            },
        ],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Downloaded:
- `raw/4DNFIUALWN8X.metadata.json` from `{FOUR_DN_FILE_JSON}`
- `raw/4DNFIUALWN8X.zip` from `{FOUR_DN_OPEN_DATA}`

Prepared:
- extracted `PVALB_human_ATAC_peaks_SPM.bed.gz` from the released archive to `staged/`
- reused local `gencode.v47.annotation.gtf.gz` for hg38 gene annotation during ATAC extraction
""",
    )


def build_4dn_chipseq() -> None:
    slug = "4dn_chipseq_peak_4DNFIGINV1VI_h3k27ac"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")
    staged_dir = ensure_dir(input_dir / "staged")
    tools_dir = ensure_dir(PROTOTYPE_DIR / "tools")

    metadata_path = raw_dir / "4DNFIGINV1VI.metadata.json"
    bigbed_path = raw_dir / "4DNFIGINV1VI.bb"
    bed_path = staged_dir / "4DNFIGINV1VI.bed"
    converter_path = tools_dir / "bigBedToBed"

    if not metadata_path.exists():
        write_json(metadata_path, _request_json(FOUR_DN_CHIP_FILE_JSON))
    _download(FOUR_DN_CHIP_DOWNLOAD, bigbed_path)
    _ensure_executable(UCSC_BIGBEDTOBED_URL, converter_path)
    if not bed_path.exists():
        subprocess.run([str(converter_path), str(bigbed_path), str(bed_path)], cwd=str(REPO_ROOT), check=True)

    card_id = "4DN_4DNFIGINV1VI__ChIPseq__H3K27ac_HCT116"
    source_meta = {
        "card_id": card_id,
        "resource_name": "4DN",
        "source_resource": "4DN",
        "source_dataset_unit": "4DN released processed ChIP-seq peak file 4DNFIGINV1VI",
        "dataset_unit_title": "4DN released human ChIP-seq peaks 4DNFIGINV1VI: H3K27ac in HCT116 RAD21-AID cells",
        "dataset_unit_type": "released_processed_peak_file",
        "organism": "human",
        "comparison_space_organism": "human",
        "modality": "chip_seq",
        "tissue_or_system": "HCT116 RAD21-AID cells",
        "contrast_label": "H3K27ac peak calls",
        "landing_page": FOUR_DN_CHIP_FILE_PAGE,
        "primary_access_url": FOUR_DN_CHIP_DOWNLOAD,
        "access_route": "4DN direct download of a released processed bigBed peak file",
        "publication_ids": [],
        "focus_node": "4DNFIGINV1VI.bed",
        "extractor_notes": "Prepared by downloading the released 4DN processed H3K27ac bigBed peak file and converting it to plain BED with the UCSC bigBedToBed utility before chipseq_peak extraction.",
        "source_files": [
            {
                "path": str(metadata_path),
                "sha256": sha256_file(metadata_path),
                "role": "dataset_metadata_json",
                "obtain_from_url": FOUR_DN_CHIP_FILE_JSON,
                "access_route": "4DN JSON metadata endpoint for the released processed ChIP-seq file",
                "landing_page_url": FOUR_DN_CHIP_FILE_PAGE,
            },
            {
                "path": str(bigbed_path),
                "sha256": sha256_file(bigbed_path),
                "role": "public_bigbed",
                "obtain_from_url": FOUR_DN_CHIP_DOWNLOAD,
                "access_route": "4DN direct download for the released processed bigBed peak file",
                "landing_page_url": FOUR_DN_CHIP_FILE_PAGE,
            },
            {
                "path": str(bed_path),
                "sha256": sha256_file(bed_path),
                "role": "converted_peak_bed",
            },
            {
                "path": str(LOCAL_HG38_GTF),
                "sha256": sha256_file(LOCAL_HG38_GTF),
                "role": "reference_gtf",
                "obtain_from_url": GENCODE_V47_URL,
                "access_route": "GENCODE release 47 human annotation",
                "landing_page_url": "https://www.gencodegenes.org/human/release_47.html",
            },
        ],
        "preparation_steps": [
            {
                "source_roles": ["public_bigbed"],
                "target_role": "converted_peak_bed",
                "label": "convert bigBed peaks to BED",
                "description": "Convert the released 4DN ChIP-seq bigBed peak file into plain BED so it can be consumed by the shared chipseq_peak extractor.",
                "command": f"{converter_path} {bigbed_path} {bed_path}",
                "script_path": str(converter_path),
                "working_directory": str(REPO_ROOT),
            },
        ],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Downloaded:
- `raw/4DNFIGINV1VI.metadata.json` from `{FOUR_DN_CHIP_FILE_JSON}`
- `raw/4DNFIGINV1VI.bb` from `{FOUR_DN_CHIP_DOWNLOAD}`

Prepared:
- downloaded `ui_test/prototype/tools/bigBedToBed` from `{UCSC_BIGBEDTOBED_URL}`
- converted the released bigBed peak file to `staged/4DNFIGINV1VI.bed`
- reused local `gencode.v47.annotation.gtf.gz` for hg38 gene annotation during ChIP-seq extraction
""",
    )


def build_gtex_grouped_bulk_rna() -> None:
    slug = "gtex_bulk_rna_grouped"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")
    staged_dir = ensure_dir(input_dir / "staged")

    tissue_specs = [
        ("Muscle - Skeletal", "Muscle.tsv.gz", "muscle_skeletal"),
        ("Liver", "Liver.tsv.gz", "liver"),
        ("Whole Blood", "Blood.tsv.gz", "whole_blood"),
        ("Lung", "Lung.tsv.gz", "lung"),
        ("Pancreas", "Pancreas.tsv.gz", "pancreas"),
    ]
    tissue_source_dir = REPO_ROOT / "rna_seq_gene_extractor" / "results" / "tmp" / "gtex_split_by_tissue"
    sample_meta_src = REPO_ROOT / "rna_seq_gene_extractor" / "prep" / "gtex_sample_metadata.tsv"
    sample_meta_path = _stage_link(sample_meta_src, raw_dir / "gtex_sample_metadata.tsv")

    staged_tissue_files: list[tuple[str, str, Path]] = []
    summary_by_slug: dict[str, pd.DataFrame] = {}
    membership_rows: list[dict[str, str]] = []

    for tissue_label, filename, tissue_slug in tissue_specs:
        src = tissue_source_dir / filename
        dst = _stage_link(src, raw_dir / filename)
        staged_tissue_files.append((tissue_label, tissue_slug, dst))
        summary_by_slug[tissue_slug] = _summarize_gene_by_sample_matrix(dst)

        header = pd.read_csv(dst, sep="\t", compression="infer", nrows=0)
        for sample_id in header.columns[2:]:
            membership_rows.append({"sample_id": str(sample_id), "tissue": tissue_label})

    selected_sample_meta = staged_dir / "selected_sample_metadata.tsv"
    pd.DataFrame(membership_rows).to_csv(selected_sample_meta, sep="\t", index=False)

    deg_long_rows: list[pd.DataFrame] = []
    for tissue_label, tissue_slug, _path in staged_tissue_files:
        focal = summary_by_slug[tissue_slug].copy()
        rest = _combine_summaries([summary_by_slug[slug] for _label, slug, _dst in staged_tissue_files if slug != tissue_slug])
        if not focal["gene_id"].equals(rest["gene_id"]):
            merged = focal.merge(
                rest[["gene_id", "n", "mean", "var"]],
                on="gene_id",
                suffixes=("", "_rest"),
                how="inner",
            )
        else:
            merged = focal.copy()
            merged["n_rest"] = rest["n"].to_numpy()
            merged["mean_rest"] = rest["mean"].to_numpy()
            merged["var_rest"] = rest["var"].to_numpy()

        mean_focal = merged["mean"].to_numpy(dtype=np.float64)
        mean_rest = merged["mean_rest"].to_numpy(dtype=np.float64)
        var_focal = np.clip(merged["var"].to_numpy(dtype=np.float64), 0.0, None)
        var_rest = np.clip(merged["var_rest"].to_numpy(dtype=np.float64), 0.0, None)
        n_focal = merged["n"].to_numpy(dtype=np.float64)
        n_rest = merged["n_rest"].to_numpy(dtype=np.float64)
        se2 = (var_focal / n_focal) + (var_rest / n_rest)
        se = np.sqrt(np.clip(se2, 1e-12, None))
        stat = (mean_focal - mean_rest) / se
        dof_num = se2 * se2
        dof_den = ((var_focal / n_focal) ** 2) / np.clip(n_focal - 1.0, 1.0, None) + ((var_rest / n_rest) ** 2) / np.clip(n_rest - 1.0, 1.0, None)
        dof = np.maximum(n_focal + n_rest - 2.0, 1.0)
        valid_dof = dof_den > 0.0
        dof[valid_dof] = dof_num[valid_dof] / dof_den[valid_dof]
        pvalue = 2.0 * stats.t.sf(np.abs(stat), df=np.clip(dof, 1.0, None))
        pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
        padj = _bh_adjust(pvalue.astype(np.float64))
        comparison_id = f"{tissue_slug}_vs_rest"
        deg_long_rows.append(
            pd.DataFrame(
                {
                    "comparison_id": comparison_id,
                    "gene_id": merged["gene_id"].astype(str),
                    "gene_symbol": merged["gene_symbol"].astype(str),
                    "log2fc": mean_focal - mean_rest,
                    "stat": stat,
                    "pvalue": pvalue,
                    "padj": padj,
                }
            )
        )

    deg_long = pd.concat(deg_long_rows, ignore_index=True)
    deg_long_path = staged_dir / "gtex_tissue_vs_rest_deg_long.tsv"
    deg_long.to_csv(deg_long_path, sep="\t", index=False)

    card_id = "GTEx_OpenAccess__BulkRNA__Grouped"
    source_files: list[dict[str, Any]] = [
        {
            "path": str(sample_meta_path),
            "sha256": sha256_file(sample_meta_path),
            "role": "gtex_sample_annotations_tsv",
            "access_route": "Downloaded from the GTEx Portal open-access datasets page",
            "landing_page_url": "https://www.gtexportal.org/home/datasets",
        },
        {
            "path": str(selected_sample_meta),
            "sha256": sha256_file(selected_sample_meta),
            "role": "gtex_selected_sample_metadata_tsv",
        },
        {
            "path": str(deg_long_path),
            "sha256": sha256_file(deg_long_path),
            "role": "comparison_space_deg_tsv",
        },
    ]
    for tissue_label, tissue_slug, path in staged_tissue_files:
        source_files.append(
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "role": f"gtex_counts_{tissue_slug}",
                "access_route": f"Downloaded GTEx open-access bulk RNA counts and subset to {tissue_label} samples",
                "landing_page_url": f"https://gtexportal.org/home/tissue/{tissue_slug}",
            }
        )

    prep_source_roles = [f"gtex_counts_{slug_name}" for _label, slug_name, _path in staged_tissue_files]
    source_meta = {
        "card_id": card_id,
        "resource_name": "GTEx",
        "source_resource": "GTEx",
        "source_dataset_unit": "GTEx open-access bulk tissue expression",
        "dataset_unit_title": "GTEx open-access grouped bulk RNA tissue-versus-rest contrasts",
        "dataset_unit_type": "grouped_bulk_tissue_expression",
        "organism": "human",
        "comparison_space_organism": "human",
        "modality": "bulk_rna",
        "tissue_or_system": "multiple canonical tissues",
        "contrast_label": "five tissue-versus-rest contrasts",
        "landing_page": "https://www.gtexportal.org/home/datasets",
        "primary_access_url": "https://www.gtexportal.org/home/datasets",
        "access_route": "GTEx open-access portal files already staged locally in the workspace",
        "publication_ids": [],
        "focus_node": "gtex_tissue_vs_rest_deg_long.tsv",
        "extractor_notes": "Built a grouped long-format DEG table for five canonical tissues by comparing each selected tissue against the pooled remainder of the selected tissue set using Welch t statistics on log2(count+1) expression.",
        "source_files": source_files,
        "preparation_steps": [
            {
                "source_roles": [*prep_source_roles, "gtex_sample_annotations_tsv"],
                "target_role": "gtex_selected_sample_metadata_tsv",
                "label": "assemble selected GTEx tissue metadata",
                "description": "Collect the sample identifiers used in the five selected GTEx tissue splits and write a compact sample metadata table for the grouped prototype run.",
                "command": f"{REPO_ROOT.parent / '.venv' / 'bin' / 'python'} ui_test/prototype/scripts/build_additional_inputs.py",
                "script_path": str((Path(__file__).resolve().parent / "build_additional_inputs.py")),
                "working_directory": str(REPO_ROOT),
            },
            {
                "source_roles": [*prep_source_roles, "gtex_selected_sample_metadata_tsv"],
                "target_role": "comparison_space_deg_tsv",
                "label": "compute tissue-versus-rest differential statistics",
                "description": "Compute a grouped long-format DEG table for five GTEx tissues by comparing each focal tissue against the pooled remainder of the selected tissue set.",
                "command": f"{REPO_ROOT.parent / '.venv' / 'bin' / 'python'} ui_test/prototype/scripts/build_additional_inputs.py",
                "script_path": str((Path(__file__).resolve().parent / "build_additional_inputs.py")),
                "working_directory": str(REPO_ROOT),
            },
        ],
        "mapping_stats": {
            "n_selected_tissues": len(staged_tissue_files),
            "n_selected_samples": len(membership_rows),
            "n_deg_rows": int(len(deg_long)),
        },
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    tissue_bullets = "\n".join(f"- `{path.name}` for `{label}`" for label, _slug, path in staged_tissue_files)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Staged from local GTEx open-access files:
{tissue_bullets}
- `gtex_sample_metadata.tsv`

Prepared:
- wrote `staged/selected_sample_metadata.tsv` for the five selected tissues
- computed Welch t-test tissue-versus-rest statistics on log2(count+1) expression
- wrote grouped extractor-ready long DE input to `staged/gtex_tissue_vs_rest_deg_long.tsv`
""",
    )


def build_tcga_brca_splicing() -> None:
    slug = "tcga_brca_splice_event_matrix"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")
    staged_dir = ensure_dir(input_dir / "staged")

    raw_matrix_src = REPO_ROOT / "splice_extractor" / "data" / "public" / "tcga_spliceseq" / "BRCA" / "BRCA_spliceseq_subset5000_norm.tsv"
    raw_annotations_src = REPO_ROOT / "splice_extractor" / "data" / "public" / "tcga_spliceseq" / "BRCA" / "BRCA_sample_annotations.tsv"
    prepared_dir = REPO_ROOT / "splice_extractor" / "results" / "runtime" / "splice_prepare_public_BRCA_subset5000_norm"

    raw_matrix = _stage_link(raw_matrix_src, raw_dir / raw_matrix_src.name)
    raw_annotations = _stage_link(raw_annotations_src, raw_dir / raw_annotations_src.name)
    psi_matrix = _stage_link(prepared_dir / "psi_matrix.tsv", staged_dir / "psi_matrix.tsv")
    sample_meta = _stage_link(prepared_dir / "sample_metadata.tsv", staged_dir / "sample_metadata.tsv")
    event_meta = _stage_link(prepared_dir / "event_metadata.tsv", staged_dir / "event_metadata.tsv")

    card_id = "TCGA_BRCA__Splicing__Tumor_vs_Normal"
    source_meta = {
        "card_id": card_id,
        "resource_name": "TCGA SpliceSeq",
        "source_resource": "TCGA SpliceSeq",
        "source_dataset_unit": "TCGA BRCA SpliceSeq public subset5000 normalized matrix",
        "dataset_unit_title": "TCGA BRCA alternative splicing: tumor versus adjacent normal",
        "dataset_unit_type": "cohort_case_control_splicing_signature",
        "organism": "human",
        "comparison_space_organism": "human",
        "modality": "splicing",
        "tissue_or_system": "breast tumor",
        "contrast_label": "tumor versus adjacent normal",
        "landing_page": "https://bioinformatics.mdanderson.org/TCGASpliceSeq/",
        "primary_access_url": "https://bioinformatics.mdanderson.org/TCGASpliceSeq/",
        "access_route": "Public TCGA SpliceSeq BRCA download with local subset5000 normalization staging",
        "publication_ids": [],
        "focus_node": "psi_matrix.tsv",
        "extractor_notes": "Uses the prepared public TCGA BRCA subset5000 normalized PSI matrix and metadata to score tumor-versus-normal alternative splicing programs with the shared splice-event matrix extractor.",
        "source_files": [
            {
                "path": str(raw_matrix),
                "sha256": sha256_file(raw_matrix),
                "role": "tcga_spliceseq_subset_matrix_tsv",
                "obtain_from_url": "https://bioinformatics.mdanderson.org/TCGASpliceSeq/",
                "access_route": "Public TCGA SpliceSeq BRCA matrix download, then local subset5000 normalization staging",
                "landing_page_url": "https://bioinformatics.mdanderson.org/TCGASpliceSeq/",
            },
            {
                "path": str(raw_annotations),
                "sha256": sha256_file(raw_annotations),
                "role": "tcga_spliceseq_sample_annotations_tsv",
                "obtain_from_url": "https://bioinformatics.mdanderson.org/TCGASpliceSeq/",
                "access_route": "Public TCGA SpliceSeq BRCA sample annotation download",
                "landing_page_url": "https://bioinformatics.mdanderson.org/TCGASpliceSeq/",
            },
            {
                "path": str(psi_matrix),
                "sha256": sha256_file(psi_matrix),
                "role": "splicing_prepared_psi_matrix_tsv",
            },
            {
                "path": str(sample_meta),
                "sha256": sha256_file(sample_meta),
                "role": "splicing_prepared_sample_metadata_tsv",
            },
            {
                "path": str(event_meta),
                "sha256": sha256_file(event_meta),
                "role": "splicing_prepared_event_metadata_tsv",
            },
        ],
        "preparation_steps": [
            {
                "source_roles": ["tcga_spliceseq_subset_matrix_tsv", "tcga_spliceseq_sample_annotations_tsv"],
                "target_role": "splicing_prepared_psi_matrix_tsv",
                "label": "prepare public SpliceSeq matrix",
                "description": "Normalize the public TCGA BRCA SpliceSeq subset and emit the extractor-ready PSI matrix.",
                "command": "geneset-extractors workflows splice_prepare_public --input_mode tcga_spliceseq ...",
                "script_path": str(REPO_ROOT / "splice_extractor" / "notes" / "alternative_splicing_analysis_agent_prompt.md"),
                "working_directory": str(REPO_ROOT),
            },
            {
                "source_roles": ["tcga_spliceseq_subset_matrix_tsv", "tcga_spliceseq_sample_annotations_tsv"],
                "target_role": "splicing_prepared_sample_metadata_tsv",
                "label": "prepare sample metadata",
                "description": "Emit extractor-ready sample metadata from the public TCGA BRCA sample annotations.",
                "command": "geneset-extractors workflows splice_prepare_public --input_mode tcga_spliceseq ...",
                "script_path": str(REPO_ROOT / "splice_extractor" / "notes" / "alternative_splicing_analysis_agent_prompt.md"),
                "working_directory": str(REPO_ROOT),
            },
            {
                "source_roles": ["tcga_spliceseq_subset_matrix_tsv"],
                "target_role": "splicing_prepared_event_metadata_tsv",
                "label": "prepare event metadata",
                "description": "Emit extractor-ready event metadata from the public TCGA BRCA SpliceSeq subset.",
                "command": "geneset-extractors workflows splice_prepare_public --input_mode tcga_spliceseq ...",
                "script_path": str(REPO_ROOT / "splice_extractor" / "notes" / "alternative_splicing_analysis_agent_prompt.md"),
                "working_directory": str(REPO_ROOT),
            },
        ],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Staged from local public-prepared BRCA splicing assets:
- `raw/{raw_matrix.name}`
- `raw/{raw_annotations.name}`
- `staged/psi_matrix.tsv`
- `staged/sample_metadata.tsv`
- `staged/event_metadata.tsv`

Prepared inputs were reused from the existing `splice_prepare_public_BRCA_subset5000_norm` workflow output to avoid duplicating the full BRCA public PSI matrix.
""",
    )


def build_prism_quick_drug_response() -> None:
    slug = "depmap_prism_quick_grouped"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")

    prism_dir = REPO_ROOT / "drug_response_extractor" / "retest_drug_response_with_bundle" / "work" / "prism_quick"
    response_path = _stage_link(prism_dir / "response_long.tsv", raw_dir / "response_long.tsv")
    groups_path = _stage_link(prism_dir / "groups.tsv", raw_dir / "groups.tsv")
    targets_path = _stage_link(prism_dir / "drug_targets.tsv", raw_dir / "drug_targets.tsv")

    card_id = "DepMap_PRISM__Grouped"
    source_meta = {
        "card_id": card_id,
        "resource_name": "DepMap PRISM",
        "source_resource": "DepMap PRISM",
        "source_dataset_unit": "DepMap PRISM quick public subset grouped by primary tissue",
        "dataset_unit_title": "DepMap PRISM grouped drug-response programs by primary tissue",
        "dataset_unit_type": "grouped_screen_programs",
        "organism": "human",
        "comparison_space_organism": "human",
        "modality": "drug_response",
        "tissue_or_system": "multiple cancer primary tissues",
        "contrast_label": "group versus rest",
        "landing_page": "https://depmap.org/portal/prism/",
        "primary_access_url": "https://depmap.org/portal/prism/",
        "access_route": "Public PRISM data staged locally into response, group, and target tables",
        "publication_ids": [],
        "focus_node": "response_long.tsv",
        "extractor_notes": "Uses the public PRISM quick subset staged as response, group, and target tables to emit grouped tissue-versus-rest drug-response gene programs.",
        "source_files": [
            {
                "path": str(response_path),
                "sha256": sha256_file(response_path),
                "role": "prism_response_long_tsv",
                "obtain_from_url": "https://depmap.org/portal/prism/",
                "access_route": "Public PRISM screen download staged into long-format response values",
                "landing_page_url": "https://depmap.org/portal/prism/",
            },
            {
                "path": str(groups_path),
                "sha256": sha256_file(groups_path),
                "role": "prism_groups_tsv",
                "obtain_from_url": "https://depmap.org/portal/prism/",
                "access_route": "Public PRISM cell-line metadata staged into primary-tissue group labels",
                "landing_page_url": "https://depmap.org/portal/prism/",
            },
            {
                "path": str(targets_path),
                "sha256": sha256_file(targets_path),
                "role": "prism_drug_targets_tsv",
                "obtain_from_url": "https://depmap.org/portal/prism/",
                "access_route": "Public PRISM treatment annotations staged into drug-target edges",
                "landing_page_url": "https://depmap.org/portal/prism/",
            },
        ],
        "preparation_steps": [],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Staged from local public-prepared PRISM assets:
- `raw/response_long.tsv`
- `raw/groups.tsv`
- `raw/drug_targets.tsv`

These files were reused from the existing `retest_drug_response_with_bundle/work/prism_quick` preparation output because they already match the shared `drug_response_screen` input contract.
""",
    )


def build_tcga_brca_cnv() -> None:
    slug = "tcga_brca_cnv_segments"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")

    seg_src = REPO_ROOT / "cnv_gene_extractor" / "cnv_validation_v3" / "work" / "brca.subset50.seg"
    seg_path = _stage_link(seg_src, raw_dir / "brca.subset50.seg")
    gtf_path = _stage_link(LOCAL_HG38_GTF, raw_dir / LOCAL_HG38_GTF.name)

    card_id = "TCGA_BRCA__CNV__Grouped"
    source_meta = {
        "card_id": card_id,
        "resource_name": "TCGA / GDC",
        "source_resource": "TCGA / GDC",
        "source_dataset_unit": "TCGA BRCA copy-number segment subset",
        "dataset_unit_title": "TCGA BRCA copy-number segments grouped sample programs",
        "dataset_unit_type": "grouped_sample_programs",
        "organism": "human",
        "comparison_space_organism": "human",
        "modality": "copy_number",
        "tissue_or_system": "breast tumor",
        "contrast_label": "sample-specific amplification and deletion programs",
        "landing_page": "https://portal.gdc.cancer.gov/projects/TCGA-BRCA",
        "primary_access_url": "https://portal.gdc.cancer.gov/projects/TCGA-BRCA",
        "access_route": "Public GDC BRCA copy-number segments staged locally as a 50-sample subset",
        "publication_ids": [],
        "focus_node": "brca.subset50.seg",
        "extractor_notes": "Uses a local 50-sample subset of public TCGA BRCA segments to emit per-sample amplification and deletion programs with the shared CNV extractor.",
        "source_files": [
            {
                "path": str(seg_path),
                "sha256": sha256_file(seg_path),
                "role": "tcga_brca_segment_subset_tsv",
                "obtain_from_url": "https://portal.gdc.cancer.gov/projects/TCGA-BRCA",
                "access_route": "Public GDC segment download staged into a 50-sample subset",
                "landing_page_url": "https://portal.gdc.cancer.gov/projects/TCGA-BRCA",
            },
            {
                "path": str(gtf_path),
                "sha256": sha256_file(gtf_path),
                "role": "reference_gtf",
                "obtain_from_url": GENCODE_V47_URL,
                "access_route": "GENCODE release 47 human annotation",
                "landing_page_url": "https://www.gencodegenes.org/human/release_47.html",
            },
        ],
        "preparation_steps": [],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Staged from local public-prepared BRCA CNV assets:
- `raw/brca.subset50.seg`
- `raw/{gtf_path.name}`

The segment subset was reused from the existing BRCA CNV validation workspace to avoid duplicating the full project download.
""",
    )


def build_jump_cellpainting_morphology() -> None:
    slug = "jump_cellpainting_u2os_morphology"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")
    staged_dir = ensure_dir(input_dir / "staged")
    bundle_dir = ensure_dir(input_dir / "bundle")

    prepared_root = REPO_ROOT / "morphology_extractor" / "morphology_validation"
    query_root = REPO_ROOT / "morphology_extractor" / "morphology_validation_retest" / "work"
    bundle_root = REPO_ROOT / "morphology_extractor" / "morphology_validation_round2" / "work" / "jump_bundle_u2os_default_dcfe8e6"

    plate_a = _stage_link(prepared_root / "data" / "prepared" / "BR00117012.profiles.tsv", raw_dir / "BR00117012.profiles.tsv")
    plate_b = _stage_link(prepared_root / "data" / "prepared" / "BR00117022.profiles.tsv", raw_dir / "BR00117022.profiles.tsv")
    exp_meta = _stage_link(prepared_root / "data" / "prepared" / "experimental_metadata.tsv", raw_dir / "experimental_metadata.tsv")
    query_candidates = _stage_link(prepared_root / "data" / "prepared" / "query_candidates.json", raw_dir / "query_candidates.json")
    query_profiles = _stage_link(query_root / "query_profiles.tsv", staged_dir / "query_profiles.tsv")
    query_metadata = _stage_link(query_root / "query_metadata.tsv", staged_dir / "query_metadata.tsv")

    bundle_files: list[tuple[str, Path]] = []
    for name in (
        "morphology_jump_target_pilot_U2OS_48_v1.bundle.json",
        "reference_profiles.tsv.gz",
        "reference_metadata.tsv.gz",
        "compound_targets.tsv.gz",
        "feature_schema.tsv.gz",
        "feature_stats.tsv.gz",
        "target_annotations.tsv.gz",
        "bundle_summary.json",
    ):
        src = bundle_root / name
        if src.exists():
            bundle_files.append((name, _stage_link(src, bundle_dir / name)))

    card_id = "JUMP_CellPainting__Morphology__Grouped"
    source_files: list[dict[str, Any]] = [
        {
            "path": str(plate_a),
            "sha256": sha256_file(plate_a),
            "role": "jump_prepared_plate_profiles_tsv",
            "obtain_from_url": JUMP_RESULTS_PAGE,
            "access_route": "Public JUMP Cell Painting U2OS plate-level morphology profiles staged locally from the consortium release",
            "landing_page_url": JUMP_RESULTS_PAGE,
        },
        {
            "path": str(plate_b),
            "sha256": sha256_file(plate_b),
            "role": "jump_prepared_plate_profiles_tsv",
            "obtain_from_url": JUMP_RESULTS_PAGE,
            "access_route": "Public JUMP Cell Painting U2OS plate-level morphology profiles staged locally from the consortium release",
            "landing_page_url": JUMP_RESULTS_PAGE,
        },
        {
            "path": str(exp_meta),
            "sha256": sha256_file(exp_meta),
            "role": "jump_experimental_metadata_tsv",
            "obtain_from_url": JUMP_TARGET_REPO,
            "access_route": "Public JUMP-Target metadata and platemap-derived experimental metadata staged locally",
            "landing_page_url": JUMP_TARGET_REPO,
        },
        {
            "path": str(query_candidates),
            "sha256": sha256_file(query_candidates),
            "role": "jump_query_candidates_json",
            "obtain_from_url": JUMP_TARGET_REPO,
            "access_route": "Locally generated list of public JUMP query candidates derived from the released target pilot metadata",
            "landing_page_url": JUMP_TARGET_REPO,
        },
        {
            "path": str(query_profiles),
            "sha256": sha256_file(query_profiles),
            "role": "morphology_query_profiles_tsv",
        },
        {
            "path": str(query_metadata),
            "sha256": sha256_file(query_metadata),
            "role": "morphology_query_metadata_tsv",
        },
    ]
    for name, path in bundle_files:
        role = {
            "reference_profiles.tsv.gz": "morphology_reference_profiles_tsv",
            "reference_metadata.tsv.gz": "morphology_reference_metadata_tsv",
            "compound_targets.tsv.gz": "morphology_compound_targets_tsv",
            "feature_schema.tsv.gz": "morphology_feature_schema_tsv",
            "feature_stats.tsv.gz": "morphology_feature_stats_tsv",
            "target_annotations.tsv.gz": "morphology_target_annotations_tsv",
            "morphology_jump_target_pilot_U2OS_48_v1.bundle.json": "morphology_reference_bundle_manifest",
            "bundle_summary.json": "morphology_bundle_summary_json",
        }.get(name)
        if role:
            source_files.append(
                {
                    "path": str(path),
                    "sha256": sha256_file(path),
                    "role": role,
                    "obtain_from_url": JUMP_RESULTS_PAGE,
                    "access_route": "Locally prepared reference-bundle component derived from public JUMP Cell Painting releases",
                    "landing_page_url": JUMP_RESULTS_PAGE,
                }
            )

    source_meta = {
        "card_id": card_id,
        "resource_name": "JUMP Cell Painting",
        "source_resource": "JUMP Cell Painting",
        "source_dataset_unit": "JUMP Cell Painting U2OS 48h public target pilot morphology queries",
        "dataset_unit_title": "JUMP Cell Painting U2OS 48h morphology target-pilot queries",
        "dataset_unit_type": "grouped_query_against_reference_bundle",
        "organism": "human",
        "comparison_space_organism": "human",
        "modality": "morphology",
        "tissue_or_system": "U2OS cell line",
        "contrast_label": "public held-out perturbation queries against reference bundle",
        "landing_page": JUMP_RESULTS_PAGE,
        "primary_access_url": JUMP_TARGET_REPO,
        "access_route": "Public JUMP Cell Painting target pilot inputs and reference bundle staged locally for morphology retrieval",
        "publication_ids": [],
        "focus_node": "query_profiles.tsv",
        "extractor_notes": "Uses a small public U2OS 48h held-out query panel staged from JUMP Cell Painting target-pilot data and scores each query against the existing public reference bundle in mechanism mode with self-matches excluded.",
        "source_files": source_files,
        "preparation_steps": [
            {
                "source_roles": [
                    "jump_prepared_plate_profiles_tsv",
                    "jump_experimental_metadata_tsv",
                    "jump_query_candidates_json",
                ],
                "target_role": "morphology_query_profiles_tsv",
                "label": "assemble held-out JUMP query profiles",
                "description": "Select a small public U2OS 48h target-pilot query panel and emit the feature table used by the morphology extractor.",
                "command": "precomputed in morphology_validation_retest/work/query_profiles.tsv",
                "script_path": str(query_root / "query_profiles.tsv"),
                "working_directory": str(REPO_ROOT),
            },
            {
                "source_roles": [
                    "jump_experimental_metadata_tsv",
                    "jump_query_candidates_json",
                ],
                "target_role": "morphology_query_metadata_tsv",
                "label": "assemble held-out JUMP query metadata",
                "description": "Emit query metadata describing perturbation IDs, perturbation types, grouping, and display labels for the public held-out panel.",
                "command": "precomputed in morphology_validation_retest/work/query_metadata.tsv",
                "script_path": str(query_root / "query_metadata.tsv"),
                "working_directory": str(REPO_ROOT),
            },
        ],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Staged from local public-prepared JUMP Cell Painting assets:
- `raw/BR00117012.profiles.tsv`
- `raw/BR00117022.profiles.tsv`
- `raw/experimental_metadata.tsv`
- `raw/query_candidates.json`
- `staged/query_profiles.tsv`
- `staged/query_metadata.tsv`
- `bundle/*` from the public U2OS 48h morphology reference bundle

The staged query panel reuses the held-out public U2OS 48h perturbations from the existing morphology validation workspace and reruns them through the shared extractor with current provenance enabled.
""",
    )


def build_calr_fig2_ontology() -> None:
    slug = "calr_fig2_ontology_grouped"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")

    calr_data = _stage_link(REPO_ROOT / "calorimetry_extractor" / "work" / "data" / "cal_repository" / "Fig2_CalR_data.csv", raw_dir / "Fig2_CalR_data.csv")
    session_csv = _stage_link(REPO_ROOT / "calorimetry_extractor" / "work" / "data" / "cal_repository" / "Fig2_CalR_Davis_Session.csv", raw_dir / "Fig2_CalR_Davis_Session.csv")

    card_id = "CalR_FIG2__Ontology__Grouped"
    source_meta = {
        "card_id": card_id,
        "resource_name": "Cal-Repository",
        "source_resource": "Cal-Repository",
        "source_dataset_unit": "Cal-Repository Figure 2 UC Davis CalR dataset",
        "dataset_unit_title": "Cal-Repository Figure 2 UC Davis indirect calorimetry dataset",
        "dataset_unit_type": "grouped_calorimetry_programs",
        "organism": "mouse",
        "comparison_space_organism": "human",
        "modality": "calorimetry",
        "tissue_or_system": "whole-animal metabolism",
        "contrast_label": "HFD and LFD group programs",
        "landing_page": CALR_LANDING_PAGE,
        "primary_access_url": CALR_REPOSITORY_URL,
        "access_route": "Public Cal-Repository Figure 2 files staged locally",
        "publication_ids": [],
        "focus_node": "Fig2_CalR_data.csv",
        "extractor_notes": "Reruns the public UC Davis Figure 2 CalR dataset through the ontology mapper using the current session parser and default packaged mouse calorimetry resources, with humanized output genes.",
        "source_files": [
            {
                "path": str(calr_data),
                "sha256": sha256_file(calr_data),
                "role": "calr_data_csv",
                "obtain_from_url": CALR_REPOSITORY_URL,
                "access_route": "Public Cal-Repository Figure 2 CalR data CSV",
                "landing_page_url": CALR_LANDING_PAGE,
            },
            {
                "path": str(session_csv),
                "sha256": sha256_file(session_csv),
                "role": "calr_session_csv",
                "obtain_from_url": CALR_REPOSITORY_URL,
                "access_route": "Public Cal-Repository Figure 2 CalR session CSV",
                "landing_page_url": CALR_LANDING_PAGE,
            },
        ],
        "preparation_steps": [],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Staged from the public Cal-Repository Figure 2 files:
- `raw/Fig2_CalR_data.csv`
- `raw/Fig2_CalR_Davis_Session.csv`

These files are rerun through `calr_ontology_mapper` with current shared provenance emission.
""",
    )


def build_calr_fig2_profile_query() -> None:
    slug = "calr_fig2_profile_query_public_grouped"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")

    calr_data = _stage_link(REPO_ROOT / "calorimetry_extractor" / "work" / "data" / "cal_repository" / "Fig2_CalR_data.csv", raw_dir / "Fig2_CalR_data.csv")
    session_csv = _stage_link(REPO_ROOT / "calorimetry_extractor" / "work" / "data" / "cal_repository" / "Fig2_CalR_Davis_Session.csv", raw_dir / "Fig2_CalR_Davis_Session.csv")
    studies_src = REPO_ROOT / "calorimetry_extractor" / "work" / "data" / "fig2_public_studies.tsv"
    studies_tsv = raw_dir / "fig2_public_studies.tsv"
    studies_df = pd.read_csv(studies_src, sep="\t")
    studies_df["calr_data_csv"] = str(calr_data)
    studies_df["session_csv"] = str(session_csv)
    studies_df.to_csv(studies_tsv, sep="\t", index=False)

    card_id = "CalR_FIG2__ProfileQuery__Grouped"
    source_meta = {
        "card_id": card_id,
        "resource_name": "Cal-Repository",
        "source_resource": "Cal-Repository",
        "source_dataset_unit": "Cal-Repository Figure 2 UC Davis public-bundle query run",
        "dataset_unit_title": "Cal-Repository Figure 2 UC Davis profile-query run against a public two-study bundle",
        "dataset_unit_type": "grouped_calorimetry_profile_query_programs",
        "organism": "mouse",
        "comparison_space_organism": "human",
        "modality": "calorimetry",
        "tissue_or_system": "whole-animal metabolism",
        "contrast_label": "HFD and LFD profile-query programs",
        "landing_page": CALR_LANDING_PAGE,
        "primary_access_url": CALR_REPOSITORY_URL,
        "access_route": "Public Cal-Repository Figure 2 files plus a small public studies manifest staged locally",
        "publication_ids": [],
        "focus_node": "fig2_public_studies.tsv",
        "extractor_notes": "Builds a tiny public calorimetry reference bundle from the Figure 2 HFD and LFD studies, then reruns the Figure 2 dataset through `calr_profile_query` with current provenance enabled.",
        "source_files": [
            {
                "path": str(calr_data),
                "sha256": sha256_file(calr_data),
                "role": "calr_data_csv",
                "obtain_from_url": CALR_REPOSITORY_URL,
                "access_route": "Public Cal-Repository Figure 2 CalR data CSV",
                "landing_page_url": CALR_LANDING_PAGE,
            },
            {
                "path": str(session_csv),
                "sha256": sha256_file(session_csv),
                "role": "calr_session_csv",
                "obtain_from_url": CALR_REPOSITORY_URL,
                "access_route": "Public Cal-Repository Figure 2 CalR session CSV",
                "landing_page_url": CALR_LANDING_PAGE,
            },
            {
                "path": str(studies_tsv),
                "sha256": sha256_file(studies_tsv),
                "role": "calr_public_studies_tsv",
                "obtain_from_url": CALR_REPOSITORY_URL,
                "access_route": "Local public studies manifest describing the Figure 2 HFD and LFD references",
                "landing_page_url": CALR_LANDING_PAGE,
            },
        ],
        "preparation_steps": [
            {
                "source_roles": ["calr_public_studies_tsv", "calr_data_csv", "calr_session_csv"],
                "target_role": "calr_public_studies_tsv",
                "label": "build public calorimetry reference bundle",
                "description": "Use the small Figure 2 public studies manifest to build a compact reference bundle for `calr_profile_query`.",
                "command": "geneset-extractors workflows calr_prepare_public --studies_tsv raw/fig2_public_studies.tsv ...",
                "script_path": str(REPO_ROOT / "calorimetry_extractor" / "work" / "data" / "fig2_public_studies.tsv"),
                "working_directory": str(REPO_ROOT),
            },
        ],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Staged from the public Cal-Repository Figure 2 files and a small public studies manifest:
- `raw/Fig2_CalR_data.csv`
- `raw/Fig2_CalR_Davis_Session.csv`
- `raw/fig2_public_studies.tsv`

The run first rebuilds a compact public reference bundle with `calr_prepare_public`, then queries that bundle with `calr_profile_query`.
""",
    )


def build_hubmap_liver_scrna() -> None:
    slug = "hubmap_liver_scrna_markers"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")
    staged_dir = ensure_dir(input_dir / "staged")

    processed_h5ad = raw_dir / "LV_processed.h5ad"
    counts_path = staged_dir / "hubmap_liver_counts.tsv"
    groups_path = staged_dir / "hubmap_liver_groups.tsv"

    _download(HUBMAP_LIVER_PROCESSED_URL, processed_h5ad)

    adata = ad.read_h5ad(processed_h5ad)
    matrix = adata.layers.get("unscaled", adata.X)
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(matrix)
    else:
        matrix = matrix.tocsr()

    barcodes = (
        adata.obs["barcode"].astype(str).to_numpy()
        if "barcode" in adata.obs.columns
        else adata.obs_names.astype(str).to_numpy()
    )
    groups = (
        adata.obs["leiden"].astype(str).to_numpy()
        if "leiden" in adata.obs.columns
        else np.full(adata.n_obs, "all", dtype=object)
    )
    gene_ids = _preferred_gene_ids(
        adata.var_names.astype(str).tolist(),
        adata.var["hugo_symbol"].astype(str).tolist() if "hugo_symbol" in adata.var.columns else adata.var_names.astype(str).tolist(),
    )

    if not counts_path.exists() or counts_path.stat().st_size == 0:
        _write_sparse_long_counts(matrix, gene_ids, np.asarray(barcodes, dtype=object), counts_path)
    pd.DataFrame({"barcode": barcodes, "group": groups}).to_csv(groups_path, sep="\t", index=False)

    card_id = "HuBMAP_Liver__scRNA__LeidenMarkers"
    source_meta = {
        "card_id": card_id,
        "resource_name": "HuBMAP",
        "source_resource": "HuBMAP",
        "source_dataset_unit": f"HuBMAP liver RNA-seq processed data product {HUBMAP_LIVER_PRODUCT_UUID}",
        "dataset_unit_title": f"HuBMAP liver RNA-seq processed data product {HUBMAP_LIVER_PRODUCT_UUID}",
        "dataset_unit_type": "single_cell_cluster_marker_collection",
        "organism": "human",
        "comparison_space_organism": "human",
        "modality": "single_cell_rna",
        "tissue_or_system": "Liver",
        "contrast_label": "published Leiden clusters from the processed HuBMAP liver RNA-seq data product",
        "landing_page": HUBMAP_LIVER_PAGE,
        "primary_access_url": HUBMAP_LIVER_PROCESSED_URL,
        "access_route": "Public HuBMAP data-products portal download",
        "publication_ids": [],
        "focus_node": "LV_processed.h5ad",
        "extractor_notes": "Uses the published HuBMAP liver processed H5AD, the provided Leiden cluster assignments, and the nonnegative `unscaled` layer to export a sparse long count table for `sc_rna_marker` without reclustering.",
        "source_files": [
            {
                "path": str(processed_h5ad),
                "sha256": sha256_file(processed_h5ad),
                "role": "hubmap_processed_h5ad",
                "obtain_from_url": HUBMAP_LIVER_PROCESSED_URL,
                "access_route": "Public HuBMAP processed H5AD download",
                "landing_page_url": HUBMAP_LIVER_PAGE,
            },
            {
                "path": str(counts_path),
                "sha256": sha256_file(counts_path),
                "role": "hubmap_marker_counts_tsv",
                "obtain_from_url": HUBMAP_LIVER_PROCESSED_URL,
                "access_route": "Locally staged sparse long count table derived from the public processed H5AD",
                "landing_page_url": HUBMAP_LIVER_PAGE,
            },
            {
                "path": str(groups_path),
                "sha256": sha256_file(groups_path),
                "role": "hubmap_marker_groups_tsv",
                "obtain_from_url": HUBMAP_LIVER_PROCESSED_URL,
                "access_route": "Locally staged Leiden-cluster table derived from the public processed H5AD",
                "landing_page_url": HUBMAP_LIVER_PAGE,
            },
        ],
        "preparation_steps": [
            {
                "source_roles": ["hubmap_processed_h5ad"],
                "target_role": "hubmap_marker_counts_tsv",
                "label": "export sparse count table",
                "description": "Export the processed HuBMAP liver H5AD `unscaled` layer into the sparse long TSV expected by the shared single-cell marker extractor.",
                "command": f"{PROTOTYPE_DIR / 'scripts' / 'build_additional_inputs.py'}",
                "script_path": str(PROTOTYPE_DIR / "scripts" / "build_additional_inputs.py"),
                "working_directory": str(REPO_ROOT),
            },
            {
                "source_roles": ["hubmap_processed_h5ad"],
                "target_role": "hubmap_marker_groups_tsv",
                "label": "export leiden groups",
                "description": "Export the published Leiden assignments from the processed HuBMAP liver H5AD into a barcode-to-group TSV.",
                "command": f"{PROTOTYPE_DIR / 'scripts' / 'build_additional_inputs.py'}",
                "script_path": str(PROTOTYPE_DIR / "scripts" / "build_additional_inputs.py"),
                "working_directory": str(REPO_ROOT),
            },
        ],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Downloaded from the public HuBMAP data-products portal:
- `raw/LV_processed.h5ad`

Staged extractor-ready files:
- `staged/hubmap_liver_counts.tsv`
- `staged/hubmap_liver_groups.tsv`

Preparation details:
- used HuBMAP liver processed RNA-seq data product `{HUBMAP_LIVER_PRODUCT_UUID}`
- reused the published `leiden` column as group labels
- exported the nonnegative `unscaled` layer into the sparse long TSV expected by `sc_rna_marker`
- retained all 4,433 processed cells and 21,218 features from the public processed object
""",
    )


def build_gse42752_methylation() -> None:
    slug = "gse42752_methylation_cpg_diff"
    input_dir = ensure_dir(PROTOTYPE_INPUTS_DIR / slug)
    raw_dir = ensure_dir(input_dir / "raw")

    cpg_path = _stage_link(
        REPO_ROOT / "dna_methylation_gene_extractors" / "analysis" / "methylation_validation" / "inputs" / "GSE42752_tumor_vs_matchednormal_cpg_diff.tsv",
        raw_dir / "GSE42752_tumor_vs_matchednormal_cpg_diff.tsv",
    )
    gtf_path = _stage_link(LOCAL_HG38_GTF, raw_dir / LOCAL_HG38_GTF.name)

    card_id = "GEO_GSE42752__Methylation__Tumor_vs_MatchedNormal"
    source_meta = {
        "card_id": card_id,
        "resource_name": "GEO",
        "source_resource": "GEO",
        "source_dataset_unit": "GSE42752 tumor versus matched normal methylation contrast",
        "dataset_unit_title": "GSE42752 differential methylation: tumor versus matched normal",
        "dataset_unit_type": "cpg_level_case_control_signature",
        "organism": "human",
        "comparison_space_organism": "human",
        "modality": "dna_methylation",
        "tissue_or_system": "head and neck tumor",
        "contrast_label": "tumor versus matched normal",
        "landing_page": GSE42752_PAGE,
        "primary_access_url": GSE42752_PAGE,
        "access_route": "Public GEO series staged locally as a differential CpG table",
        "publication_ids": [],
        "focus_node": "GSE42752_tumor_vs_matchednormal_cpg_diff.tsv",
        "extractor_notes": "Uses the staged public GSE42752 tumor-versus-matched-normal differential CpG table from the methylation validation workspace and reruns it through the shared methylation extractor with current provenance enabled.",
        "source_files": [
            {
                "path": str(cpg_path),
                "sha256": sha256_file(cpg_path),
                "role": "methylation_cpg_diff_tsv",
                "obtain_from_url": GSE42752_PAGE,
                "access_route": "Public GEO methylation differential table staged locally from the validation workspace",
                "landing_page_url": GSE42752_PAGE,
            },
            {
                "path": str(gtf_path),
                "sha256": sha256_file(gtf_path),
                "role": "reference_gtf",
                "obtain_from_url": GENCODE_V47_URL,
                "access_route": "GENCODE release 47 human annotation",
                "landing_page_url": "https://www.gencodegenes.org/human/release_47.html",
            },
        ],
        "preparation_steps": [],
    }
    write_json(input_dir / "source_meta.json", source_meta)
    _write_overlay(input_dir, source_meta)
    _write_readme(
        input_dir / "README.md",
        f"""
# {slug}

Staged from a public GEO methylation contrast already prepared locally:
- `raw/GSE42752_tumor_vs_matchednormal_cpg_diff.tsv`
- `raw/{gtf_path.name}`

The CpG differential table is rerun through `methylation_cpg_diff` with explicit `delta_beta`, `pval`, and `qval` column mappings.
""",
    )


def main() -> None:
    build_gtex_grouped_bulk_rna()
    build_motrpac_wat_proteomics()
    build_motrpac_wat_ptm()
    build_4dn_atac()
    build_4dn_chipseq()
    build_tcga_brca_splicing()
    build_prism_quick_drug_response()
    build_tcga_brca_cnv()
    build_jump_cellpainting_morphology()
    build_calr_fig2_ontology()
    build_calr_fig2_profile_query()
    build_hubmap_liver_scrna()
    build_gse42752_methylation()


if __name__ == "__main__":
    main()
