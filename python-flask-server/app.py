from __future__ import annotations

import json
import math
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import pandas as pd
from flask import Flask, abort, render_template, request, send_file


APP_DIR = Path(__file__).resolve().parent
SRC_DIR = APP_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from retrieval import RevealRetrievalIndex, normalize_gene_list  # noqa: E402


app = Flask(__name__)

GENESET_API_ROOT = "https://translator.broadinstitute.org/genetics_provider/geneset_extractor"
GENESET_LIST_URL = f"{GENESET_API_ROOT}/gene-sets"
GENESET_DETAIL_URL = f"{GENESET_API_ROOT}/gene-set"
GENESET_PROVENANCE_URL = f"{GENESET_API_ROOT}/gene_set_provenance"


@lru_cache(maxsize=1)
def load_index() -> RevealRetrievalIndex:
    return RevealRetrievalIndex(APP_DIR / "data")


@lru_cache(maxsize=512)
def load_json(path: str) -> dict[str, Any]:
    return json.loads(resolve_bundle_path(path).read_text(encoding="utf-8"))


def resolve_bundle_path(path: str) -> Path:
    raw = Path(path)
    return raw if raw.is_absolute() else (APP_DIR.parent / raw)


def _extract_gene_set_records(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]
    if isinstance(payload, dict):
        for key in ("gene_sets", "genesets", "items", "results", "data"):
            value = payload.get(key)
            if isinstance(value, list):
                return [item for item in value if isinstance(item, dict)]
        return [payload]
    return []


@lru_cache(maxsize=1)
def load_gene_sets() -> list[dict[str, Any]]:
    payload = fetch_remote_json(GENESET_LIST_URL)
    return _extract_gene_set_records(payload)


@lru_cache(maxsize=512)
def fetch_remote_json(url: str) -> dict[str, Any] | list[Any]:
    with urlopen(url, timeout=20) as response:
        return json.loads(response.read().decode("utf-8"))


@lru_cache(maxsize=512)
def load_gene_set_detail(gene_set_id: int) -> dict[str, Any]:
    payload = fetch_remote_json(f"{GENESET_DETAIL_URL}?gene_set_id={gene_set_id}")
    if not isinstance(payload, dict):
        raise ValueError(f"Unexpected detail payload for gene_set_id={gene_set_id}")
    return payload


@lru_cache(maxsize=512)
def load_gene_set_provenance(gene_set_id: int) -> dict[str, Any]:
    payload = fetch_remote_json(f"{GENESET_PROVENANCE_URL}?gene_set_id={gene_set_id}")
    if not isinstance(payload, dict):
        raise ValueError(f"Unexpected provenance payload for gene_set_id={gene_set_id}")
    return payload


def _first_present(record: dict[str, Any], keys: tuple[str, ...]) -> str | None:
    for key in keys:
        value = record.get(key)
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return text
    return None


def build_gene_set_summary(gene_sets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for record in gene_sets:
        name = _first_present(
            record,
            ("standard_name", "name", "label", "title", "geneset_name", "gene_set_name", "gene_set_id", "id"),
        )
        if not name:
            continue
        tags = record.get("tags")
        if isinstance(tags, list):
            context = ", ".join(str(tag) for tag in tags if str(tag).strip()) or "Unknown"
        else:
            context = _first_present(record, ("context", "tissue_or_system", "system", "organism", "tags")) or "Unknown"
        summaries.append(
            {
                "gene_set_id": int(record.get("gene_set_id")) if record.get("gene_set_id") is not None else None,
                "name": name,
                "source": _first_present(record, ("collection_name", "source", "resource_name", "resource", "provider")) or "Unknown",
                "category": _first_present(record, ("license_code", "category", "modality", "type")) or "Unknown",
                "context": context,
            }
        )
    summaries.sort(key=lambda item: item["name"].lower())
    return summaries


def build_filter_options(gene_sets: list[dict[str, Any]]) -> tuple[list[str], list[str], list[str]]:
    resource_options = sorted({item["source"] for item in gene_sets if item["source"] != "Unknown"})
    modality_options = sorted({item["category"] for item in gene_sets if item["category"] != "Unknown"})
    tissue_options = sorted({item["context"] for item in gene_sets if item["context"] != "Unknown"})
    return resource_options, modality_options, tissue_options


def build_graph_view(knowledge_graph: dict[str, Any]) -> dict[str, Any]:
    nodes = knowledge_graph.get("nodes", [])
    edges = knowledge_graph.get("edges", [])
    width = 900
    height = 520
    cx = width / 2
    cy = height / 2
    radius = max(120, min(width, height) * 0.34)
    count = max(len(nodes), 1)
    positions: dict[str, tuple[float, float]] = {}

    for idx, node in enumerate(nodes):
        angle = (2 * math.pi * idx / count) - (math.pi / 2)
        positions[str(node["id"])] = (
            cx + radius * math.cos(angle),
            cy + radius * math.sin(angle),
        )

    svg_edges: list[dict[str, Any]] = []
    for edge in edges:
        source = positions.get(str(edge.get("source")))
        target = positions.get(str(edge.get("target")))
        if not source or not target:
            continue
        svg_edges.append(
            {
                "x1": source[0],
                "y1": source[1],
                "x2": target[0],
                "y2": target[1],
                "label": str(edge.get("label", "")),
            }
        )

    svg_nodes: list[dict[str, Any]] = []
    for node in nodes:
        node_id = str(node["id"])
        x, y = positions[node_id]
        node_type = str(node.get("type", "Node"))
        if node_type == "GeneSet":
            color = "#ff6600"
        elif "Analysis" in node_type:
            color = "#35669a"
        else:
            color = "#7c757d"
        svg_nodes.append(
            {
                "id": node_id,
                "label": str(node.get("name") or node.get("label") or node_id),
                "type": node_type,
                "description": str(node.get("description", "")),
                "x": x,
                "y": y,
                "color": color,
            }
        )

    return {
        "width": width,
        "height": height,
        "nodes": svg_nodes,
        "edges": svg_edges,
    }


def build_gene_set_page(detail: dict[str, Any], provenance: dict[str, Any]) -> dict[str, Any]:
    knowledge_graph = provenance.get("knowledge_graph")
    if not isinstance(knowledge_graph, dict):
        knowledge_graph = detail.get("knowledge_graph", {})

    genes = detail.get("gene_symbols", [])
    gene_symbols = [str(item.get("symbol")) for item in genes if isinstance(item, dict) and item.get("symbol")]

    return {
        "gene_set_id": detail.get("gene_set_id"),
        "standard_name": detail.get("standard_name") or detail.get("card_id"),
        "collection_name": detail.get("collection_name") or detail.get("resource_name"),
        "license_code": detail.get("license_code"),
        "card_id": detail.get("card_id"),
        "dataset_unit_title": detail.get("dataset_unit_title"),
        "contrast_label": detail.get("contrast_label"),
        "organism": detail.get("organism"),
        "comparison_space_organism": detail.get("comparison_space_organism"),
        "modality": detail.get("modality"),
        "gene_symbols": gene_symbols,
        "graph": build_graph_view(knowledge_graph),
        "graph_nodes": knowledge_graph.get("nodes", []),
        "graph_edges": knowledge_graph.get("edges", []),
        "detail_json": detail,
        "provenance_json": provenance,
        "detail_url": f"{GENESET_DETAIL_URL}?gene_set_id={detail.get('gene_set_id')}",
        "provenance_url": f"{GENESET_PROVENANCE_URL}?gene_set_id={detail.get('gene_set_id')}",
    }


def format_species(row: pd.Series) -> str:
    native = str(row["organism"])
    comparison = str(row["comparison_space_organism"])
    if native == comparison:
        return native
    return f"{native} source, {comparison} comparison space"


def badge_values(row: pd.Series) -> list[str]:
    return [
        str(row["resource_name"]),
        str(row["modality"]),
        str(row["tissue_or_system"]),
        format_species(row),
    ]


def build_dataset_rationale(row: pd.Series) -> dict[str, str]:
    overlap_genes = [gene for gene in str(row.get("overlapping_genes", "")).split(",") if gene][:5]
    overlap_text = ", ".join(overlap_genes) if overlap_genes else "broad shared signal without strong direct overlap"
    relevance = (
        f"This original {row['resource_name']} dataset contains {row['modality']} signal from "
        f"{row['tissue_or_system']} for the contrast '{row['contrast_label']}', which aligns with your query through "
        f"{overlap_text}."
    )
    reuse = (
        f"Use the original dataset entry point at {row['resource_name']} first. The portal exposes "
        f"{row['source_file_count']} tracked file artifacts behind this card and links them back to the source landing page."
    )
    evidence = (
        f"The ranking itself is driven by a derived {int(row['n_signature_genes'])}-gene signature in "
        f"{row['comparison_space_organism']} comparison space."
    )
    return {
        "relevance": relevance,
        "reuse": reuse,
        "evidence": evidence,
    }


def build_download_payload(row: pd.Series, overlay: dict[str, Any], graph: dict[str, Any]) -> dict[str, Any]:
    payload = row.to_dict()
    payload["provenance"] = overlay
    payload["provenance_graph"] = graph
    payload["shared_genes"] = [gene for gene in str(row.get("overlapping_genes", "")).split(",") if gene]
    return payload


def build_card_view(row: pd.Series) -> dict[str, Any]:
    overlay = load_json(str(row["provenance_path"]))
    graph = load_json(str(row["provenance_graph_path"]))
    rationale = build_dataset_rationale(row)
    node_map = {str(node["id"]): node for node in graph["nodes"]}
    edge_map = {str(edge["id"]): edge for edge in graph["edges"]}

    default_node_id = str(graph.get("default_node_id") or next(iter(node_map), ""))
    default_edge_id = str(graph.get("default_edge_id") or next(iter(edge_map), ""))

    return {
        "card_id": str(row["card_id"]),
        "title": str(row["dataset_unit_title"]),
        "badges": badge_values(row),
        "rationale": rationale,
        "landing_page": str(row["landing_page"]),
        "access_route": str(row["access_route"]),
        "primary_access_url": overlay.get("primary_access_url"),
        "scores": {
            "overlap_score": float(row["overlap_score"]),
            "enrichment_score": float(row["enrichment_score"]),
            "latent_score": float(row["latent_score"]),
            "overlapping_genes": str(row["overlapping_genes"] or "none"),
        },
        "overlay": overlay,
        "graph": graph,
        "default_node": node_map.get(default_node_id),
        "default_edge": edge_map.get(default_edge_id),
        "nodes": list(node_map.values()),
        "edges": list(edge_map.values()),
        "downloads": {
            "signature": str(row["signature_path"]),
            "overlay": overlay,
            "graph": graph,
            "card": build_download_payload(row, overlay, graph),
        },
    }


def extract_query_text() -> str:
    gene_text = request.form.get("gene_text", "")
    upload = request.files.get("gene_file")
    if upload and not gene_text.strip():
        return upload.read().decode("utf-8")
    return gene_text


@app.route("/", methods=["GET", "POST"])
def index() -> str:
    gene_set_summaries: list[dict[str, str]] = []
    resource_options: list[str] = []
    modality_options: list[str] = []
    tissue_options: list[str] = []

    gene_text = ""
    query_genes: list[str] = []
    selected_resource = "All"
    selected_modality = "All"
    selected_tissue = "All"
    cards: list[dict[str, Any]] = []
    result_rows: list[dict[str, Any]] = []
    message: tuple[str, str] | None = None

    try:
        gene_set_summaries = build_gene_set_summary(load_gene_sets())
        resource_options, modality_options, tissue_options = build_filter_options(gene_set_summaries)
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError) as exc:
        message = ("warning", f"Could not load gene sets from {GENESET_LIST_URL}: {exc}")

    if request.method == "POST":
        gene_text = extract_query_text()
        query_genes = normalize_gene_list(gene_text)
        selected_resource = request.form.get("resource", "All")
        selected_modality = request.form.get("modality", "All")
        selected_tissue = request.form.get("tissue", "All")

        if not query_genes:
            message = ("info", "Enter or upload a gene list to rank relevant datasets.")
        else:
            filters = {
                "resource_name": None if selected_resource == "All" else selected_resource,
                "modality": None if selected_modality == "All" else selected_modality,
                "tissue_or_system": None if selected_tissue == "All" else selected_tissue,
            }
            results = load_index().score_query(query_genes, filters=filters, top_k=20)
            if results.empty:
                message = ("warning", "No cards matched the current filters.")
            else:
                cards = [build_card_view(row) for _, row in results.iterrows()]
                result_rows = results[
                    [
                        "card_id",
                        "resource_name",
                        "tissue_or_system",
                        "overlap_score",
                        "enrichment_score",
                        "latent_score",
                        "final_score",
                    ]
                ].to_dict(orient="records")

    return render_template(
        "index.html",
        resource_options=resource_options,
        modality_options=modality_options,
        tissue_options=tissue_options,
        gene_text=gene_text,
        query_genes=query_genes,
        selected_resource=selected_resource,
        selected_modality=selected_modality,
        selected_tissue=selected_tissue,
        cards=cards,
        result_rows=result_rows,
        gene_set_summaries=gene_set_summaries,
        gene_set_list_url=GENESET_LIST_URL,
        message=message,
    )


@app.route("/gene-sets/<int:gene_set_id>")
def gene_set_detail(gene_set_id: int) -> str:
    try:
        detail = load_gene_set_detail(gene_set_id)
        provenance = load_gene_set_provenance(gene_set_id)
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError, ValueError) as exc:
        abort(502, description=f"Could not load gene set {gene_set_id}: {exc}")

    return render_template(
        "gene_set_detail.html",
        page=build_gene_set_page(detail, provenance),
    )


@app.route("/download/<card_id>/<kind>")
def download(card_id: str, kind: str):
    idx = load_index()
    matches = idx.cards[idx.cards["card_id"].astype(str) == card_id]
    if matches.empty:
        abort(404)

    row = matches.iloc[0]
    overlay = load_json(str(row["provenance_path"]))
    graph = load_json(str(row["provenance_graph_path"]))

    if kind == "signature":
        path = resolve_bundle_path(str(row["signature_path"]))
        return send_file(path, as_attachment=True, download_name=f"{card_id}.geneset.tsv")

    payloads = {
        "overlay": (overlay, f"{card_id}.overlay.json"),
        "graph": (graph, f"{card_id}.graph.json"),
        "card": (build_download_payload(row, overlay, graph), f"{card_id}.card.json"),
    }
    if kind not in payloads:
        abort(404)

    payload, filename = payloads[kind]
    return app.response_class(
        json.dumps(payload, indent=2, default=str),
        mimetype="application/json",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8000, debug=True)
