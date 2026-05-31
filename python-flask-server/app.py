from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import pandas as pd
import streamlit as st


APP_DIR = Path(__file__).resolve().parent
SRC_DIR = APP_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from retrieval import RevealRetrievalIndex, normalize_gene_list  # noqa: E402

try:
    from streamlit_agraph import Edge, Node, _agraph

    HAVE_AGRAPH = True
except ImportError:
    HAVE_AGRAPH = False


st.set_page_config(page_title="CFDE-REVEAL Prototype", layout="wide")

st.markdown(
    """
    <style>
    :root {
      --portal-orange: #ff6600;
      --portal-blue: #35669a;
      --portal-gray-100: #f8f8f8;
      --portal-gray-150: #f1f1f1;
      --portal-gray-200: #eeeeee;
      --portal-gray-300: #dddddd;
      --portal-gray-500: #777777;
      --portal-border: #d9d9d9;
      --portal-text: #1f2933;
    }
    html, body, [class*="css"] {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
      font-size: 14px;
    }
    .stApp {
      background: #ffffff;
      color: var(--portal-text);
    }
    .block-container {
      padding-top: 2rem;
      padding-bottom: 3rem;
    }
    h1, h2, h3 {
      font-family: inherit;
      letter-spacing: normal;
      color: var(--portal-text);
    }
    .portal-hero {
      display: flex;
      flex-direction: column;
      gap: 0.35rem;
      margin-bottom: 1.5rem;
    }
    .portal-brand {
      font-size: 1.75em;
      font-weight: 700;
      line-height: 1em;
      color: var(--portal-text);
      max-width: 100%;
      overflow-wrap: anywhere;
    }
    .hero-kicker {
      font-size: 0.95rem;
      font-weight: 700;
      color: var(--portal-orange);
      margin: 0;
    }
    .hero-title {
      font-size: 1.2em;
      font-weight: 700;
      color: var(--portal-text);
      margin: 0;
    }
    .hero-subtitle {
      font-size: 1rem;
      color: var(--portal-text);
      max-width: 56rem;
      margin: 0;
    }
    .chip-row {
      display: flex;
      gap: 0.5rem;
      flex-wrap: wrap;
      margin: 0.5rem 0 1rem 0;
    }
    .chip {
      background: var(--portal-gray-200);
      border: 0.5px solid var(--portal-gray-300);
      border-radius: 10px;
      padding: 0.15rem 0.65rem;
      font-size: 0.95rem;
    }
    .result-title {
      font-size: 1.2rem;
      font-weight: 700;
      margin-bottom: 0.35rem;
    }
    .section-lead {
      color: var(--portal-orange);
      font-size: 1.2em;
      font-weight: 700;
      margin: 0 0 0.35rem 0;
    }
    .portal-note {
      display: block;
      color: var(--portal-gray-500);
      margin: 4px 0 8px 0;
      padding: 4px 8px;
      background-color: var(--portal-gray-100);
      border-left: 3px solid #7c757d;
    }
    div[data-testid="stVerticalBlock"] > div[data-testid="stContainer"] {
      background: #ffffff;
      border: 1px solid var(--portal-border);
      border-radius: 8px;
      padding: 0.8rem 1rem 1rem;
      box-shadow: none;
    }
    div[data-testid="stExpander"] {
      border: none;
    }
    div[data-testid="stExpander"] details {
      border: none;
      background: transparent;
    }
    div[data-testid="stExpander"] summary {
      background: var(--portal-gray-200);
      border-radius: 5px;
      padding: 10px 12px;
    }
    div[data-testid="stExpander"] summary:hover {
      background: var(--portal-gray-300);
    }
    div[data-testid="stExpander"] details[open] summary {
      margin-bottom: 0.75rem;
    }
    div.stButton > button, div.stDownloadButton > button {
      background: var(--portal-blue);
      color: #ffffff;
      border: 1px solid var(--portal-blue);
      border-radius: 0.25rem;
    }
    div.stButton > button:hover, div.stDownloadButton > button:hover {
      background: #2b5680;
      border-color: #2b5680;
      color: #ffffff;
    }
    div[data-testid="stTextArea"] textarea, div[data-testid="stFileUploader"] section {
      border-color: var(--portal-border);
    }
    </style>
    """,
    unsafe_allow_html=True,
)


GRAPH_COLORS = {
    "source": "#35669a",
    "intermediate": "#ff6600",
    "output": "#7c757d",
}

EDGE_COLORS = {
    "source_preparation": "#ff6600",
    "extractor_conversion": "#35669a",
}


@st.cache_resource
def load_index() -> RevealRetrievalIndex:
    return RevealRetrievalIndex(APP_DIR / "data")


@st.cache_data
def load_manifest() -> pd.DataFrame:
    return pd.read_parquet(APP_DIR.parent / "data" / "card_manifest.parquet")


def resolve_bundle_path(path: str) -> Path:
    raw = Path(path)
    return raw if raw.is_absolute() else (APP_DIR.parent / raw)


@st.cache_data
def load_json(path: str) -> dict[str, Any]:
    return json.loads(resolve_bundle_path(path).read_text(encoding="utf-8"))


def _format_species(row: pd.Series) -> str:
    native = str(row["organism"])
    comparison = str(row["comparison_space_organism"])
    if native == comparison:
        return native
    return f"{native} source, {comparison} comparison space"


def _badge_values(row: pd.Series) -> list[str]:
    return [
        str(row["resource_name"]),
        str(row["modality"]),
        str(row["tissue_or_system"]),
        _format_species(row),
    ]


def _badge_row(values: list[str]) -> str:
    return "".join(f'<span class="chip">{value}</span>' for value in values if value)


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


def _download_payload(row: pd.Series, overlay: dict[str, Any], graph: dict[str, Any]) -> dict[str, Any]:
    payload = row.to_dict()
    payload["provenance"] = overlay
    payload["provenance_graph"] = graph
    payload["shared_genes"] = [gene for gene in str(row.get("overlapping_genes", "")).split(",") if gene]
    return payload


def _node_title(node: dict[str, Any]) -> str:
    lines = [
        node["label"],
        f"Stage: {node.get('stage', 'unknown')}",
        f"Role(s): {node.get('role_summary', '')}",
    ]
    if node.get("path"):
        lines.append(str(node["path"]))
    return "\n".join(lines)


def _graph_nodes(graph: dict[str, Any]) -> list[Node]:
    out: list[Node] = []
    for node in graph["nodes"]:
        stage = str(node.get("stage", "source"))
        out.append(
            Node(
                id=str(node["id"]),
                label=str(node["label"]),
                title=_node_title(node),
                shape="box" if stage != "output" else "ellipse",
                size=24 if stage != "output" else 22,
                color=GRAPH_COLORS.get(stage, "#6b7280"),
            )
        )
    return out


def _graph_edges(graph: dict[str, Any]) -> list[Edge]:
    out: list[Edge] = []
    dense_graph = len(graph["edges"]) > 40
    for edge in graph["edges"]:
        out.append(
            Edge(
                source=str(edge["source"]),
                target=str(edge["target"]),
                label="" if dense_graph else str(edge["label"]),
                color=EDGE_COLORS.get(str(edge.get("edge_kind", "")), "#7c8793"),
                smooth=not dense_graph,
            )
        )
    return out


def _graph_config(graph: dict[str, Any]) -> dict[str, Any]:
    dense_graph = len(graph["edges"]) > 40
    return {
        "width": "1100px",
        "height": "550px",
        "autoResize": True,
        "layout": {
            "hierarchical": {
                "enabled": dense_graph,
                "direction": "LR",
                "sortMethod": "directed",
                "nodeSpacing": 180,
                "levelSeparation": 220,
            },
            "improvedLayout": True,
            "randomSeed": 17,
        },
        "interaction": {
            "hover": True,
            "tooltipDelay": 120,
            "dragView": True,
            "zoomView": True,
        },
        "physics": {
            "enabled": False,
            "stabilization": {
                "enabled": False,
            },
        },
        "edges": {
            "font": {"size": 10},
        },
    }


def _render_agraph(graph: dict[str, Any], key: str) -> str | None:
    if not HAVE_AGRAPH:
        return None
    nodes_data = [node.to_dict() for node in _graph_nodes(graph)]
    edges_data = [edge.to_dict() for edge in _graph_edges(graph)]
    config_json = json.dumps(_graph_config(graph))
    data_json = json.dumps({"nodes": nodes_data, "edges": edges_data})
    return _agraph(data=data_json, config=config_json, key=key)


def _select_defaults(graph: dict[str, Any], key_prefix: str, clicked: str | None, node_map: dict[str, Any], edge_map: dict[str, Any]) -> tuple[str, str]:
    node_key = f"{key_prefix}_node"
    edge_key = f"{key_prefix}_edge"
    if node_key not in st.session_state:
        st.session_state[node_key] = graph.get("default_node_id") or next(iter(node_map))
    if edge_key not in st.session_state:
        st.session_state[edge_key] = graph.get("default_edge_id") or next(iter(edge_map))
    if clicked in node_map:
        st.session_state[node_key] = clicked
    elif clicked in edge_map:
        st.session_state[edge_key] = clicked
    return st.session_state[node_key], st.session_state[edge_key]


def _render_node_detail(node: dict[str, Any], key_prefix: str) -> None:
    st.markdown(f"**{node['label']}**")
    st.write(node.get("description", ""))
    st.write(f"Stage: {node.get('stage', 'unknown')}")
    st.write(f"Role(s): {node.get('role_summary', 'n/a')}")
    st.code(str(node.get("path", "")), language="text")
    if node.get("sha256"):
        st.caption(f"sha256: {node['sha256']}")
    if node.get("obtain_from_url"):
        st.markdown(f"[Obtain this file]({node['obtain_from_url']})")
    if node.get("landing_page_url") and node.get("landing_page_url") != node.get("obtain_from_url"):
        st.markdown(f"[Dataset page]({node['landing_page_url']})")
    st.caption(node.get("access_route", ""))


def _render_edge_detail(edge: dict[str, Any], key_prefix: str) -> None:
    st.markdown(f"**{edge['label']}**")
    st.write(edge.get("description", ""))
    st.code(str(edge.get("command", "")), language="bash")
    if edge.get("working_directory"):
        st.caption(f"Working directory: {edge['working_directory']}")
    if edge.get("script_path"):
        st.markdown(f"Script path: `{edge['script_path']}`")
        script_path = resolve_bundle_path(str(edge["script_path"]))
        if script_path.exists():
            with st.expander("Preview script source", expanded=False):
                try:
                    st.code(script_path.read_text(encoding="utf-8"), language="python")
                except UnicodeDecodeError:
                    st.write("Script exists but could not be decoded as UTF-8.")
    if edge.get("notebook_url"):
        st.markdown(f"[Notebook reference]({edge['notebook_url']})")
    if edge.get("parameters"):
        with st.expander("Command parameters", expanded=False):
            st.json(edge["parameters"])


def render_provenance_panel(graph: dict[str, Any], key_prefix: str) -> None:
    node_map = {str(node["id"]): node for node in graph["nodes"]}
    edge_map = {str(edge["id"]): edge for edge in graph["edges"]}

    clicked: str | None = None
    if HAVE_AGRAPH:
        clicked = _render_agraph(graph, key=f"{key_prefix}_agraph")
        st.caption("Use the graph to focus a node. Step details are available in the processing-step selector below.")
    else:
        st.info("Install `streamlit-agraph` to enable the interactive graph. The selector-based provenance details remain available.")

    default_node_id, default_edge_id = _select_defaults(graph, key_prefix, clicked, node_map, edge_map)

    node_ids = list(node_map)
    edge_ids = list(edge_map)
    selected_node_id = st.selectbox(
        "File node",
        node_ids,
        index=node_ids.index(default_node_id),
        key=f"{key_prefix}_node",
        format_func=lambda value: f"{node_map[value]['label']} [{node_map[value]['stage']}]",
    )
    selected_edge_id = st.selectbox(
        "Processing step",
        edge_ids,
        index=edge_ids.index(default_edge_id),
        key=f"{key_prefix}_edge",
        format_func=lambda value: edge_map[value]["label"],
    )

    left, right = st.columns(2)
    with left:
        _render_node_detail(node_map[selected_node_id], key_prefix)
    with right:
        _render_edge_detail(edge_map[selected_edge_id], key_prefix)


def main() -> None:
    idx = load_index()
    manifest = load_manifest()

    st.markdown(
        """
        <div class="portal-hero">
          <div class="portal-brand">CFDE REVEAL dataset search</div>
          <div class="hero-kicker">Turn your research question into reusable dataset leads</div>
          <div class="hero-title">Find original datasets linked to your gene list, then inspect the evidence and provenance behind the match.</div>
          <div class="hero-subtitle">Explore dataset-first matches from GTEx, MoTrPAC, and GEO. Each result foregrounds the source dataset and exposes the tracked file lineage used to derive the supporting signature.</div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="chip-row"><div class="chip">GTEx</div><div class="chip">MoTrPAC</div><div class="chip">GEO</div><div class="chip">Dataset-first ranking</div><div class="chip">File-level provenance graph</div></div>',
        unsafe_allow_html=True,
    )

    with st.sidebar:
        st.subheader("Filters")
        resource = st.selectbox("Resource", ["All"] + sorted(manifest["resource_name"].unique().tolist()))
        modality = st.selectbox("Modality", ["All"] + sorted(manifest["modality"].unique().tolist()))
        tissue = st.selectbox("Tissue / system", ["All"] + sorted(manifest["tissue_or_system"].unique().tolist()))

    st.markdown('<div class="section-lead">Tell us what you\'re studying or curious about:</div>', unsafe_allow_html=True)
    left, right = st.columns([3, 2])
    with left:
        gene_text = st.text_area("Paste gene list", height=180, placeholder="ACTA1\nCKM\nMYH1\n...")
    with right:
        upload = st.file_uploader("Upload text file", type=["txt", "tsv", "csv"])
        if upload is not None and not gene_text.strip():
            gene_text = upload.getvalue().decode("utf-8")
            st.caption("Loaded genes from uploaded file.")

    query_genes = normalize_gene_list(gene_text)
    st.markdown(f'<div class="portal-note">Detected genes: {len(query_genes)}</div>', unsafe_allow_html=True)

    filters = {
        "resource_name": None if resource == "All" else resource,
        "modality": None if modality == "All" else modality,
        "tissue_or_system": None if tissue == "All" else tissue,
    }

    if not query_genes:
        st.info("Enter or upload a gene list to rank relevant datasets.")
        return

    results = idx.score_query(query_genes, filters=filters, top_k=20)
    if results.empty:
        st.warning("No cards matched the current filters.")
        return

    st.markdown('<div class="section-lead">Ranked original datasets</div>', unsafe_allow_html=True)
    for _, row in results.iterrows():
        overlay = load_json(str(row["provenance_path"]))
        graph = load_json(str(row["provenance_graph_path"]))
        rationale = build_dataset_rationale(row)
        payload = _download_payload(row, overlay, graph)

        with st.container(border=True):
            st.markdown(f'<div class="result-title">{row["dataset_unit_title"]}</div>', unsafe_allow_html=True)
            st.markdown(f'<div class="chip-row">{_badge_row(_badge_values(row))}</div>', unsafe_allow_html=True)
            st.markdown(f"**Why this original dataset is relevant**  {rationale['relevance']}")
            st.markdown(f"**What you can reuse here**  {rationale['reuse']}")
            st.markdown(
                "**Source entry point**  "
                f"[Dataset page]({row['landing_page']})"
            )
            if overlay.get("primary_access_url"):
                st.markdown(
                    "**Primary data access**  "
                    f"[Open acquisition URL]({overlay['primary_access_url']})"
                )
            st.caption(str(row["access_route"]))

            with st.expander("Shared biological signal and derived evidence", expanded=False):
                st.write(f"Overlap score: {row['overlap_score']:.4f}")
                st.write(f"Enrichment-like score: {row['enrichment_score']:.4f}")
                st.write(f"Latent signature similarity: {row['latent_score']:.4f}")
                st.write(f"Overlapping genes: {row['overlapping_genes'] or 'none'}")
                st.caption(rationale["evidence"])

            with st.expander("Provenance graph", expanded=False):
                render_provenance_panel(graph, key_prefix=f"card_{row['card_id']}")

            dl1, dl2, dl3, dl4 = st.columns(4)
            with dl1:
                st.download_button(
                    "Signature TSV",
                    resolve_bundle_path(str(row["signature_path"])).read_bytes(),
                    file_name=f"{row['card_id']}.geneset.tsv",
                    key=f"dl_sig_{row['card_id']}",
                )
            with dl2:
                st.download_button(
                    "Overlay JSON",
                    json.dumps(overlay, indent=2).encode("utf-8"),
                    file_name=f"{row['card_id']}.overlay.json",
                    key=f"dl_overlay_{row['card_id']}",
                )
            with dl3:
                st.download_button(
                    "Graph JSON",
                    json.dumps(graph, indent=2).encode("utf-8"),
                    file_name=f"{row['card_id']}.graph.json",
                    key=f"dl_graph_{row['card_id']}",
                )
            with dl4:
                st.download_button(
                    "Card JSON",
                    json.dumps(payload, indent=2, default=str).encode("utf-8"),
                    file_name=f"{row['card_id']}.card.json",
                    key=f"dl_card_{row['card_id']}",
                )

    with st.expander("Result table", expanded=False):
        st.dataframe(
            results[
                [
                    "card_id",
                    "resource_name",
                    "tissue_or_system",
                    "overlap_score",
                    "enrichment_score",
                    "latent_score",
                    "final_score",
                ]
            ],
            width="stretch",
        )

    with st.expander("Full provenance viewer", expanded=False):
        chosen = st.selectbox("Select card", results["card_id"].tolist())
        selected = results[results["card_id"] == chosen].iloc[0]
        graph = load_json(str(selected["provenance_graph_path"]))
        render_provenance_panel(graph, key_prefix="full_viewer")


if __name__ == "__main__":
    main()
