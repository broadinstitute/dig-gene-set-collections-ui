# Frozen Bundle Plan

## Goal

Create a portable frozen handoff directory under `ui_test/` that is sufficient to:

- run the CFDE REVEAL Streamlit app locally without the rest of this workspace
- include the minimum generated data required for the current 48-card portal
- provide enough context and scripts for a collaborator to add more downloaded gene-set cards later and rebuild the local index

## Bundle contents

- root `README.md` describing what is included, how to run, and how to extend
- `app/` with the runnable Streamlit entrypoint, `src/`, and `requirements.txt`
- `data/cards/` with normalized per-card extracted artifacts:
  - `geneset.tsv`
  - `geneset.meta.json`
  - `geneset.provenance.json`
- `data/source_meta/` with one `source_meta.json` per card
- `data/provenance/` with overlay JSON
- `data/provenance_graph/` with graph JSON
- `data/card_manifest.parquet`
- `data/signature_index.parquet`
- `app/data/` with the prebuilt DuckDB and latent-model artifacts needed to run the app
- `inventory/cards.tsv` describing every bundled card and where its source metadata and extracted artifacts live
- `scripts/` with the minimum local rebuild path for manifest and retrieval index refresh
- `context/` with concise notes, plans, and run inventory sufficient for a new maintainer to understand the prototype scope

## Implementation steps

1. Build a portable bundle layout and normalize each card into a per-card directory with local relative paths.
2. Copy the app code and adapt it to resolve bundle-local paths rather than absolute workspace paths.
3. Add bundle-local rebuild scripts that regenerate manifest, graph overlays, DuckDB, and latent index from the bundled inventory.
4. Copy the minimal context/history files that explain the prototype scope, current corpus, and extension workflow.
5. Validate that the bundled app starts locally from the new directory and document the exact run command in the bundle root `README.md`.
