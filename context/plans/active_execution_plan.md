# UI Test Execution Milestones

## Goal

Build a local-only CFDE-REVEAL v1.0 prototype in `ui_test/` that recommends reusable datasets from GTEx, MoTrPAC, and GEO based on an uploaded or pasted gene list, with provenance and evidence exposed per card.

## Scope constraints

- Work entirely inside `ui_test/`.
- Keep dataset recommendations and reuse rationale as the primary UI.
- Use GTEx, MoTrPAC, and GEO as the initial three sources.
- Prefer overlay-first provenance; patch `dig-gene-set-extractors` only if required.
- Target 12 to 20 cards for the first demo corpus.
- Build the UI in Streamlit.

## Milestones

### Milestone 1: Reconnaissance and output contract

Objective:
Establish how to reuse `dig-gene-set-extractors` without rewriting extraction logic.

Work:
- Inspect available extractor families and the RNA-oriented workflow paths.
- Identify the standard output contract for signatures, metadata, and provenance.
- Confirm whether provenance overlay support already exists and how it is emitted.
- Note any Harmonizome-related modes as comparison baselines only.

Deliverables:
- `ui_test/prototype/notes/repo_recon.md`
- initial assumptions for extractor invocation and output normalization

Exit criteria:
- At least one viable extractor workflow is identified for each planned source family.
- The required per-card fields are mapped to either existing outputs or overlay fields.

### Milestone 2: Source catalog and corpus selection

Objective:
Define the prototype corpus before implementation work expands.

Work:
- Create a resource catalog for GTEx, MoTrPAC, and GEO.
- Enumerate candidate dataset units and their access paths.
- Select a first demo corpus totaling 12 to 20 cards.
- Record friction points, especially around MoTrPAC processed outputs and GEO study cleanliness.

Deliverables:
- `ui_test/prototype/resource_catalog.yaml`
- `ui_test/prototype/notes/demo_corpus_selection.md`

Exit criteria:
- The selected corpus includes GTEx, MoTrPAC, and GEO.
- Every selected unit has a stable landing page, a clear contrast or tissue label, and a plausible extraction path.

### Milestone 3: Thin-slice extraction pipeline

Objective:
Produce a small but real set of signatures and provenance outputs end to end.

Work:
- Stand up the minimal directory layout under `ui_test/prototype/`.
- Run extraction for a thin slice from each source, prioritizing RNA-only paths.
- Emit standard extractor outputs for each run.
- Add provenance overlay JSON and normalize one manifest row per dataset unit.

Deliverables:
- `ui_test/prototype/signatures/`
- `ui_test/prototype/card_manifest.parquet`
- `ui_test/prototype/signature_index.parquet`
- `ui_test/prototype/provenance/`

Exit criteria:
- At least one successful dataset-unit extraction exists for GTEx, MoTrPAC, and GEO.
- Each successful unit has a signature file, metadata, provenance output, and a normalized card row.

### Milestone 4: Provenance audit and minimal patch decision

Objective:
Decide whether the current provenance model is sufficient for the UI.

Work:
- Audit whether extractor provenance plus overlay can recover all required UI fields.
- Validate landing page, access route, organism, comparison-space organism, source files, command, and publication identifiers.
- If the contract is incomplete, define the smallest backward-compatible patch.

Deliverables:
- `ui_test/prototype/notes/provenance_audit.md`
- optional `ui_test/prototype/patches/dig_gse_provenance.patch`
- optional `ui_test/prototype/patches/pr_summary.md`

Exit criteria:
- Either overlay-first provenance is documented as sufficient, or a minimal patch plan is written with exact missing fields and compatibility notes.

### Milestone 5: Local retrieval index and scoring

Objective:
Turn the extracted corpus into a searchable recommendation backend.

Work:
- Build a local DuckDB or SQLite index for cards, signatures, and provenance pointers.
- Implement gene normalization and comparison-space handling.
- Compute overlap-based, enrichment-style, and latent signature similarity scores.
- Combine scores into a final ranking contract for the UI.

Deliverables:
- `ui_test/prototype/app/data/` index artifacts
- `ui_test/prototype/app/src/retrieval.py`
- `ui_test/prototype/notes/scoring_contract.md`

Exit criteria:
- A pasted gene list can be scored against all prototype cards in one local command.
- Every result includes enough structured evidence for a card view and an evidence drawer.

### Milestone 6: Streamlit prototype UI

Objective:
Expose the prototype as a local app centered on dataset recommendations.

Work:
- Build the landing/input, ranked results, evidence drawer, and provenance viewer areas.
- Keep cards as the default ranked-results presentation.
- Use template-based rationale text driven by metadata and evidence.
- Support pasted gene lists, uploaded text files, and core filters.

Deliverables:
- `ui_test/prototype/app/`
- `ui_test/prototype/README.md`

Exit criteria:
- The app starts with one documented command.
- Results render as recommendation cards with evidence and provenance available per card.
- The top-level language emphasizes relevant datasets and reuse value rather than generic enrichment hits.

### Milestone 7: Corpus completion and UX hardening

Objective:
Expand from the thin slice to a convincing demo corpus and stabilize the UI.

Work:
- Fill the corpus to 12 to 20 cards.
- Refine card copy, filters, and provenance summaries.
- Add per-card downloads for signature TSV, provenance JSON, and card JSON.
- Capture screenshots of the main happy-path flows.

Deliverables:
- completed `ui_test/prototype/card_manifest.parquet`
- `ui_test/prototype/screenshots/`

Exit criteria:
- The corpus reaches the target size with all three sources represented.
- Each card has a working landing page, provenance JSON, and reuse-focused rationale.

### Milestone 8: Smoke-test validation and handoff

Objective:
Verify the prototype meets the stated acceptance criteria.

Work:
- Run muscle, immune/inflammation, and mitochondrial smoke tests.
- Confirm cross-resource plausibility in the returned cards.
- Verify startup steps, download behavior, and provenance visibility.
- Document known gaps and deferred work.

Deliverables:
- `ui_test/prototype/eval/smoke_tests.md`
- final updates to `ui_test/prototype/README.md`

Exit criteria:
- The prototype satisfies the three-source, 12 to 20 card, local-run, evidence-backed, provenance-aware acceptance criteria.
- Remaining limitations are documented explicitly enough for the next implementation pass.

## Immediate next step

If the next iteration proceeds, focus on two linked tracks:

1. provenance graph patching
2. dataset-first result-card and provenance-viewer redesign

See `ui_test/plans/provenance_graph_and_dataset_first_ui_plan.md` for the proposed scope and sequencing.

Implementation note:
- The embedded `dig-gene-set-extractors` repo now emits extractor-native `lineage` metadata with explicit file nodes, process-step edges, and invocation records. The next prototype iteration should consume that structure instead of inferring the conversion step from overlay text.

Current implementation sequence:
1. regenerate prototype signatures against the patched extractor with explicit RNA legacy-mode settings
2. build normalized per-card provenance graph JSON by combining source preparation files with extractor-native lineage
3. refactor the Streamlit app so result cards foreground dataset reuse and the provenance panel renders the normalized graph plus node/edge detail panes
4. align the Streamlit prototype visual language with the `cfde-main` DIG/DUG mechanism discovery component so later porting preserves header, palette, and section-chrome expectations

Current status:
- Completed the corpus rebuild against the patched extractor and regenerated overlay metadata.
- Completed normalized `prototype/provenance_graph/*.graph.json` generation for all demo cards.
- Completed a dataset-first Streamlit redesign with interactive file-node graph rendering.
- Completed visual alignment of the Streamlit prototype with `broadinstitute/dig-dug-portal` `cfde-main` `cfdeMechanismDiscovery.vue`.
- Remaining refinement: edge details are available through the processing-step selector rather than direct edge-click capture from the graph component.

## Additional Extraction Track

### Goal

Expand the demo beyond the original GTEx, MoTrPAC RNA, and GEO thin slice so the UI exposes at least one real run for more extractor families and more CFDE DCCs, following `ui_test/plans/additional_extractions.txt`.

### Constraints

- Keep all new prototype orchestration inside `ui_test/`.
- Preserve the current UI-facing artifact contract:
  - `card_manifest.parquet`
  - `signature_index.parquet`
  - `provenance/*.overlay.json`
  - `provenance_graph/*.graph.json`
- Add the requested run bookkeeping under:
  - `ui_test/prototype_inputs/`
  - `ui_test/prototype_outputs/`
  - `ui_test/prototype_runs/run_log.tsv`
- Prefer the smallest generalization of the current scripts over a full rewrite.

### Execution sequence

1. Generalize the prototype build scripts so new dataset slugs can be staged in `prototype_inputs/`, extracted into `prototype_outputs/`, and registered in a shared run log.
2. Land the highest-confidence new public runs first:
   - GTEx grouped bulk RNA
   - MoTrPAC proteomics
   - MoTrPAC PTM
   - 4DN ATAC
3. If public processed inputs are straightforward, extend to:
   - 4DN ChIP-seq
   - one public single-cell run from HuBMAP, SPARC, Bridge2AI, or SenNet
4. Rebuild the manifest/index/UI over both the legacy prototype corpus and the newly registered runs.
5. Record blocked datasets explicitly when the public file path, preprocessing burden, or extractor-input contract is not yet clean enough for a same-turn implementation.

### Success criteria

- At least three new public dataset runs are available through the Streamlit UI.
- Every new run has:
  - extractor-standard outputs
  - provenance overlay JSON
  - a run-log entry
  - a normalized card row in the rebuilt manifest
- Grouped outputs retain enriched `manifest.tsv` columns and surface individual cards in the UI.

### Current status

- Completed the shared additional-run layout under `ui_test/prototype_inputs/`, `ui_test/prototype_outputs/`, and `ui_test/prototype_runs/`.
- Completed the first successful additional public runs:
  - `motrpac_wat_proteomics_diff`
  - `motrpac_wat_ptm`
  - `4dn_atac_bulk_4DNFIUALWN8X_pvalb`
- Completed the second successful additional public tranche:
  - `gtex_bulk_rna_grouped` with five GTEx tissue-versus-rest cards
  - `tcga_brca_splice_event_matrix`
  - `depmap_prism_quick_grouped` with selected lung, pancreas, and breast cards
  - `tcga_brca_cnv_segments` with representative amplification and deletion cards
- Completed the third successful additional public tranche:
  - `jump_cellpainting_u2os_morphology` with four selected held-out query cards
  - `calr_fig2_ontology_grouped` with selected HFD and LFD global cards
  - `calr_fig2_profile_query_public_grouped` with selected HFD and LFD global cards
  - `gse42752_methylation_cpg_diff`
- Completed a fourth additional public run:
  - `4dn_chipseq_peak_4DNFIGINV1VI_h3k27ac` via `chipseq_peak`
- Completed a fifth additional public run:
  - `hubmap_liver_scrna_markers` via `sc_rna_marker`, using the published HuBMAP liver processed RNA-seq data product and exported Leiden groups
- Completed manifest/index integration for the expanded tranche so the UI corpus now includes 48 cards total.
- Remaining work after this tranche:
  - final inventory/reporting polish
  - explicit blocked-dataset notes for SPARC, Bridge2AI/CM4AI, and SenNet unless a smaller direct-download public processed matrix is identified
- Current blocker assessment for the remaining CFDE DCC sections:
  - SPARC: the best Pennsieve-hosted scRNA candidate discovered so far exposes public metadata and file listings, but the dataset is over 5 GB and Pennsieve documents AWS requester-pays download as the supported path for larger public datasets
  - Bridge2AI / CM4AI: the current public Dataverse beta release exposes release metadata and provenance sidecars, but not a directly downloadable processed perturb-seq matrix suitable for same-turn extractor staging
  - SenNet: the public portal exposes ontology and entity helpers, but this pass did not uncover a stable unauthenticated processed scRNA/scATAC download path that could be scripted into the extractor pipeline
