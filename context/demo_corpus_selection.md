# Milestone 2: Demo Corpus Selection

## Target size

The first demo corpus is fixed at 18 cards:

- 6 GTEx tissue cards
- 6 MoTrPAC transcriptomics contrast cards
- 6 GEO Kang IFN-beta contrast cards

This stays inside the requested 12 to 20 card range while covering all three required sources.

## Selection rules applied

Each selected unit satisfies the Milestone 2 filter:

- stable public landing page
- clear biological label
- clear assay or modality
- plausible extractor path
- provenance can name source files and extraction method

## Selected GTEx cards

Chosen for broad biological coverage and clean tissue interpretation:

- `GTEx_V8__Muscle_Skeletal__tissue_specific_expression`
- `GTEx_V8__Liver__tissue_specific_expression`
- `GTEx_V8__Whole_Blood__tissue_specific_expression`
- `GTEx_V8__Lung__tissue_specific_expression`
- `GTEx_V8__Pancreas__tissue_specific_expression`
- `GTEx_V8__Heart_Left_Ventricle__tissue_specific_expression`

Why these:

- They cover metabolic, immune, respiratory, and muscle systems.
- They are easy to explain in card language.
- They create plausible wins for the planned smoke tests, especially muscle and mitochondrial lists.

## Selected MoTrPAC cards

Chosen to keep the first slice RNA-only while retaining clear training-vs-control semantics:

- `MoTrPAC__SKMGN__RNA__male_8w_vs_control`
- `MoTrPAC__SKMGN__RNA__female_8w_vs_control`
- `MoTrPAC__LIVER__RNA__male_8w_vs_control`
- `MoTrPAC__LIVER__RNA__female_8w_vs_control`
- `MoTrPAC__HEART__RNA__male_8w_vs_control`
- `MoTrPAC__HEART__RNA__female_8w_vs_control`

Why these:

- They align with the prototype narrative around exercise-trained tissues.
- They reuse a single public data package with a consistent table format.
- They create immediate complementarity to GTEx tissue cards.

Notes:

- `SKMGN` is gastrocnemius skeletal muscle in the MoTrPAC package documentation.
- These cards require an explicit cross-species overlay because the source identifiers are rat.

## Selected GEO cards

The first GEO slice uses one well-understood study with multiple contrasts:

- `GEO_GSE96583__B_cells__ifnb_stim_vs_ctrl`
- `GEO_GSE96583__CD14_Monocytes__ifnb_stim_vs_ctrl`
- `GEO_GSE96583__CD4_T_cells__ifnb_stim_vs_ctrl`
- `GEO_GSE96583__CD8_T_cells__ifnb_stim_vs_ctrl`
- `GEO_GSE96583__Dendritic_cells__ifnb_stim_vs_ctrl`
- `GEO_GSE96583__FCGR3A_Monocytes__ifnb_stim_vs_ctrl`

Why this study:

- A local workflow and extracted outputs already exist in this workspace.
- The biological contrasts are obvious and useful for immune-focused queries.
- It reduces acquisition risk while still using a real GEO source.

Tradeoff:

- This is single-cell PBMC data rather than bulk tissue. The UI and provenance summary must name the cell-type focus explicitly.

## Deferred candidates

Deferred for later corpus expansion:

- GTEx adipose and pancreas-adjacent tissues
- MoTrPAC blood, lung, and kidney transcriptomics contrasts
- GEO Kang `NK_cells_stim_vs_ctrl`
- GEO Kang `Megakaryocytes_stim_vs_ctrl`
- additional human bulk GEO case/control series once the first thin slice is working

## Milestone 2 decision

Milestone 2 is complete enough to proceed. The initial corpus is fixed and scriptable, and no source requires browser automation or new credentials.
