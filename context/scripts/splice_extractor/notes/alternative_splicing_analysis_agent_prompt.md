# Analysis prompt: validate and dissect the alternative splicing extractor implementation in dig-gene-set-extractors

You are reviewing a branch or local checkout of `flannick/dig-gene-set-extractors` that should contain a new splicing assay family.

Your job is to:
1. clone and install the repo
2. run the test suite, especially the new splicing tests
3. run the splicing converters and workflows on included toy fixtures
4. run the splicing converters and workflows on at least one real or public dataset if available
5. run `pigean` and `eaggl` on the resulting outputs
6. summarize the top genes, gene sets, and factors
7. judge whether the outputs reflect plausible biology or likely artifacts
8. write a concise but evidence-rich report

Do not invent successful runs. If a command fails or a tool is missing, record that clearly and continue with the rest of the review.

## Assumptions

- The repo is `https://github.com/flannick/dig-gene-set-extractors`
- The splicing implementation includes:
  - `splice_event_diff`
  - `splice_event_matrix`
  - `workflows splice_prepare_public`
  - `workflows splice_prepare_reference_bundle`
- The repo includes splicing docs under:
  - `docs/assays/splicing/guide.md`
  - `docs/assays/splicing/reference_bundle.md`
  - `docs/assays/splicing/methods.tex`
- The repo includes targeted tests:
  - `tests/test_splice_event_diff_converter.py`
  - `tests/test_splice_event_matrix_converter.py`
  - `tests/test_splice_prepare_public_workflow.py`
  - `tests/test_splice_prepare_reference_bundle_workflow.py`
- The repo includes splicing toy fixtures under `tests/data`, including:
  - `toy_splice_event_diff.tsv`
  - `toy_splice_matrix.tsv`
  - `toy_splice_sample_metadata.tsv`
  - `toy_splice_coverage_matrix.tsv`
  - `toy_spliceseq_public.tsv`
  - `toy_spliceseq_sample_annotations.tsv`
  - packaged tiny bundle tables such as `splice_event_aliases_human_v1.tsv.gz`
- `pigean` and `eaggl` may or may not already be installed in the environment

## Step 0: choose datasets before you run anything

Always do both of these layers:

### A. Regression layer: must run on repo fixtures

These are not biological truth sets. They are parser, scoring, metadata, and workflow regression checks.

Run at minimum:
- the direct differential converter on `tests/data/toy_splice_event_diff.tsv`
- the matrix converter on `tests/data/toy_splice_matrix.tsv` plus `tests/data/toy_splice_sample_metadata.tsv`
- the public staging workflow on `tests/data/toy_spliceseq_public.tsv` plus `tests/data/toy_spliceseq_sample_annotations.tsv`
- the reference-bundle workflow using the staged `bundle_source_row.tsv`
- one bundle-backed run and one no-bundle run

Expected regression outcomes from the current toy fixtures:
- `toy_splice_event_diff.tsv`
  - should validate successfully
  - should emit signed GMT outputs with both positive and negative sets
  - with the packaged toy bundle, the top-ranked gene should include `KCNN4`
  - metadata and run summaries should record all three bundle resources when they are used
  - changing `impact_mode` or `event_dup_policy` should materially change at least one gene score, especially `G_KCNN4`
  - changing `ambiguous_gene_policy` should control whether ambiguous toy genes such as `GENEA` and `GENEB` appear
- `toy_splice_matrix.tsv`
  - single-contrast mode should emit exactly one program
  - grouped mode such as `condition_within_group` should emit two child outputs plus `manifest.tsv`
  - top genes should include at least one of `KCNN4`, `MAPK1`, or `GENEA`
  - `welch_t` versus `mean_diff` should change at least one score, especially `G_KCNN4`
  - missingness filters should drop some events when `min_present_per_condition` is strict
- `toy_spliceseq_public.tsv` plus `toy_spliceseq_sample_annotations.tsv`
  - public staging should report 4 samples and 3 events
  - staged output should include `bundle_source_row.tsv`
  - `bundle_source_row.tsv` should have 12 rows in the toy case
  - a bundle built from the staged toy source should report 3 canonical events
  - the bundle should auto-resolve at runtime when `--resources_dir` points at it
- warning-path fixtures dominated by retained introns or low-confidence events should trigger warnings in stderr and in `run_summary.json`

Treat any deviation from those expectations as a possible implementation issue unless the branch intentionally changed the test contract.

### B. Biology layer: run at least one real or public dataset if feasible

Priority order:
1. real user-provided compatible splicing inputs
2. public TCGA SpliceSeq tumor-vs-normal data
3. already prepared LeafCutter, MAJIQ, or Whippet differential tables
4. if none of the above is available, stop after the regression layer and ask the user for inputs using the exact request block below

Do not download raw FASTQ and rerun splicing callers from scratch unless the user explicitly asks for that. This review is about the extractor layer, so prefer already prepared PSI tables, delta-PSI tables, or differential event tables.

### Compatible real-input families

If the user gave you inputs, map them as follows:
- `tcga_spliceseq`: PSI table plus sample annotations -> `workflows splice_prepare_public`, then `splice_event_matrix`, and optionally `splice_prepare_reference_bundle`
- `leafcutter`: event or effect table, optionally cluster stats -> `splice_event_diff --tool_family leafcutter`
- `majiq`: differential event table -> `splice_event_diff --tool_family majiq`
- `whippet`: differential event table -> `splice_event_diff --tool_family whippet`
- `generic`: event-diff TSV or PSI matrix plus sample metadata -> `splice_event_diff` or `splice_event_matrix`

### Public defaults if the user did not provide real inputs

Prefer `TCGA SpliceSeq` first because the repo has a dedicated `splice_prepare_public --input_mode tcga_spliceseq` workflow for that family.

#### Exact public-data acquisition recipe for TCGA SpliceSeq

Use this recipe first unless the user already gave you compatible real inputs.

##### A. Current PSI download page

Use the current MD Anderson TCGA SpliceSeq page:

- `https://bioinformatics.mdanderson.org/TCGASpliceSeq/PSIdownload.jsp`

Do not use the older `projects.insilico.us.com` URL commonly cited in older papers. That older URL is reported broken in the SpliceSeq GitHub issue tracker, while the current MD Anderson public-software page points to the `bioinformatics.mdanderson.org` site.

Important operational notes:
- treat the PSI acquisition as a browser-form download, not as a documented API call
- I could verify the current form URL and the fields shown on the page, but I could not verify a stable direct scripted POST or GET endpoint for automating the PSI download itself
- do not spend large amounts of review time reverse-engineering the form if it is not trivial in your environment
- the page says the result is a tab-delimited event-by-sample matrix and warns that large whole-cohort downloads such as all-BRCA may time out
- the older FAQ text says the returned file may be comma-separated, so always inspect delimiter and normalize to TSV before using the repo workflow
- if a whole-cohort download times out or returns an empty or partial file, retry by downloading one splice-event type at a time and merge the per-type files after normalizing delimiter

##### B. Exact form settings for a compatible whole-cohort PSI matrix

Use one cohort at a time. Start with `BRCA`, then try `LUAD`, `PRAD`, `KIRC`, and `THCA` if needed.

For the MD Anderson form, use:
- TCGA disease type: the cohort code, for example `BRCA`
- Gene HUGO symbols: leave blank
- TCGA Samples: leave blank
- Splice event types: select all seven types `AA`, `AD`, `ES`, `RI`, `AP`, `AT`, `ME`
- Percentage of Samples with PSI Value: `10`
- Minimum Average Expression Percentage: `0`
- Minimum PSI Range (delta across samples): `0`
- Minimum PSI Standard Deviation: `0`
- Include a Gene RPKM Matrix: `No`
- Include Gene Descriptions: `No`
- Include Clinical Data: `No`

Save the raw download as:
- `data/public/tcga_spliceseq/<COHORT>/<COHORT>_spliceseq_raw.txt`

If the whole-cohort download fails, download one event type at a time with the same settings and save them as:
- `data/public/tcga_spliceseq/<COHORT>/<COHORT>_AA_raw.txt`
- `data/public/tcga_spliceseq/<COHORT>/<COHORT>_AD_raw.txt`
- `data/public/tcga_spliceseq/<COHORT>/<COHORT>_ES_raw.txt`
- `data/public/tcga_spliceseq/<COHORT>/<COHORT>_RI_raw.txt`
- `data/public/tcga_spliceseq/<COHORT>/<COHORT>_AP_raw.txt`
- `data/public/tcga_spliceseq/<COHORT>/<COHORT>_AT_raw.txt`
- `data/public/tcga_spliceseq/<COHORT>/<COHORT>_ME_raw.txt`

Normalize each raw file to TSV before you do anything else.

Exact delimiter-normalization command:

```bash
RAW=data/public/tcga_spliceseq/BRCA/BRCA_spliceseq_raw.txt
NORM=data/public/tcga_spliceseq/BRCA/BRCA_spliceseq.tsv
python - "$RAW" "$NORM" <<'PY'
import csv
import sys
inp, out = sys.argv[1], sys.argv[2]
with open(inp, 'r', newline='') as f:
    sample = f.read(65536)
dialect = csv.Sniffer().sniff(sample, delimiters=',\t')
with open(inp, 'r', newline='') as fin, open(out, 'w', newline='') as fout:
    reader = csv.reader(fin, dialect)
    writer = csv.writer(fout, delimiter='\t', lineterminator='\n')
    for row in reader:
        writer.writerow(row)
PY
```

Exact merge command for per-event-type downloads:

```bash
python - \
  data/public/tcga_spliceseq/BRCA/BRCA_spliceseq.tsv \
  data/public/tcga_spliceseq/BRCA/BRCA_AA.tsv \
  data/public/tcga_spliceseq/BRCA/BRCA_AD.tsv \
  data/public/tcga_spliceseq/BRCA/BRCA_ES.tsv \
  data/public/tcga_spliceseq/BRCA/BRCA_RI.tsv \
  data/public/tcga_spliceseq/BRCA/BRCA_AP.tsv \
  data/public/tcga_spliceseq/BRCA/BRCA_AT.tsv \
  data/public/tcga_spliceseq/BRCA/BRCA_ME.tsv <<'PY'
import csv
import sys
out = sys.argv[1]
inputs = sys.argv[2:]
header = None
with open(out, 'w', newline='') as fout:
    writer = None
    for path in inputs:
        with open(path, 'r', newline='') as fin:
            reader = csv.reader(fin, delimiter='\t')
            this_header = next(reader)
            if header is None:
                header = this_header
                writer = csv.writer(fout, delimiter='\t', lineterminator='\n')
                writer.writerow(header)
            elif this_header != header:
                raise SystemExit('Header mismatch: %s' % path)
            for row in reader:
                writer.writerow(row)
PY
```

##### C. Exact sample-annotation source and URLs

Use UCSC Xena clinical matrices as the sample-annotation source because the cohort metadata explicitly lists these datasets and exposes `sample_type` for the relevant TCGA cohorts.

Cohort-specific sample-annotation URLs:
- BRCA: `https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix.gz`
- LUAD: `https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix.gz`
- PRAD: `https://tcga.xenahubs.net/download/TCGA.PRAD.sampleMap/PRAD_clinicalMatrix.gz`
- KIRC: `https://tcga.xenahubs.net/download/TCGA.KIRC.sampleMap/KIRC_clinicalMatrix.gz`
- THCA: `https://tcga.xenahubs.net/download/TCGA.THCA.sampleMap/THCA_clinicalMatrix.gz`

Exact download pattern:

```bash
COHORT=BRCA
ROOT="data/public/tcga_spliceseq/${COHORT}"
mkdir -p "$ROOT"
curl -L -A "Mozilla/5.0" \
  "https://tcga.xenahubs.net/download/TCGA.${COHORT}.sampleMap/${COHORT}_clinicalMatrix.gz" \
  -o "$ROOT/${COHORT}_clinicalMatrix.tsv.gz"
```

If `curl` gets blocked in your environment, open the same URL in a browser and save it manually. Keep the same filename.

##### D. Exact tumor-vs-adjacent-normal selection rules

Use these rules for the primary tumor-versus-adjacent-normal contrast:
- keep only one TCGA cohort per run
- retain only samples that appear in both the PSI matrix header and the cohort-matched clinical matrix
- keep `tumor` samples when `sample_type` is `Primary Tumor` or `Primary Solid Tumor`, or when the TCGA barcode sample code is `01`
- keep `adjacent_normal` samples when `sample_type` is `Solid Tissue Normal`, or when the TCGA barcode sample code is `11`
- drop every other sample type, including recurrent tumor, metastatic tumor, blood-derived normal, blood-derived cancer, cell lines, xenografts, and unlabeled samples
- prefer sample-level matching over patient-level matching
- normalize IDs in this order and pick the first unambiguous match mode that works against the PSI header: full barcode, first 15 characters, first 12 characters
- if matching only works at the first-12-character patient level and any patient maps to both tumor and adjacent-normal samples, treat that cohort as ambiguous and skip it rather than guessing
- if a cohort has zero matched adjacent-normal samples after filtering, skip that cohort and try the next one
- for the main biology demonstration, prefer cohorts with at least 20 tumor and 5 adjacent-normal samples after matching; smaller cohorts can still be used as parser checks, but label them as low power

##### E. Exact command to build a minimal sample-annotation TSV matched to the PSI header

This command downloads no data on its own. It takes the already downloaded cohort PSI TSV and Xena clinical matrix, infers the matching ID mode, applies the tumor-versus-adjacent-normal rules above, and writes a minimal sample-annotation table suitable for the extractor workflow.

```bash
COHORT=BRCA
ROOT="data/public/tcga_spliceseq/${COHORT}"
PSI="$ROOT/${COHORT}_spliceseq.tsv"
CLIN="$ROOT/${COHORT}_clinicalMatrix.tsv.gz"
ANN="$ROOT/${COHORT}_sample_annotations.tsv"
python - "$COHORT" "$PSI" "$CLIN" "$ANN" <<'PY'
import csv
import gzip
import sys
cohort, psi_path, clin_path, out_path = sys.argv[1:5]

def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', newline='')
    return open(path, 'r', newline='')

with open_text(psi_path) as f:
    header = next(csv.reader(f, delimiter='\t'))
psi_samples = [x.strip() for x in header if x.strip().startswith('TCGA-')]
if not psi_samples:
    raise SystemExit('No TCGA-like sample columns found in PSI header')

clinical_rows = []
with open_text(clin_path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        sid = (row.get('sampleID') or row.get('sample') or row.get('Sample') or '').strip()
        st = (row.get('sample_type') or row.get('_sample_type') or '').strip()
        if sid and st:
            clinical_rows.append((sid, st))
if not clinical_rows:
    raise SystemExit('Could not parse sampleID and sample_type from clinical matrix')

def code_from_sid(sid):
    return sid[13:15] if len(sid) >= 15 else ''

candidates = [
    ('full', lambda s: s),
    ('sample15', lambda s: s[:15]),
    ('patient12', lambda s: s[:12]),
]

best_name = None
best_matches = -1
for name, fn in candidates:
    keyset = {fn(sid) for sid, _ in clinical_rows}
    matches = sum(1 for s in psi_samples if s in keyset)
    if matches > best_matches:
        best_name = name
        best_matches = matches

if best_matches <= 0:
    raise SystemExit('No clinical IDs match PSI header')

fn = dict(candidates)[best_name]
sample_to_rows = {}
for sid, st in clinical_rows:
    key = fn(sid)
    if key not in psi_samples:
        continue
    code = code_from_sid(sid)
    if st in ('Primary Tumor', 'Primary Solid Tumor') or code == '01':
        condition = 'tumor'
    elif st in ('Solid Tissue Normal',) or code == '11':
        condition = 'adjacent_normal'
    else:
        continue
    sample_to_rows.setdefault(key, []).append((sid, st, condition))

ambiguous = [k for k, v in sample_to_rows.items() if len({x[2] for x in v}) > 1]
if ambiguous and best_name == 'patient12':
    raise SystemExit(
        'Ambiguous patient-level matching: %d patient IDs map to both tumor and normal. '
        'Skip this cohort or use a PSI table with sample-level IDs.' % len(ambiguous)
    )

with open(out_path, 'w', newline='') as out:
    writer = csv.writer(out, delimiter='\t', lineterminator='\n')
    writer.writerow(['sample_id', 'patient_id', 'sample_type', 'condition', 'cohort', 'id_mode', 'barcode_full'])
    written = 0
    for key in sorted(sample_to_rows):
        sid, st, condition = sorted(sample_to_rows[key])[0]
        writer.writerow([key, sid[:12], st, condition, cohort, best_name, sid])
        written += 1

tumor = 0
normal = 0
with open(out_path, 'r', newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if row['condition'] == 'tumor':
            tumor += 1
        elif row['condition'] == 'adjacent_normal':
            normal += 1
print({'cohort': cohort, 'id_mode': best_name, 'matched_samples': written, 'tumor': tumor, 'adjacent_normal': normal})
PY
```

##### F. Exact workflow command once the files exist

After you have the normalized PSI TSV and the matched sample-annotation TSV, run the repo workflow like this:

```bash
COHORT=BRCA
geneset-extractors workflows splice_prepare_public \
  --input_mode tcga_spliceseq \
  --psi_tsv data/public/tcga_spliceseq/${COHORT}/${COHORT}_spliceseq.tsv \
  --sample_annotations_tsv data/public/tcga_spliceseq/${COHORT}/${COHORT}_sample_annotations.tsv \
  --out_dir tests/tmp/splice_prepare_public_${COHORT} \
  --organism human \
  --genome_build hg38 \
  --study_id TCGA_${COHORT} \
  --study_label tcga_spliceseq_${COHORT}
```

Record in your report:
- the exact MD Anderson page URL you used for PSI download
- whether the PSI acquisition needed split-by-event-type downloads
- the Xena clinicalMatrix URL used for sample annotations
- the inferred ID matching mode: `full`, `sample15`, or `patient12`
- the final tumor and adjacent-normal sample counts
- whether any cohort had to be skipped because of missing normals, ambiguous ID matching, or broken downloads

##### G. Practical failure rules

If the current environment cannot obtain the PSI table from the MD Anderson form after one whole-cohort attempt and, if needed, one split-by-event-type attempt, do not burn the rest of the review time on scraping. Instead:
- finish the regression-layer review on repo fixtures
- record the exact failure mode for the public download attempt
- if no compatible public PSI table is available locally, then send the exact input request block from this prompt to ask the user for real inputs

Remember the TCGA SpliceSeq FAQ warning: adjacent normal tissues are imperfect controls because of tissue-composition differences and possible field effects. Use them for extractor validation and broad biology checks, but do not overstate tumor-versus-normal results as purely tumor-intrinsic splicing biology.

Try at least one of the following tumor-vs-adjacent-normal settings once you have both the PSI table and matching sample annotations:

- `BRCA` tumor vs adjacent normal
  - Why this is a good default: common carcinoma, often gives a clear epithelial versus mesenchymal alternative-splicing axis
  - What you should hope to see: RNA-splicing or RNA-processing signal, epithelial differentiation, cell-cell junction or cytoskeletal programs, and ESRP-linked epithelial-splicing themes
  - Plausible top genes or event-target families: `CD44`, `CTNND1`, `ENAH`, `FGFR2`, or related epithelial-program genes
  - Common artifacts: subtype mixing, tumor-purity differences, and stromal admixture masquerading as EMT

- `LUAD` tumor vs adjacent normal
  - Why this is useful: carcinoma with frequent splicing dysregulation and a wide range of tumor-state variation
  - What you should hope to see: RNA-processing or spliceosome signal, proliferation or cell-cycle signal, and epithelial-versus-mesenchymal state changes
  - Plausible special case: if the chosen cohort or subset is enriched for `RBM10` disruption, you may see exon-skipping style signal or translation-related downstream programs
  - Common artifacts: smoking-related RNA quality shifts, retained intron inflation, and tumor-purity effects

- `PRAD` tumor vs adjacent normal
  - Why this is useful: often gives a strong epithelial and secretory-state signal and is a good place to see ESRP-linked epithelial splicing
  - What you should hope to see: epithelial, junction, cell-adhesion, secretory, or androgen-linked programs more than generic immune noise
  - Common artifacts: stromal contamination and broad tissue-composition effects

- `KIRC` or `THCA` tumor vs adjacent normal
  - Why this is useful: good stress test for whether the extractor overcalls broad tumor-state shifts as splicing-specific biology
  - What you should hope to see: broad tumor-normal separation, tissue identity, metabolism, or hypoxia-like signal rather than a very crisp canonical splicing-factor program
  - Interpretation rule: if the signal is broad and tissue-state driven, do not oversell it as spliceosome-specific biology

### Tool-family specific public data, if already available as tables

If you have prepared `LeafCutter`, `MAJIQ`, or `Whippet` differential tables, prioritize contrasts with known splicing perturbation instead of random case-control comparisons.

Best positive-control contrast types:
- splicing-factor mutant versus wild-type, such as `SF3B1`, `U2AF1`, `SRSF2`, or `RBM10`
- EMT or epithelial-mesenchymal state contrasts
- regulator perturbations involving `ESRP1` or `ESRP2`

What to expect at the gene-set and factor level:
- `SF3B1`: alternative 3-prime splice-site style biology and broad RNA-processing or spliceosome signal
- `RBM10`: exon-skipping or cassette-exon programs plus coherent downstream proliferation or translation effects
- `ESRP1` or `ESRP2`, or EMT contrasts: epithelial-versus-mesenchymal splicing programs, cell junction and cytoskeletal gene sets, and genes such as `CD44`, `ENAH`, `FGFR2`, or `CTNND1`
- `U2AF1` or `SRSF2`: splice-site choice or RNA-processing signal rather than generic DE-like pathway output

If you need a non-cancer parser sanity control because you cannot find processed cancer tables for a given tool family, a strong tissue-separation splicing dataset is acceptable. A heart-versus-brain or epithelial-versus-mesenchymal prepared table is a better parser check than a weak cancer contrast. Make it explicit that this is a parser or scoring sanity check, not the main cancer-biology demonstration.

### If no compatible real or public tables are available

Stop after the toy-fixture regression layer and send the user exactly this:

```text
Send one of these:

  - the paths to your real splicing inputs
  - a dataset name/accession you want me to target
  - the tool family and file types you have, for example:
      - tcga_spliceseq: PSI table + sample annotations
      - leafcutter: event/effect table, optionally cluster stats
      - majiq: differential event table
      - whippet: differential event table
      - generic: event diff TSV or PSI matrix + sample metadata
```

## Step 1: clone and install

Run something like:

```bash
git clone https://github.com/flannick/dig-gene-set-extractors.git
cd dig-gene-set-extractors
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

Record:
- commit SHA
- branch name
- python version

## Step 2: repo-level smoke checks

Run:

```bash
geneset-extractors list
pytest -q
```

If the full suite is too slow, still run the targeted splicing tests below.

Record:
- whether the new converters appear in `geneset-extractors list`
- whether the new docs seem present
- whether the full suite passes

## Step 3: targeted splicing tests

Run at least:

```bash
pytest -q tests/test_splice_event_diff_converter.py
pytest -q tests/test_splice_event_matrix_converter.py
pytest -q tests/test_splice_prepare_public_workflow.py
pytest -q tests/test_splice_prepare_reference_bundle_workflow.py
```

If a test file is missing, record that explicitly.

For each failing test:
- capture the exact test name
- capture the traceback or core error
- state whether the failure looks like parsing, scoring, resource handling, grouped output, metadata, workflow handoff, or test-fixture drift

## Step 4: inspect CLI and docs

Run:

```bash
geneset-extractors convert splice_event_diff --help
geneset-extractors convert splice_event_matrix --help
geneset-extractors workflows splice_prepare_public --help
geneset-extractors workflows splice_prepare_reference_bundle --help
```

Inspect:
- `docs/assays/splicing/guide.md`
- `docs/assays/splicing/reference_bundle.md`
- `docs/assays/splicing/methods.tex`

Summarize:
- whether the CLI matches the documentation
- whether the docs explain missing-bundle behavior, defaults, and input contracts
- whether the methods note matches the implemented scoring choices
- whether the docs tell users what types of datasets should produce strong biology versus weak or artifact-prone output

## Step 5: run the public staging workflow on the included toy TCGA SpliceSeq-like data

Use the actual fixture paths from `tests/data`:

```bash
geneset-extractors workflows splice_prepare_public \
  --input_mode tcga_spliceseq \
  --psi_tsv tests/data/toy_spliceseq_public.tsv \
  --sample_annotations_tsv tests/data/toy_spliceseq_sample_annotations.tsv \
  --out_dir tests/tmp/splice_prepare_public \
  --organism human \
  --genome_build hg38 \
  --study_id TCGA_TOY \
  --study_label toy_tcga_splice
```

Then inspect:
- `psi_matrix.tsv`
- `sample_metadata.tsv`
- `event_metadata.tsv`
- `bundle_source_row.tsv`
- `prepare_summary.json`

Check:
- row and column counts
- whether the toy workflow reports 4 samples and 3 events
- whether `bundle_source_row.tsv` has 12 rows
- whether event ids are standardized
- how missing PSI is handled
- whether sample and event id maps look sensible

If the workflow supports it and a real or public TCGA SpliceSeq dataset is available, run the same workflow on one real tumor-vs-normal cohort as well.

For a real or public TCGA SpliceSeq cohort, record:
- tumor type
- how many tumor and normal samples were used
- whether sample annotations were complete or had to be repaired
- whether the prepared event metadata looks plausible

## Step 6: build the reference bundle from the staged source rows

Use the staged `bundle_source_row.tsv` to make a small local bundle.

Example pattern:

```bash
printf "path\tsource_dataset\n%s\tTCGA_TOY\n" \
  "$(pwd)/tests/tmp/splice_prepare_public/bundle_source_row.tsv" \
  > tests/tmp/splice_sources.tsv

geneset-extractors workflows splice_prepare_reference_bundle \
  --sources_tsv tests/tmp/splice_sources.tsv \
  --out_dir tests/tmp/splice_bundle \
  --organism human \
  --bundle_id toy_splice_bundle_v1
```

Inspect:
- `splice_event_aliases_human_v1.tsv.gz`
- `splice_event_ubiquity_human_v1.tsv.gz`
- `splice_event_impact_human_v1.tsv.gz`
- `bundle_provenance.json`
- `local_resources_manifest.json`

Check:
- whether the toy bundle reports 3 canonical events
- whether the bundle can be auto-resolved by `--resources_dir`
- whether the alias table has deterministic mappings
- whether ubiquity counts look numerically sensible
- whether impact priors stay conservative rather than extreme

If you also built a real or public bundle from one or more staged TCGA SpliceSeq cohorts, compare the toy bundle and real bundle behavior:
- toy bundle should mainly prove wiring and resource auto-resolution
- real bundle should mildly improve harmonization and specificity, not dominate the ranking

## Step 7: run `splice_event_diff` on included toy fixtures and any real differential tables

### 7A. Mandatory toy generic differential run

Run the direct converter on the included generic toy event-diff table.

```bash
geneset-extractors convert splice_event_diff \
  --splice_tsv tests/data/toy_splice_event_diff.tsv \
  --out_dir tests/tmp/splice_event_diff_generic \
  --organism human \
  --genome_build hg38 \
  --resources_dir tests/tmp/splice_bundle
```

Also run a no-bundle comparison:

```bash
geneset-extractors convert splice_event_diff \
  --splice_tsv tests/data/toy_splice_event_diff.tsv \
  --out_dir tests/tmp/splice_event_diff_nobundle \
  --organism human \
  --genome_build hg38 \
  --use_reference_bundle false
```

Validate outputs with:

```bash
geneset-extractors validate tests/tmp/splice_event_diff_generic
geneset-extractors validate tests/tmp/splice_event_diff_nobundle
```

Inspect:
- `geneset.tsv`
- `geneset.full.tsv`
- `genesets.gmt`
- `geneset.meta.json`
- `run_summary.json`
- `run_summary.txt`

Record:
- top 20 genes by absolute score
- top 20 genes by weight
- whether `KCNN4` is near the top in the bundled toy run
- whether the bundle changes rankings in intuitive ways
- whether the bundle resources are listed as used in metadata

### 7B. Additional toy stress tests if time allows

Repeat one or more toy runs to mirror the test suite behavior:
- change `impact_mode` from `conservative` to `none`
- change `event_dup_policy` from `highest_confidence` to `max_abs`
- change `ambiguous_gene_policy` from `drop` to `split_equal`
- create or reuse a low-confidence retained-intron-heavy toy file and confirm warnings appear

Expected outcomes:
- at least one score should change when impact or dedup settings change
- ambiguous toy genes such as `GENEA` and `GENEB` should only appear under non-drop handling
- retained-intron or low-confidence dominance should produce warnings rather than a silent clean run

### 7C. Real differential tables, if available

Only exercise `LeafCutter`, `MAJIQ`, or `Whippet` modes if you have a compatible real or public table. Do not invent family-specific fixtures if the repo does not include them.

For each real table you do have:
- run one bundle-backed configuration if the table can be matched to the local bundle
- run one no-bundle configuration
- record the actual `--tool_family` used
- if `LeafCutter` is used and a cluster-stats file exists, test both with and without `--cluster_stats_tsv`

What to expect from real differential tables:
- strong splicing-factor or EMT contrasts should give coherent RNA-processing, spliceosome, epithelial, junction, or EMT programs
- weak contrasts may legitimately return noisy or broad biology, but they still should not collapse entirely to giant genes with many measurable events
- the bundle should mildly denoise or harmonize; it should not completely rewrite the biology

## Step 8: run `splice_event_matrix` on toy PSI matrices and any staged real/public matrices

### 8A. Mandatory toy matrix run

Use the toy matrix and toy sample metadata. If `event_metadata_tsv` is optional for the toy generic case, follow the tested path. Include the toy coverage matrix if supported.

Example pattern:

```bash
geneset-extractors convert splice_event_matrix \
  --psi_matrix_tsv tests/data/toy_splice_matrix.tsv \
  --sample_metadata_tsv tests/data/toy_splice_sample_metadata.tsv \
  --coverage_matrix_tsv tests/data/toy_splice_coverage_matrix.tsv \
  --study_contrast condition_a_vs_b \
  --condition_a case \
  --condition_b control \
  --out_dir tests/tmp/splice_event_matrix_single \
  --organism human \
  --genome_build hg38 \
  --resources_dir tests/tmp/splice_bundle
```

Then test a grouped mode if supported, for example:
- `condition_within_group`
- `group_vs_rest`

Validate outputs and inspect:
- whether grouped runs produce `manifest.tsv`
- whether the toy grouped run emits two child outputs
- whether per-program outputs validate
- whether `welch_t` versus `mean_diff` materially changes results
- whether low-present events are filtered as expected
- whether top genes include one of `KCNN4`, `MAPK1`, or `GENEA`

### 8B. Matrix runs on staged public data

If Step 5 produced a real or public prepared TCGA SpliceSeq cohort, run `splice_event_matrix` on it using:
- `psi_matrix.tsv`
- `sample_metadata.tsv`
- `event_metadata.tsv`

Start with a simple tumor-versus-normal contrast.

What to expect from a successful public matrix run:
- one dominant biology axis, not a random gene list
- coherent enrichments in RNA splicing, RNA processing, epithelial or EMT state, junctions, cytoskeleton, cell cycle, or tissue identity depending on the cohort
- similar broad story between the matrix front-end and any compatible differential-table front-end

## Step 9: discover and run `pigean` and `eaggl`

Do not guess CLI syntax. First discover what is available.

Run:

```bash
command -v pigean || true
command -v eaggl || true
python - <<'PY'
import importlib.util
mods = ["pigean", "eaggl"]
for m in mods:
    print(m, bool(importlib.util.find_spec(m)))
PY
```

Then probe usage with whichever path exists:

```bash
pigean --help || true
eaggl --help || true
python -m pigean --help || true
python -m eaggl --help || true
```

Use the actual installed syntax, not an assumed syntax.

### What to analyze with `pigean`

Run `pigean` on:
- one bundle-backed `splice_event_diff` output
- one no-bundle `splice_event_diff` output
- one `splice_event_matrix` output
- one real or public dataset output if available

Record:
- top enriched gene sets
- posterior or score if reported
- whether the same broad biology appears across direct and matrix front-ends
- whether bundle-backed enrichments look mildly more specific or instead look over-biased

Expected high-level themes by dataset type:
- toy runs are smoke tests only and should not be overinterpreted biologically
- BRCA or PRAD-like carcinoma runs should often show epithelial, junction, EMT, or RNA-processing signal
- LUAD-like runs may add proliferation or cell-cycle signal and sometimes stronger RNA-processing signal
- KIRC or THCA-like runs may be dominated by tissue-state or metabolism themes, which is acceptable if reported honestly

### What to analyze with `eaggl`

Run `eaggl` on the same outputs, using the installed interface to infer latent factors or factor-like summaries.

Record:
- top factors
- factor weights or scores
- whether factors correspond to coherent splicing biology, broad RNA processing, epithelial or EMT state, proliferation, or obvious nuisance structure

If either tool is unavailable:
- state that clearly
- do not invent enrichment or factor results
- continue with the rest of the repo review

## Step 10: biology-versus-artifact review

For each major output, assess whether the signal looks like real biology or a likely artifact.

### Evidence for real biology

Look for patterns like:
- coherent spliceosome or RNA-binding-factor signal
- epithelial-versus-mesenchymal splicing programs in carcinoma datasets
- neuron-like, epithelial-like, or tissue-specific splicing programs in strong positive-control datasets
- specific exon-skipping or junction-switching programs concentrated in plausible pathway genes
- consistent signal across direct and matrix entrypoints
- stronger specificity with the bundle than without the bundle, without collapsing to a tiny prior-driven gene list

Dataset-specific interpretation rules:
- BRCA or PRAD tumor-vs-normal: lack of any epithelial, junction, EMT, or RNA-processing structure is suspicious
- LUAD tumor-vs-normal: some proliferation or broad tumor-state signal is normal, but pervasive retained introns and low-confidence events without coherent biology is suspicious
- KIRC or THCA tumor-vs-normal: broad tissue-state or metabolism signal may be real, but do not oversell it as a canonical splicing-factor program
- strong parser sanity controls such as tissue-separation splicing contrasts should yield crisp, interpretable biology; weak incoherent output there suggests a parsing or scoring issue

### Common artifacts to check for

Check explicitly for:
- genes ranking highly only because they have many measured events
- retained intron dominance suggestive of RNA quality or pre-mRNA contamination
- heavy dependence on low-support or low-confidence events
- top genes dominated by housekeeping, ribosomal, or mitochondrial signal
- splicing factors ranking simply because the fixture was built around them rather than because the event-to-gene map is informative
- bundle-induced over-connection, where priors overwhelm the raw data
- event alias collapse errors that merge unrelated events
- sign inconsistencies between `delta_psi`, `stat`, and emitted positive or negative GMTs
- grouped outputs that differ only because of bookkeeping bugs rather than biology

### Bundle-specific review

Compare bundle-backed and no-bundle runs.

Check whether the bundle:
- improves harmonization across tool-family fixtures or public datasets
- mildly downweights ubiquitous events
- applies conservative impact priors
- stays shrink-to-neutral when evidence is weak
- avoids dramatic re-ranking that looks prior-dominated

## Step 11: produce a final report

Write a report with these sections.

### 1. Environment and commit
- repo URL
- branch
- commit SHA
- python version
- whether `pigean` and `eaggl` were available

### 2. Dataset inventory and rationale
- which toy fixtures were run
- which real or public datasets were run
- why those datasets were chosen
- what biology you expected before running them

### 3. Implementation completeness
- were all four entrypoints present
- did docs exist
- did specs and CLI help look complete

### 4. Test summary
- full suite status
- targeted splicing test status
- notable failures

### 5. Functional smoke tests
- toy public staging workflow results
- bundle build results
- direct converter results
- matrix converter results
- validation status
- whether the observed toy outcomes match the expected toy contract

### 6. Top biological outputs
For each representative real or public run:
- top genes
- top gene sets from `pigean`
- top factors from `eaggl`
- one-paragraph interpretation
- whether observed biology matched what you expected for that dataset

### 7. Artifact assessment
For each representative run:
- likely real biology
- likely artifacts
- confidence level
- whether the warning system caught the issue

### 8. Recommendations
Provide concrete recommendations, ranked by impact:
- blocking bugs
- scoring or prior adjustments
- metadata or summary improvements
- documentation gaps
- public-dataset examples that should be added to docs
- future extensions that are reasonable after v1 is stable

## Deliverable style

Be specific and empirical.
Quote exact file names and commands.
Do not hide failures.
Do not overclaim biology from toy fixtures.
Prefer a clear distinction between:
- what the implementation does correctly
- what looks suspicious
- what remains untested
