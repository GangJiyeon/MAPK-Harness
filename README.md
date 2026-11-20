# MAPK-Harness

MAPK-Harness studies two related questions: how the measured activity of a compound
differs across KRAS variants, and how differences in signaling protein abundance
could change simulated ERK suppression between cells carrying the same mutation.
The project also evaluates the provenance, reproducibility, and failure handling
of the tools used to explore these questions.

This is an engineering and simulation project. Biochemical potency, cellular pERK,
and ODE output are different quantities; none alone establishes a patient's
response. The emphasis is on explicit assumptions, comparison against baselines,
and failure detection.

## Current status and 30-day plan

**Implemented:** BindingDB extraction and assay classification, mechanism-specific
datasets, scaffold/patent split preparation, five-seed Random Forest baseline,
and a G12D–WT exchange-assay selectivity table (120 shared molecules). The
corresponding scripts are in `notebooks/`. The generated data and result CSVs
under `data/processed/` are not committed; rerun the scripts with the source data
to reproduce them.

**Planned:** GCN evaluation on comparable splits; a verified Kholodenko-based
ODE; an explicitly assumed potency-to-ODE input mapping; paired virtual-cell
simulations and sensitivity analysis; provenance-aware prediction, simulation,
and drug-likeness tools; and a budget-matched agent harness with fault injection.
The current `backend/` modules are empty placeholders, and `frontend/` is an
unmodified Vite starter, not a working simulator.

The work is grouped into ten three-day cycles: study the relevant concepts,
implement one slice, then verify and record the result.

| Cycle | Deliverable |
|---|---|
| 01 | Define data flow, units, provenance, and the biochemical-to-cellular boundary |
| 02 | Freeze the evaluation candidate pool and oracle using the existing selectivity table; check train/test overlap and lookup rules |
| 03 | Train a GCN and compare it with RF on the same split and metrics |
| 04 | Reproduce the reference MAPK ODE and check trajectories, units, and conservation |
| 05 | Define and stress-test the assumed drug-effect mapping to ODE inputs |
| 06 | Simulate paired untreated/treated virtual cells with varied protein abundance |
| 07 | Calculate AUC-based responses and sensitivity under stated distribution assumptions |
| 08 | Connect prediction, lookup, simulation, and drug-likeness tools with explicit error and provenance fields |
| 09 | Compare equal-budget candidate search baselines and test fault detection |
| 10 | Reproduce the pipeline, document limitations, and run a small H₂ or LiH VQE-to-classical energy comparison |

The day-30 quantum calculation is a small, separate computational experiment. It
does not predict KRAS binding, IC50, or cell response. It is included as an
extension to the learning plan, not a dependency of the MAPK simulator or harness.

## The question

> Under the same KRAS mutation and the same drug exposure, how much heterogeneity
> in ERK suppression comes from cell-to-cell differences in signaling protein
> abundance, and which components contribute most?

**Genetic axis** — The same compound can have different measured activity against
different KRAS variants. G12D–WT selectivity is one assay-level comparison; it
does not by itself quantify clinical toxicity.

**Non-genetic axis** — Cells carrying the same mutation still respond differently.
This project varies protein abundance across simulated cells to examine one
possible extrinsic source of heterogeneity. It does not assume that this source
dominates every KRAS G12D cell line or treatment condition.

Two papers ground this framing:
- Filippi et al., *Cell Reports* (2016), PMID 27264188 — in MEK-ERK, the dominant
  source of cell-to-cell variation in ERK phosphorylation is extrinsic, particularly
  upstream of MEK.
- PMC5130170 — protein-abundance differences reproduce heterogeneous EGF responses,
  with specific proteins acting as sensitive nodes.

Both carry conditions. Neither claims to explain drug response in KRAS G12D cancer
in full, and the second uses a model with ERK–nuclear-pore regulation that the
simplified Kholodenko model here does not include.

## Planned architecture

```
Tools
  molecular model      SMILES + variant     → assay-level pIC50 estimate
  MAPK ODE             inhibition          → ERK-PP trajectory
  cell population      protein abundances  → response distribution
  sensitivity          simulation output   → contribution ranking

Validation
  applicability domain   max Tanimoto to training set
  numerical stability    negative concentrations, mass conservation
  self-consistency       CV→0 must reproduce the single-cell ODE
  fault injection        5 fault classes, measured detection rate

Harness
  ReAct loop over three proposed tools
  task: find compounds that inhibit G12D while sparing wild-type
  baselines: random / potency ranking / budget-matched search / ReAct
```

The planned ODE uses a reduced Kholodenko (2000) cascade. Its species, parameters,
feedback, and KRAS-to-RAF input mapping still need to be implemented and checked
against a reference before biological interpretation.

### Proposed tool contract

```
predict_pic50(smiles, variant) -> {
  pic50: float | null,
  source: "gcn" | "rf" | "lookup",
  assay_endpoint: string | null,
  model_or_data_version: string | null,
  in_domain: bool,
  error: string | null
}

simulate_population(pic50, dose_nM, n_cells, seed) -> {
  response_rate: float,
  auc_reduction: {mean, std, p10, p90},
  mapping_assumptions: object,
  converged: bool,
  error: string | null
}

check_druglike(smiles) -> {mw, logp, qed, passes, error}
```

The intended interface distinguishes measured lookup values, model predictions,
and assumption-based simulation outputs. Errors are returned explicitly so the
harness can test failure handling.

## Data

BindingDB release 202611 (2026-10-31), md5 verified.

### Filtering

| Step | Records |
|---|---|
| All measurements | 3,243,660 |
| KRAS in target name | 30,278 |
| G12D variant | 8,284 |
| Single target (complexes removed) | 4,956 |
| IC50 with `=` relation | 3,934 |
| Outliers removed (pIC50 > 10) | **3,881** |

Complexes were removed because records naming SOS1, RAF or PI3K alongside KRAS
measure protein-protein interaction inhibition, a different mechanism. Records with
`>` (2,620) and `<` (1,529) are kept separately — they carry information but are not
values.

### Assay mechanisms are not interchangeable

Mean pIC50 differs by up to three log units depending on what the assay measured.
Training on the pooled set would teach the model which experiment a record came
from rather than how tightly a compound binds.

| Mechanism | Molecules (G12D) | Mean pIC50 |
|---|---|---|
| exchange (SOS1-mediated nucleotide exchange) | 350 | 8.59 |
| binding | 233 | 8.26 |
| cellular (ERK phosphorylation in cells) | 921 | 7.30 |
| covalent labeling | 332 | — |
| no assay description | 437 | — |

**422 records were excluded after reading the assay text.** They measure covalent
labeling by mass spectrometry over up to 24 hours — a reaction rate, not a
reversible binding affinity — and the description specifies K-Ras (1-169) **G12C**,
contradicting the G12D label on those rows. A patent covering several variants had
its assay description generalized across all of them. This is the clearest example
in the dataset of why the validation layer exists: the label alone is not reliable.

### Datasets

| Name | Molecules | Scaffolds | Patents | Mean pIC50 |
|---|---|---|---|---|
| cellular | 1,450 | 361 | 8 | 7.16 |
| biochem (exchange + binding) | 652 | 233 | 13 | 8.47 |
| mixed | 2,102 | 594 | 21 | 7.57 |

The RF baseline has been run on all three datasets. GCN comparison on matched
splits is planned. Pooling assay mechanisms remains an experimental comparison,
not evidence that they are interchangeable.

### Selectivity dataset

Exchange assays run purified protein against several variants, so the same compound
often appears under more than one. Cellular assays use variant-specific cell lines
and barely overlap (G12D ∩ WT: 17 molecules).

| Pair, exchange assays | Shared molecules |
|---|---|
| G12D ∩ WT | **120** |
| G12D ∩ G12C | 46 |

```
selectivity = pIC50(G12D) − pIC50(WT)
range −1.55 to 2.42, median 1.43, sd 0.76, 13 negative
```

Potency and selectivity have a near-zero correlation in this table (**0.041**), and
the top five by each measure share **no molecules**. Ranking by potency alone does
not solve the task.

### Splits

Scaffold split with shuffled scaffold groups, verified to share no scaffold between
train and test. Patent split is reported as a secondary evaluation.

Two earlier split implementations were discarded. Assigning the largest scaffold
groups to the test set left only three scaffolds there; assigning the smallest left
a test set made almost entirely of singletons. Shuffling the groups before filling
the training set to 80% gives a test set whose scaffold sizes range from 1 to 71.

## Results

### Random Forest baseline

Morgan fingerprint, radius 2, 2048 bits. Five seeds.

| Dataset | Scaffold split R² | Random split R² | Mean predictor R² |
|---|---|---|---|
| biochem | 0.392 ± 0.136 | 0.600 ± 0.072 | −0.018 |
| mixed | 0.371 ± 0.247 | 0.612 ± 0.015 | −0.066 |
| cellular | 0.296 ± 0.165 | 0.493 ± 0.048 | −0.144 |

Random split sits roughly 0.20 above scaffold split across all three, which
quantifies how much of the apparent performance is scaffold memorization. RMSE
of 0.88–0.96 means predictions are off by about an order of magnitude in IC50 —
useful for triage, not for precision.

The variance on `mixed` is large because a single 140-molecule scaffold may or may
not land in the test set. GCN comparisons will use ten seeds.

### GCN, cell population, sensitivity, harness

In progress.

## Limitations

- **Single source.** All G12D records come from US patents; none are peer reviewed.
  Applicants have an incentive to report their own compounds favorably.
- **Applicant concentration.** Ten patents account for 90.7% of molecules, so the
  model may only generalize within the chemical space those groups explored. This is
  the direct motivation for applicability-domain checks in the harness.
- **No experimental validation.** The cell population model is mechanistic and
  unvalidated. It is a hypothesis explored in simulation.
- **Evaluation target.** The candidate-pool oracle and train/test separation are
  still to be fixed. A tool-derived exhaustive ranking measures tool use under a
  call budget, not drug efficacy; a held-out measured oracle answers a different
  question and must be labeled separately.
- **No literature anchor.** MRTX1133 and other reference inhibitors come from ChEMBL
  and are absent here. Acquiring them separately is pending.
- **Mechanism coverage.** Selectivity rests on 120 exchange-assay molecules;
  the binding set (34) was too small to include.

## Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python notebooks/10_env_check.py
```

```
notebooks/02_extract_kras.py      BindingDB → KRAS records
notebooks/03_classify.py          variant, complex status, IC50 relation
notebooks/05_patent_check.py      scaffold and patent distribution
notebooks/06_assay_join.py        join assay descriptions
notebooks/07_assay_classify.py    classify by mechanism
notebooks/08_build_datasets.py    three datasets
notebooks/09_featurize_split.py   fingerprints and scaffold splits
notebooks/11_rf_baseline.py       RF baseline
notebooks/12_sanity_checks.py     outliers, fingerprints, prediction spread
notebooks/13_mutation_survey.py   all variants through the pipeline
notebooks/14_mutation_pipeline.py selectivity dataset
```

Decision records are in `docs/adr/`.

## Citation

Gilson MK, Liu T, Baitaluk M, Nicola G, Hwang L, Chong J. BindingDB in 2015: A
public database for medicinal chemistry, computational chemistry and systems
pharmacology. *Nucleic Acids Research* 44:D1045–D1053.
doi:[10.1093/nar/gkv1072](https://doi.org/10.1093/nar/gkv1072)
