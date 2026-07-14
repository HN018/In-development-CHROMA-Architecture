# CHROMA v6 — System Design

**C**holesterol **H**omeostasis **R**egulator discovery via **O**rchestrated **M**ulti-method **A**nalysis

**v6 · July 2026**
**Author:** Haoning Yang · Saheki Lab, LKCMedicine, NTU Singapore
**Status:** Design document — supersedes CHROMA v5
**Supersedes:** CHROMA v5 (comparability-gated convergence)

---

## 0. Why v6 exists, and what it removes

v5 was a convergence engine. Its central premise — Principle 2, *"convergence is the unit
of evidence"* — was that many orthogonal screens each provide weak evidence, and
consistency across them is the signal. Layers 2a (RRA) and 4.0 (the Comparability Audit)
were the heart of the system.

**A screen-by-screen audit of the lab's actual data showed there is no convergence to
find.**

The audit is documented in §14. In summary:

| Dataset | v5 assumed | v6 finds |
|---|---|---|
| Minyoung ddFP CRISPR | one of N screens | **the screen.** Bidirectional, 3 replicates. |
| Dylan AmphoB CRISPR | co-equal screen | n=1, no replication, unidirectional |
| Golgi-IP proteomics | orthogonal screen | different compartment, different quantity |
| 4× bulk RNA-seq | orthogonal screens | **not screens.** Response readouts. |
| Public screens (Pfeffer, ALOD4, Trinh, van den Boomen) | RRA partners | none phenotype- and format-compatible |

The system was built to aggregate evidence the project does not have.

**The critical finding is not that v5 would break — it is that v5 would succeed at the
wrong thing.** Registering the four RNA-seq datasets as screens would produce strong
convergence on SREBP2 targets (HMGCR, LDLR, SQLE, INSIG1), because those genes are the
cell's *response* to perturbed cholesterol. LOPO would look excellent. Canonical genes
would recover into P1. Every acceptance criterion in v5 §11 would pass. The pipeline
would report a triumphant sensitivity check while having discovered nothing.

v5's `validation_class` does not catch this. It flags genes that are *already known*. It
does not flag genes that are *downstream consequences*. That is a distinct failure mode,
and closing it is the core of v6.

### 0.1 What v6 removes

| Removed | Reason |
|---|---|
| **Layer 2a — RRA** | Requires ≥3 comparable ranked lists. There is one screen. |
| **§4.0 Comparability Audit** | Its only function was gating entry to RRA. |
| **`pyrra` dependency + R parity test** | No longer used. |
| **MOFA+ (§4.6)** | Multi-omics factors across datasets measuring different things in different compartments. |
| **Public benchmark path (§11 criterion 10)** | No public screen is comparable. Recovering GGT7 on an incomparable screen validates nothing about this pipeline. |
| **Fixed-weight blending** | Already removed in v3; stays removed. |

### 0.2 What v6 adds

| Added | Purpose |
|---|---|
| **Responder/Regulator Discrimination (Layer 2)** | **The core contribution.** Separates genes that *regulate* PM cholesterol from genes the cell *turns on when cholesterol drops*. |
| **`downstream_responder_risk` concern flag** | Travels with the gene; never deletes it. |
| **Bidirectional analysis** | The screen's `high` arm (negative regulators) was unused. |
| **Orthogonal Filter Bank (Layer 3)** | Replaces RRA. Each filter removes one named false-positive class. |
| **Dataset registration classes** | `screen` / `compartment_feature` / `response_signature`. Prevents a non-screen being ranked as one. |
| **Rounds-asymmetry handling** | The screen's own design biases against essential genes. Must be corrected upstream. |

### 0.3 What v6 keeps, unchanged

Every design *principle* from v2–v5 survives. Only the architecture changes.

- Nothing is silently dropped; exclusion means a labelled, recoverable bin.
- Screens are understood before they are scored (`.desc.md` required).
- The LLM advises; it never decides.
- Adversarial debate is useful when grounded.
- Heavy methods are optional, not assumed.
- Transparent P1–P4 bins; no hard vetoes.
- Benchmark recovery is not discovery.
- Provenance-aware annotation (`documentation_provenance`).
- Discovery must outrun the known-positive set (LOPO, label-shuffle).

---

## 1. Design Principles

1. **Nothing is silently dropped.** Every gene keeps a continuous score and a record to
   the end. "Exclusion" means deprioritized into a labelled bin with a logged, recoverable
   reason.

2. **One screen, analyzed properly, beats seven screens aggregated badly.** *(NEW,
   replaces v5 Principle 2)* Convergence is only evidence when the things converging are
   independent measurements of the same quantity. Where that condition fails, apparent
   convergence is an artifact of shared structure, not a signal.

3. **A regulator is not a responder.** *(NEW — the central principle of v6)* A gene that
   changes when cholesterol changes is not thereby a gene that *controls* cholesterol.
   Transcriptional readouts of a perturbed cell measure the response, and mistaking that
   response for a set of causes is the dominant false-positive mode in this field. The two
   are separated explicitly, by construction, and never conflated.

4. **Every dataset declares what kind of thing it is.** *(NEW)* A dataset is a `screen`, a
   `compartment_feature`, or a `response_signature`. Only screens rank genes. This is
   asserted by a human in `.desc.md`, never inferred from file format.

5. **Filters are named by what they remove.** *(NEW)* Each element of the Filter Bank
   exists to eliminate one specific, named false-positive class. A filter that cannot say
   what it removes does not belong in the pipeline.

6. **Screens are understood before they are scored.** Each dataset carries a `.desc.md`
   context object, read before normalization and passed to every downstream agent.

7. **The LLM advises; it never decides.** Grounded, cited, archived, hypothesis-generating
   annotation on top of a fully reproducible computational core.

8. **Discovery must outrun the known-positive set.** Validated on the ability to recover
   withheld biology (LOPO), not on re-ranking canonical genes.

9. **Heavy methods are optional, not assumed.** The classical core runs on a laptop.

10. **Benchmark recovery is not discovery.** Unchanged from v5, but note: with the public
    benchmark path removed, there is no benchmark. The only bar is a real
    `discovery_candidate` from this lab's own data.

---

## 2. Architecture Overview

Five layers. Layers 1–3 and 5 deterministic and reproducible; Layer 4 advisory.

```
┌─────────────────────────────────────────────────────────────────────┐
│ LAYER 0 — REGISTRATION (deterministic)                              │
│   Every dataset declares registration_class:                        │
│     screen | compartment_feature | response_signature               │
│   Only `screen` may produce a ranked gene list.                     │
│   Reader → DatasetContext (+ documentation_provenance)              │
└──────────────────────────────┬──────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│ LAYER 1 — PRIMARY SCREEN ANALYSIS (deterministic)                   │
│   The ddFP screen, analyzed properly. This is the whole ranking.    │
│   1a  MAGeCK MLE, batch-aware (S1 | S2,S3), replicate-native        │
│   1b  Bidirectional: `low` (positive reg.) AND `high` (negative reg.)│
│   1c  Primary contrast = high vs low (rounds-matched, drift cancels)│
│   1d  Secondary = vs main (rounds-asymmetric — essentiality-flagged) │
│   → ScoreVector per gene, per direction                             │
└──────────────────────────────┬──────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│ LAYER 2 — RESPONDER / REGULATOR DISCRIMINATION (NEW; the heart)     │
│   2a  Build cholesterol-response signature from CLEAN genetic       │
│       contrasts only (ORP9 KO, GRAMD1 TKO, QKO vs WT)               │
│   2b  Score every candidate for signature membership                │
│   2c  → downstream_responder_risk flag (annotates; never deletes)   │
│   A screen hit that is ALSO a signature member is a suspected       │
│   feedback node, not a regulator.                                   │
└──────────────────────────────┬──────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│ LAYER 3 — ORTHOGONAL FILTER BANK (deterministic)                    │
│   Replaces RRA. Each filter removes one named false-positive class. │
│   3a  AmphoB dose-concordance    → removes ddFP reporter artifacts  │
│   3b  DepMap essentiality        → removes rounds-asymmetry artifacts│
│   3c  Golgi-IP compartment prior → positive prior (not a filter)    │
│   3d  RWR network propagation    → positive prior; rescues sub-thresh│
│   3e  Co-essentiality proximity  → positive prior                   │
│   → ConvergenceVector per gene                                      │
└──────────────────────────────┬──────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│ LAYER 4 — PRIORITIZATION + ANNOTATION                               │
│   4a  PriorityScore (PU learner) + NoveltyScore  [deterministic]    │
│       LOPO across ≥3 withheld modules → rank distribution           │
│   4b  LLM annotation [advisory]: Screen Reader, Senior Postdoc,     │
│       Proposer, Skeptic, Orchestrator, Narrative Generator          │
│       Debate-value ablation decides single-agent vs. full debate    │
└──────────────────────────────┬──────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────┐
│ LAYER 5 — TRIAGE & OUTPUT (deterministic)                           │
│   P1–P5 bins × direction (positive | negative regulator)            │
│   Concern bins incl. downstream_responder_risk                      │
│   Nothing dropped.                                                  │
└─────────────────────────────────────────────────────────────────────┘
```

### Key structural changes from v5

- **Layer 2 is new and is the system's actual contribution.** v5's Layer 2 was
  convergence; v6's Layer 2 is discrimination. This is the layer that makes the project
  publishable, because it addresses a problem every cholesterol screen has and nobody
  handles well.
- **The primary screen is now a layer, not an input.** v5 treated it as one of N. It is
  N.
- **Every gene now carries a direction** (positive or negative regulator). v5's
  single-axis design could not represent this.

---

## 3. Layer 0 — Registration

### 3.1 Registration classes

Every dataset in `chroma_config.yaml` must declare a `registration_class`. This is a
human assertion in `.desc.md`. It is never inferred from file format, and the system will
not run without it.

| Class | May rank genes? | Meaning |
|---|---|---|
| `screen` | **Yes** | Perturbs many genes, reads one phenotype. Assigns a phenotype to a perturbed gene. |
| `compartment_feature` | No | Measures protein/gene presence at a compartment. Feeds §3c prior. |
| `response_signature` | No | Reads many transcripts after perturbing one gene/state. Feeds Layer 2. |

**Only `screen` datasets enter Layer 1.** A `response_signature` registered as a `screen`
is the failure mode described in §0; the class system exists to make it impossible rather
than merely discouraged.

### 3.2 Required `.desc.md` fields

Carried from v5, plus three new required fields:

```python
{
  "dataset_id": str,
  "registration_class": "screen" | "compartment_feature" | "response_signature",  # NEW
  "cell_type": str,
  "perturbation": str,
  "perturbation_class": "null_ko" | "hypomorph" | "knockdown" |                  # NEW
                        "overexpression" | "genotype_panel" | "pharmacological",
  "phenotype_axis": str,   # NEW — what compartment/pool, and sign direction
  "target_compartment": str,
  "n_replicates": int,
  "replicate_independence": str,   # NEW — are the replicates actually independent?
  "phenotype": str,
  "score_interpretation": {"positive": str, "negative": str, "zero": str},
  "positive_controls": List[str],
  "known_limitations": str,
  "confounder_flags": List[str],
  "documentation_provenance": "contemporaneous" | "reconstructed",
  "reconstruction_note": Optional[str],
  "raw_description": str
}
```

**`phenotype_axis` and `perturbation_class` are new and are load-bearing.** v5's
Comparability Audit tested *statistical* comparability (gene-universe overlap, rank
distributions) but not *semantic* or *perturbational* comparability. It would have passed
the Pfeffer PFO* screen — which measures lysosomal cholesterol with a sign that flips
relative to PM — and the ALOD4 gene-trap screen, whose hypomorphs are not comparable to
null KOs. Both would have corrupted a rank aggregation while passing every automated
check. These two human-asserted fields close that gap at negligible cost.

**`replicate_independence` is new** because "n_replicates: 3" is a lie when two of the
three are parallel dishes from one batch.

---

## 4. Layer 1 — Primary Screen Analysis

**This layer is the ddFP screen, and the ddFP screen is the project.** It deserves the
attention v5 gave to the entire convergence engine.

### 4.1 The screen

`crispr_hela_ddfp_ldl_restoration` (Minyoung). HeLa + ddFP-PM-RA / B-GRAM-W reporter,
GeCKO v2 pooled KO, LPDS + mevastatin depletion → 50 µg/ml LDL 48 h → FACS.

Three replicates (S1; S2, S3). Three pools per replicate: `main` (1 sorting round),
`low` and `high` (0.7% gates, 3 sorting rounds). Nine libraries, all sequenced
separately.

### 4.2 Bidirectional analysis — NEW

The original analysis used only the `low` arm. **The `high` arm is equally valid and
reports the opposite phenotype.** Roughly half the screen's information was unused.

| Arm | Phenotype | Gene class |
|---|---|---|
| `low`-enriched | KO **fails to restore** PM cholesterol on LDL | **positive regulator** |
| `high`-enriched | KO **overshoots** PM cholesterol on LDL | **negative regulator** |

**Positive and negative regulators are never aggregated into one list, one score, or one
bin.** They are separate rankings, separately validated, separately triaged. A gene strong
in `low` and a gene strong in `high` are both real and both interesting; they are not
convergent with each other.

### 4.3 Contrasts and the rounds-asymmetry problem — CRITICAL

`main` was collected after **one** sorting round; `high` and `low` after **three**. No
round-3 `main` exists. This is unfixable in the data.

| Contrast | Rounds matched | Status |
|---|---|---|
| `high vs low` | **3 vs 3 — yes** | **PRIMARY.** Drift cancels. |
| `low vs main` | 3 vs 1 — no | Secondary. Essentiality-flagged. |
| `high vs main` | 3 vs 1 — no | Secondary. Essentiality-flagged. |

**Why this matters more than it looks.** The `vs main` contrasts confound the cholesterol
phenotype with two extra rounds of culture, sorting stress, and proliferation-driven sgRNA
dropout experienced only by the `high`/`low` arms. Essential genes drop out of the 3-round
arms *for reasons unrelated to cholesterol*, and therefore appear artifactually depleted.

**The screen's own design is quietly vetoing essential genes** — precisely the gene class
that CHROMA's no-veto principle and P3 bin exist to protect. The bias is introduced
upstream of anything the pipeline can see. It must be corrected here, in Layer 1, or it
cannot be corrected at all.

**Handling:**
- `high vs low` is the primary ranking. Use it wherever possible.
- `vs main` contrasts are retained (they catch genes moving in only one direction) but
  every gene sourced from them carries a mandatory `rounds_asymmetry_risk` flag and is
  cross-checked against DepMap essentiality (§3b) before entering P1/P2.

### 4.4 Batch structure

S1 was run at a different time, with a different cell batch. S2 and S3 were parallel
dishes from one batch. **Effective n ≈ 2, not 3.**

- Run **MAGeCK MLE with batch as a covariate** (S1 | S2,S3).
- Do **not** average S1/S2/S3 into one list before analysis. A gene supported only by
  S2+S3 is a batch-effect candidate, and averaging hides exactly that.
- Report, per gene, whether support is S1-inclusive or S2/S3-only. Surface this in the
  dossier.

### 4.5 Output

`ScoreVector` per gene, per direction (positive / negative), with:
`beta`, `se`, `fdr`, `contrast_source` (high_vs_low | low_vs_main | high_vs_main),
`batch_support` (S1_inclusive | S2_S3_only), `rounds_asymmetry_risk` (bool).

No cutoff is applied. No hit/no-hit decision is made here.

---

## 5. Layer 2 — Responder / Regulator Discrimination

**This is the heart of v6 and the project's actual scientific contribution.**

### 5.1 The problem

A CRISPR screen for PM cholesterol will recover two kinds of gene:

- **Regulators** — genes that control PM cholesterol. What the project wants.
- **Responders / feedback nodes** — genes the cell activates *because* cholesterol
  changed. SREBP2 targets are the canonical case. Knock one out and PM cholesterol may
  well change, because you have broken the feedback loop — but the gene is a consequence
  of cholesterol sensing, not an upstream controller of PM cholesterol levels.

Every cholesterol screen has this problem. Most do not address it. **The lab has the
resource to address it, and did not realize it.**

### 5.2 The resource — the RNA-seq panel, reframed

Four bulk RNA-seq datasets exist. v5 would have registered them as screens, which is the
failure mode of §0. They are not screens.

**They are the answer key, not the exam.** They describe what a cholesterol-perturbed cell
*looks like* transcriptionally. That is precisely the tool needed to identify responders.

### 5.3 Building the signature (2a)

Use **only the genetically clean cholesterol perturbations**:

| Dataset | Contrast | Use |
|---|---|---|
| `rnaseq_orp9_gramd1_panel` | ORP9 KO, GRAMD1 TKO, QKO vs WT | **PRIMARY** |
| `rnaseq_ipsc_astrocyte` | GRAMD1 TKO vs WT, *within cell type* | Supporting |
| `rnaseq_myotube_differentiation` | — | **EXCLUDED** |
| `rnaseq_hela_trim28_ko` | — | **EXCLUDED — circular** |

**Exclusions are not optional:**

- **Myotube differentiation** — differentiation changes thousands of genes. The cholesterol
  signature is buried inside a myogenic program and cannot be separated from it with this
  design. There is no cholesterol-specific perturbation, only a *correlation* between
  differentiation state and PM cholesterol.

- **TRIM28 KO — circularity.** TRIM28 was identified *by the ddFP screen*. Feeding
  TRIM28-KO transcriptomics into the evidence layer means genes downstream of TRIM28
  receive convergence credit tracing back to the screen that nominated TRIM28. It would
  present as independent cross-modality confirmation. **It is not.** It is the screen's
  own output, laundered through a second assay. v5 §5.2 guards against circularity in the
  *label* set; this is circularity in the *evidence* set and is otherwise uncaught.

  This dataset is **barred from the evidence layer**, not down-weighted. It is available
  for TRIM28 mechanism at Phase 5 (human review) only.

Expect the resulting signature to be **SREBP2-dominated**. That is the point, not a bug.

### 5.4 Scoring and flagging (2b, 2c)

For every gene in the screen's ranked list, compute `responder_score`: membership strength
in the cholesterol-response signature (e.g. signed rank of its differential expression
across the clean contrasts).

```
downstream_responder_risk = HIGH   if gene is a strong signature member
                          = MEDIUM if weak/inconsistent signature member
                          = LOW    if absent from signature
```

**Interpretation, stated plainly:**

| Screen hit? | In signature? | Reading |
|---|---|---|
| Yes | **No** | **Strongest candidate.** A regulator the cell does not transcriptionally respond with. |
| Yes | Yes | **Suspected feedback node.** Flag `downstream_responder_risk: HIGH`. Not deleted. |
| No | Yes | A responder. Correctly not a candidate. |

**This is a concern flag, not a veto.** Per Principle 1, the gene keeps its score and its
record. A real regulator can also be transcriptionally responsive — SREBF2 itself is both.
The flag says *"this gene's screen phenotype may be a feedback artifact; look harder,"*
not *"this gene is wrong."* The Skeptic (Layer 4b) is required to address it; the human
reviewer (Phase 5) adjudicates it.

### 5.5 Why this is the contribution

Layer 2 is what makes v6 a paper rather than a pipeline run. It answers a question the
field has: *how do you tell a cholesterol regulator from a cholesterol responder?* The
answer — use transcriptomic response data as a negative filter rather than as evidence —
is generalizable beyond this lab and beyond cholesterol.

---

## 6. Layer 3 — Orthogonal Filter Bank

Replaces RRA. **Each filter is named for the false-positive class it removes.** A filter
that cannot name what it removes does not belong here (Principle 5).

### 6.1 Filters (remove false positives)

| # | Filter | Removes | Source |
|---|---|---|---|
| 3a | **AmphoB dose-concordance** | **ddFP reporter artifacts** — genes affecting probe folding/maturation rather than cholesterol | Dylan's screen |
| 3b | **DepMap essentiality cross-check** | **Rounds-asymmetry artifacts** (§4.3) — essential genes artifactually depleted in 3-round arms | Public (DepMap) |

**3a — AmphoB dose-concordance.** Dylan's screen has n=1 per arm and cannot rank. It can,
however, *corroborate*: a completely different assay (survival under a cholesterol-binding
cytolysin) with completely different failure modes, reading the same phenotype in the same
direction. That is exactly what an underpowered second screen is good for.

Require dose-concordance before use: a gene must be enriched at **both** 400 and 500,
within a library, to count as corroborated. This converts the dose axis from a weak
pseudo-replicate into a noise filter against amphotericin's stochastic bottleneck.

libA and libB are genuinely independent (different sgRNAs, different transduction) —
corroboration in both is stronger than in one.

*This screen is used as a filter, never as a ranking input.* Rank aggregation gives every
input list equal authority; giving an unreplicated single-dish selection the same
authority as a three-replicate FACS screen would be indefensible. RRA has no weighting
mechanism. A filter does not need one.

### 6.2 Priors (promote true positives — never veto)

| # | Prior | Adds | Source |
|---|---|---|---|
| 3c | **Golgi-IP compartment concordance** | Membrane-trafficking plausibility | Golgi-IP proteomics |
| 3d | **RWR network propagation** (α ≈ 0.85) | Rescues sub-threshold genes with strong network support | STRING, BioGRID, DoRothEA/CollecTRI, Reactome FI, DepMap co-ess. |
| 3e | **Co-essentiality proximity to sterol module** | Functional proximity to SCAP/SREBF2/HMGCR/NPC1/LDLR etc. | DepMap |

**3c note.** The Golgi-IP measures a *different compartment* (Golgi, not PM) and a
*different quantity* (protein abundance, not phenotype). It is a **prior**, not evidence.
A protein whose Golgi abundance is cholesterol-responsive is interesting; it is not thereby
a PM cholesterol regulator. Also: OSW-1 causes Golgi fragmentation, which changes what the
HA pulldown physically captures for reasons independent of cholesterol. Declare this; do
not normalize it away.

**3d note.** RWR is the layer that keeps a single screen from being a single point of
failure. Seeded from Layer 1's ScoreVectors, it surfaces genes that are individually
sub-threshold but sit in a strongly-supported neighborhood. This is where CHROMA's
"convergence" premise still lives — as network convergence, not screen convergence.

---

## 7. Layer 4 — Prioritization and Annotation

### 7.1 Prioritization (deterministic) — largely unchanged from v5

- **PriorityScore** — calibrated, SHAP-explained bagging-PU classifier (`pulearn`).
- **NoveltyScore** — label-independent; convergence minus literature footprint.
- Never collapsed into one number.
- **Run separately per direction** (positive / negative regulator). NEW.

**Validation** (unchanged from v5 §5.2):
- **LOPO across ≥3 independently withheld modules** (SREBP-processing, LDLR-trafficking,
  lysosomal cholesterol export). Report the **distribution** of recovery ranks (median,
  IQR), never a single run.
- **Held-out positive recovery** — hide 50%, report rank distribution.
- **Label-shuffle null** — must collapse to chance, or the run is rejected as leakage.

**Metrics:** Precision@50/@100 and LOPO recovery distribution are primary. AUROC is
report-only, never gated (with ~35 positives in ~20k genes it is dominated by easy
negatives).

### 7.2 Annotation (LLM; advisory) — unchanged from v5, with one addition

Agents 0–5 as in v5: Screen Reader, Senior Postdoc (knowledge bank), Proposer, Skeptic,
Orchestrator, Narrative Generator, plus the Fabrication Checker. All grounding rules
carry over: every factual claim needs a verifiable database identifier or is explicitly
labelled a hypothesis. Fabrication ceiling 2%. The LLM cannot delete a gene.

**Provenance discount (v5 §6.4) retained:** a Skeptic concern resting on a
`documentation_provenance: reconstructed` field is downgraded one severity level unless
independently corroborated by a database record. This matters here — Dylan's screen and
three of the four RNA-seq datasets are reconstructed.

**NEW — the Skeptic must address `downstream_responder_risk`.** For any candidate flagged
HIGH, the Skeptic is required to argue the responder case explicitly, and the Proposer to
rebut it. This is the single most consequential concern in v6 and cannot be left implicit.

**Debate-value ablation (v5 §6.9) retained.** Single-agent grounded annotator vs. full
debate on a 30–50 gene sample; compare fabrication rate, concern precision, cost/latency.
The 6-agent structure ships as default only if it earns it.

---

## 8. Layer 5 — Triage and Output

### 8.1 Bins — now two-dimensional

Every gene is binned **within its direction**. There are two independent triage tables:
positive regulators and negative regulators.

- **P1 — Convergent-canonical.** High PriorityScore, resembles known biology.
  Confidence-building; not the discovery payload.
- **P2 — Convergent-novel.** High screen signal AND high novelty AND
  `downstream_responder_risk: LOW`. **The primary discovery bin.**
- **P3 — Orphan / essential / compartment-resident.** Uncharacterized, possibly essential,
  localizing to a relevant membrane. Explicitly promoted, never vetoed. *Expect this bin to
  matter more in v6 than v5 — the rounds asymmetry (§4.3) actively suppresses essential
  genes, so anything surfacing here despite that headwind deserves attention.*
- **P4 — Network-rescued.** Sub-threshold in the screen but strongly supported by RWR /
  co-essentiality. Awaiting replication.
- **P5 — Suspected feedback node.** NEW. Strong screen signal but
  `downstream_responder_risk: HIGH`. **Not a rejection bin.** These genes are real screen
  hits whose phenotype may be a feedback artifact. Some will be genuine regulators that are
  also transcriptionally responsive. Requires human adjudication.

### 8.2 Concern flags (orthogonal, sortable)

`downstream_responder_risk` **(NEW)** · `rounds_asymmetry_risk` **(NEW)** ·
`reporter_artifact_risk` **(NEW)** · `batch_support: S2_S3_only` **(NEW)** ·
`essential_gene_risk` · `fitness_artifact_risk` · `cell_line_artifact_risk` ·
`prior_clinical_failure`

A concern never deletes a gene. It annotates it, travels with it, and is sortable.

### 8.3 `validation_class` — retained, simplified

- `benchmark_recovery` — gene already has a peer-reviewed characterization.
  **Sensitivity check only. Never cited as evidence of discovery capability.**
- `discovery_candidate` — Senior Postdoc finds no prior characterization.
  **This is the class that matters.**

With the public benchmark path removed, there is no external benchmark. Canonical-gene
recovery (LDLR, NPC1, SCAP, MBTPS1) from *this lab's own screen* is the sanity check.

---

## 9. Data Layer

```yaml
scope:
  name: "pm_cholesterol"
  canonical_module: [SCAP, SREBF2, SREBF1, MBTPS1, MBTPS2, INSIG1, INSIG2,
                     HMGCR, LDLR, LDLRAP1, NPC1, NPC2, MYLIP]
  profile: "lite"

datasets:

  # ---- SCREENS (may rank genes) ----
  - id: "crispr_hela_ddfp_ldl_restoration"
    registration_class: "screen"           # NEW
    role: "primary"                        # NEW
    source: "private"
    format: "mageck_mle"
    perturbation_class: "null_ko"
    phenotype_axis: "pm_accessible_cholesterol_ldl_stimulated"
    bidirectional: true                    # NEW
    replicates: {S1: {batch: A}, S2: {batch: B}, S3: {batch: B}}
    pools: [main, low, high]
    contrasts:
      primary:   ["high_vs_low"]           # rounds-matched
      secondary: ["low_vs_main", "high_vs_main"]   # rounds-asymmetric
    rounds_asymmetry: true                 # NEW — main=1 round, high/low=3 rounds
    description: "data/private/crispr_hela_ddfp_ldl_restoration.desc.md"

  - id: "crispr_hela_amphotericin_survival"
    registration_class: "screen"
    role: "filter_only"                    # NEW — corroborates, never ranks
    source: "private"
    format: "mageck_test"
    perturbation_class: "null_ko"
    phenotype_axis: "pm_accessible_cholesterol_steady_state"
    bidirectional: false                   # cannot detect negative regulators
    libraries: [libA, libB]                # independent
    arms: [DMSO, ampho_400, ampho_500]     # n=1 each
    require_dose_concordance: true         # NEW — must enrich at BOTH doses
    description: "data/private/crispr_hela_amphotericin_survival.desc.md"

  # ---- COMPARTMENT FEATURE (may not rank genes) ----
  - id: "proteomics_golgi_ip"
    registration_class: "compartment_feature"
    target_compartment: "golgi"
    conditions: [WT, QKO, mevastatin, OSW1]
    caveat: "OSW-1 causes Golgi fragmentation — confounded, declare not correct"
    description: "data/private/proteomics_golgi_ip.desc.md"

  # ---- RESPONSE SIGNATURE (may not rank genes) ----
  - id: "rnaseq_orp9_gramd1_panel"
    registration_class: "response_signature"
    signature_role: "primary"
    contrasts: ["ORP9KO_vs_WT", "GRAMD1TKO_vs_WT", "QKO_vs_WT"]

  - id: "rnaseq_ipsc_astrocyte"
    registration_class: "response_signature"
    signature_role: "supporting"
    contrasts: ["TKO_vs_WT_within_celltype"]   # cross-celltype contrast EXCLUDED

  - id: "rnaseq_myotube_differentiation"
    registration_class: "response_signature"
    signature_role: "excluded"
    exclusion_reason: "Differentiation program confounds; no cholesterol-specific perturbation"

  - id: "rnaseq_hela_trim28_ko"
    registration_class: "response_signature"
    signature_role: "excluded"
    exclusion_reason: "CIRCULAR — TRIM28 was nominated by the primary screen"
    barred_from_evidence_layer: true          # NEW — hard bar, not a down-weight

responder_discrimination:                     # NEW — Layer 2
  enabled: true
  signature_source: ["rnaseq_orp9_gramd1_panel", "rnaseq_ipsc_astrocyte"]
  flag: "downstream_responder_risk"
  action: "flag_only"                         # NEVER delete (Principle 1)
  expect_srebp2_dominated: true

filter_bank:                                  # NEW — Layer 3, replaces RRA
  filters:
    - {id: amphob_dose_concordance, removes: "ddfp_reporter_artifact"}
    - {id: depmap_essentiality,     removes: "rounds_asymmetry_artifact"}
  priors:
    - {id: golgi_compartment,  source: proteomics_golgi_ip, as_veto: false}
    - {id: rwr_propagation,    alpha: 0.85,
       layers: {string: 0.30, biogrid: 0.15, dorothea: 0.20,
                reactome: 0.10, depmap_coess: 0.25}}
    - {id: coessentiality,     anchor: scope.canonical_module}

# REMOVED in v6: rank_aggregation (RRA), comparability_audit, mofa

prioritization:
  priority_model: pu_bagging
  novelty_score: true
  per_direction: true                         # NEW
  primary_metrics: [precision_at_50, lopo_recovery_distribution, holdout_median_rank]
  auroc: report_only

validation:
  lopo: true
  lopo_modules: ["srebp_processing", "ldlr_trafficking", "lysosomal_cholesterol_export"]
  holdout_fraction: 0.5
  label_shuffle_null: true

annotation:
  model: "TODO — re-verify availability/pricing before each run"
  temperature: 0
  top_n_candidates: 200
  senior_postdoc_rescue_n: 20
  require_verifiable_ids: true
  provenance_discount_rule: true
  skeptic_must_address_responder_risk: true   # NEW
  debate_mode: "ablation_pending"
  fabrication_rate_ceiling: 0.02
  can_exclude_genes: false

triage:
  bins: [P1_convergent_canonical, P2_convergent_novel,
         P3_orphan_essential_compartment, P4_network_rescued,
         P5_suspected_feedback_node]         # NEW
  per_direction: true                        # NEW
  validation_class: [benchmark_recovery, discovery_candidate]
  hard_vetoes: none
  silent_filtering: forbidden
```

---

## 10. Software Architecture

```
chroma/
├── chroma_config.yaml
├── pyproject.toml
│
├── chroma/
│   ├── cli.py
│   │
│   ├── workers/
│   │   ├── dataset_registrar.py            # NEW — enforces registration_class
│   │   ├── screen_normalizer.py
│   │   ├── mageck_runner.py                # NEW — batch-aware MLE, bidirectional
│   │   ├── responder_signature.py          # NEW — Layer 2a/2b: THE CORE
│   │   ├── filter_bank.py                  # NEW — Layer 3, replaces rra_aggregator
│   │   ├── network_propagator.py
│   │   ├── depmap_coessentiality.py
│   │   ├── compartment_concordance.py
│   │   ├── pu_learner.py
│   │   ├── novelty_scorer.py
│   │   └── validator.py                    # multi-module LOPO distribution
│   │
│   ├── agents/                             # unchanged from v5
│   │   ├── screen_reader.py
│   │   ├── senior_postdoc.py
│   │   ├── proposer.py
│   │   ├── skeptic.py                      # CHANGED: must address responder_risk
│   │   ├── single_agent_annotator.py
│   │   ├── orchestrator.py
│   │   ├── narrative_generator.py
│   │   └── fabrication_checker.py
│   │
│   ├── bus/  · tools/  · utils/
│
├── results/
│   ├── scores/  ·  signature/               # NEW
│   ├── filter_report.json                   # NEW — what each filter removed, and why
│   ├── responder_report.md                  # NEW — signature + flagged genes
│   ├── candidate_dossiers/
│   ├── ranked_candidates_positive.csv       # NEW — per direction
│   ├── ranked_candidates_negative.csv       # NEW
│   ├── validation_report.md
│   └── audit_log/
│
└── tests/
    ├── test_dataset_registrar.py            # NEW — a response_signature CANNOT rank
    ├── test_responder_signature.py          # NEW
    ├── test_filter_bank.py                  # NEW
    ├── test_mageck_runner.py                # NEW — batch + bidirectional
    ├── test_agents.py  ·  test_fabrication_checker.py
    └── fixtures/

# DELETED from v5:
#   workers/rra_aggregator.py
#   workers/screen_comparability_auditor.py
#   workers/mofa_integrator.py
#   tests/test_rra.py, tests/test_rra_r_parity.py, tests/test_comparability_auditor.py
```

### Dependencies — removed

```diff
- "pyrra>=0.1",        # RRA — no longer used
- "mofapy2>=0.7",      # MOFA+ — removed
```

Retained: pandas, numpy, scipy, scikit-learn, shap, `pulearn>=0.3` (correct package —
provides ElkanotoPuClassifier + bagging variant), networkx, openai, pydantic, biopython.

Plus MAGeCK (external binary — `mageck mle` / `mageck test`).

---

## 11. CLI

```bash
# Validate registration classes — a response_signature must not be registered as a screen
chroma validate --config chroma_config.yaml

# Layer 1: analyze the primary screen (batch-aware, bidirectional)
chroma screen --config chroma_config.yaml --direction both

# Layer 2: build the cholesterol-response signature and flag responders  [NEW]
chroma signature --config chroma_config.yaml

# Layer 3: run the filter bank
chroma filter --config chroma_config.yaml

# Layers 1-3 + prioritization, no LLM cost  (recommended first run)
chroma run --config chroma_config.yaml --workers-only

# Knowledge bank, ablation, agents  (unchanged from v5)
chroma run --config chroma_config.yaml --knowledge-bank-only
chroma ablate-debate --config chroma_config.yaml --sample-size 40
chroma run --config chroma_config.yaml --agents-only

chroma report --config chroma_config.yaml
chroma review --config chroma_config.yaml
```

---

## 12. Validation and Acceptance Criteria

The system is acceptable when all of the following hold:

1. **Registration is enforced.** No `response_signature` or `compartment_feature` dataset
   contributes a ranked gene list. `test_dataset_registrar.py` passes. **NEW**

2. **The primary screen is analyzed bidirectionally.** Both positive and negative
   regulator rankings are produced, and are never merged. **NEW**

3. **Rounds asymmetry is handled.** `high vs low` is the primary contrast. Every gene
   sourced from a `vs main` contrast carries `rounds_asymmetry_risk` and has been
   cross-checked against DepMap. **NEW**

4. **Batch structure is respected.** MAGeCK MLE run with batch covariate; every gene
   reports `batch_support`. No gene reaches P1/P2 on S2/S3-only support without that being
   surfaced. **NEW**

5. **The responder signature is built from clean contrasts only.** TRIM28-KO and myotube
   RNA-seq are excluded, and the exclusions are logged. `barred_from_evidence_layer` is
   enforced for TRIM28. **NEW — this is the criterion that prevents the §0 failure.**

6. **Canonical genes recover.** LDLR, LDLRAP1, NPC1, SCAP, MBTPS1 appear in P1 of the
   positive-regulator table, with concern flags attached, never removed. MYLIP/IDOL is
   expected in the negative-regulator table.

7. **LOPO works across ≥3 independently withheld modules**, reported as a distribution
   (median, IQR), not a single run.

8. **Label-shuffle collapses to chance.** No leakage.

9. **Precision@50 is reported and tracked.** No AUROC gate.

10. **No silent exclusion.** Every down-weighted gene is recoverable from the audit log
    with its reason.

11. **Every filter reports what it removed.** `filter_report.json` names, for each filter,
    the genes it removed and the false-positive class it was removing. **NEW**

12. **Fabrication rate < 2%** on a re-queried sample, for whichever annotation mode ships.

13. **The Skeptic addresses every `downstream_responder_risk: HIGH` candidate explicitly.**
    **NEW**

14. **At least one genuine `discovery_candidate` survives to a testable hypothesis.** A
    gene with no prior peer-reviewed characterization, from this lab's own screen, reaching
    P2 or P3, with `downstream_responder_risk: LOW`, a suggested validation experiment, and
    correct compartment concordance.

    **This is the only bar for the project's stated goal.** With the public benchmark path
    removed, there is no other.

15. **Computational core is bit-reproducible.** LLM layer archived and advisory.

---

## 13. Phased Rollout

**Phase 0 — Registration and screen QC.**
Write `.desc.md` for every dataset with `registration_class`, `phenotype_axis`,
`perturbation_class`, `replicate_independence`. Get Minyoung to sign off on the ddFP file
(`documentation_provenance: contemporaneous`). Resolve open TODOs: AmphoB units, Golgi-IP
replicate count and quantification method, Dylan's timings.
**Exit gate:** `chroma validate` passes; every dataset's class is human-asserted and
reviewed.

**Phase 1 — Primary screen, done properly.**
`chroma screen --direction both`. Batch-aware MAGeCK MLE. Both directions. All three
contrasts, with `high vs low` primary.
**Exit gate:** canonical genes (LDLR, NPC1, SCAP) recover in the positive-regulator
ranking; MYLIP-class genes appear in the negative ranking; the `high` arm produces a
biologically sane list, not noise. *If the `high` arm is pure noise, say so and drop it —
do not force it.*

**Phase 2 — The responder signature.**
`chroma signature`. Build from ORP9/GRAMD1/QKO. Confirm it is SREBP2-dominated (it should
be — that is the validation that it is measuring what it claims).
**Exit gate:** signature recovers known SREBP2 targets; TRIM28 and myotube data confirmed
excluded; `responder_report.md` written.

**Phase 3 — Filter bank + prioritization.**
`chroma run --workers-only`. Zero LLM cost.
**Exit gate:** LOPO distribution across ≥3 modules is defensible, not merely favorable on
one run. `filter_report.json` shows what each filter removed. The P2 bin is populated with
`downstream_responder_risk: LOW` genes.

**Phase 4 — Knowledge bank, ablation, agents.**
As v5 §16 Phases 2–4. Unchanged.

**Phase 5 — Human review and wet-lab hand-off.**
Review P2/P3 `discovery_candidate` dossiers. Adjudicate P5 (suspected feedback nodes) —
some will be real. Select ≥1 candidate for a concrete validation experiment.
**This is the phase that answers the project's question.** Everything prior is
infrastructure in service of getting here without missing a real signal or trusting a
fabricated one.

---

## 14. Explicit Assumptions

- The ddFP screen is the primary and effectively only ranking evidence. **The project's
  power is bounded by that one screen.** This is stated, not hidden.
- Effective replication in the primary screen is n ≈ 2 (S1 | S2+S3), not 3.
- The rounds asymmetry (`main` = 1 round, `high`/`low` = 3) is **unfixable in the data**.
  No round-3 `main` exists.
- Dylan's screen has no biological replication and cannot detect negative regulators. It
  corroborates; it does not rank.
- The Golgi-IP measures a different compartment (Golgi) than the screens (PM). It is a
  prior.
- The RNA-seq datasets measure *response*, not *regulation*. This is the foundational
  assumption of Layer 2. If it is wrong, Layer 2 is wrong.
- No public screen is phenotype- and format-compatible with the ddFP screen. Audited:
  Pfeffer/Lu (casTLE, permeabilized/lysosomal, sign-discordant), Ma ALOD4 (gene-trap
  hypomorphs, no fixed gene universe), Trinh (LDL uptake), van den Boomen (SREBP reporter),
  Chu (amphotericin — closest, worth re-examining).
- HeLa only. Cell-line artifacts are a flag, not a veto.
- Co-essentiality reflects cancer cell-line fitness — a known limitation, carried into
  every dossier.
- System designed for human protein-coding genes (~20,598, HGNC).

---

## 15. Open Questions

1. **Is the `high` arm real?** It has never been analyzed. It may be noise (0.7% gate,
   three rounds, small population). Phase 1 answers this. If it is noise, v6's
   bidirectionality collapses to v5's single axis — and that is a legitimate finding, not a
   failure.

2. **How strong is the responder confound, actually?** Layer 2 assumes SREBP2 targets will
   contaminate the screen hit list. This is well-motivated but untested *in this screen*.
   Phase 2 should quantify it: what fraction of the top 200 are signature members? If the
   answer is ~2%, Layer 2 is insurance. If it is ~40%, Layer 2 is the paper.

3. **Should Dylan's screen be re-run with replicates?** It is currently a filter because it
   cannot support more. Three biological replicates would make it a genuine second screen
   and would restore a real convergence layer. **This is the single highest-value new
   experiment available.** Estimated cost is low relative to what it buys.

4. **Chu et al. 2015 (amphotericin, PM cholesterol delivery)** is the closest public screen
   in *design* to Dylan's. It was not examined in detail. Worth a look — it may be the one
   public dataset that is genuinely comparable, and would give Dylan's arm a partner.

5. **Is TRIM28 a regulator or a feedback node?** It came from the screen; it has its own
   RNA-seq; it is the obvious first test case for Layer 2. Note the circularity bar: the
   TRIM28 RNA-seq cannot be used to *evaluate* TRIM28's own responder status without
   care.

6. **Does the Golgi-IP separate cholesterol from PI4P?** QKO raises both. The mevastatin
   arm partially deconvolves (cholesterol down, PI4P slightly up), but no condition cleanly
   isolates cholesterol. Any Golgi prior derived from this carries that ambiguity.

---

## 16. What v6 deliberately does not do

It does not rescue the convergence architecture by lowering the bar for what counts as a
screen. The four RNA-seq datasets *could* be made to produce ranked gene lists. Those lists
would aggregate beautifully. They would converge on SREBP2 targets, every metric would
pass, and the result would be worthless.

v6 removes the convergence layer because the data does not support it, and replaces it with
the discrimination layer that the data *does* support — and which, as it happens, addresses
a harder and more interesting problem than the one v5 was solving.

The honest version of this project is smaller and sharper than v5 imagined: **one carefully
analyzed bidirectional CRISPR screen, filtered by orthogonal evidence, with a principled
method for separating regulators from responders.**

That is publishable. The seven-dataset convergence engine was not, because the convergence
was not real.
