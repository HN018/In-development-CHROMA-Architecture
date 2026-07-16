# CHROMA
### Cholesterol Homeostasis Regulator discovery via Orchestrated Multi-method Analysis

A self-initiated design study in **multi-omics data integration** for target discovery — and, as much as anything, a worked example of *stress-testing a computational design against real data and revising it when the evidence demands it.*

> **Status: architecture + design docs (in development).** This repository is a system design, not a deployed discovery pipeline. It documents both the intended architecture and — importantly — what changed when the design met real datasets.

---

## The problem

Pooled CRISPR screens, proteomics, and transcriptomics each carry partial, noisy evidence about which genes regulate plasma-membrane cholesterol. No single dataset is decisive. The original idea behind CHROMA was to integrate several such datasets so that genes supported across independent lines of evidence rise to the top — with a deterministic, reproducible core and a large-language-model layer confined to *hypothesis generation and annotation only*, never to selecting or excluding genes.

## Intended architecture (5 layers)

| 1 — Normalization | Harmonize heterogeneous score types into comparable, provenance-tagged ranks without silently dropping genes tested in only one modality | Designed |
| 2 — Integration | Combine ranked evidence across modalities into a single prioritization | Redesigned (see below) |
| 3 — Network context | Propagate support across a protein-interaction graph so genes linked to strong hits gain weight | Designed |
| 4 — Grounded LLM annotation | Advisory only: annotate uncharacterized ("orphan") candidates against existing literature; **no gene selection or exclusion** | Designed |
| 5 — Prioritization & reporting | Produce a ranked, evidence-attributed shortlist for experimental follow-up | Designed |

## What changed, and why

The original design assumed the available datasets were **replicated screens measuring the same thing**, so that convergence across them would be meaningful. A systematic audit of the actual datasets showed that assumption did not hold:

- Only one dataset qualified as a genuinely replicated functional-genomics screen.
- The transcriptomic datasets measure the cell's **response** to cholesterol perturbation (e.g. SREBP2 feedback targets), not the **regulators** of cholesterol — so aggregating them would have produced strong, statistically clean convergence on the *wrong* class of genes, while passing every acceptance check.

That is the dangerous failure mode: not a pipeline that breaks, but one that **succeeds convincingly at the wrong task.** The response to that was to stop, and redesign the integration around **distinguishing regulators from downstream responders** — treating response signatures as a *filter to remove false positives* rather than as positive evidence.

This is the part of the project I consider most valuable: the design is only as trustworthy as the data feeding it, and recognizing when it isn't is the actual work.

## Design principles

- **Deterministic core, advisory LLM.** Every gene-level decision is made by reproducible, inspectable computation. The LLM annotates and proposes; it never selects or excludes.
- **No silent exclusion.** Genes measured in only some datasets are handled explicitly, not dropped without trace.
- **Provenance-aware.** Every score carries where it came from, so integration never blends incomparable evidence unknowingly.
- **Benchmark ≠ discovery.** Method behavior is validated on controlled benchmarks separately from any real-data discovery claim.

## What this repository is (and isn't)

- **Is:** a documented architecture, the design rationale, and an honest record of a design-vs-data course correction.
- **Isn't:** a validated system that has discovered novel regulators. It hasn't — the available data could not yet support that claim, which is itself the finding that drove the redesign.


