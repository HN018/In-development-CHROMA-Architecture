# CHROMA
**Cholesterol Homeostasis Regulator discovery via Orchestrated Multi-method Analysis**

*Status: Active Development / Architecture Phase*

## Overview
CHROMA is an in-development 5-layer multi-agent computational framework designed to surface novel metabolic and organelle regulators. 

Standard single-method screens (e.g., isolated CRISPR knockouts or standalone proteomics) frequently discard critical, sub-threshold biological targets due to rigid statistical cutoffs. CHROMA solves this "silent exclusion" problem by integrating heterogeneous datasets-CRISPR KO screens, RNA-seq, and proteomics-into a unified pipeline. It identifies genes that are consistently supported across orthogonal screens, even if they are statistically weak in any single assay, to discover truly novel regulators of cholesterol homeostasis and lysosomal integrity.

## Core Philosophy
1. **No Silent Exclusion:** "Dropouts" and hard cutoffs are removed. Every target maintains a continuous score and an audit trail to the end. "Exclusion" simply means being routed to a deprioritized, documented data bin.
2. **Convergence as Evidence:** The system rewards targets consistently supported across multiple screens via Robust Rank Aggregation (RRA) and multi-layer network propagation.
3. **Grounded AI:** Large Language Models (LLMs) are used strictly as an advisory annotation layer. LLM agents engage in an adversarial debate (Proposer vs. Skeptic) to generate mechanistic hypotheses, but all factual claims must be cited against explicit database identifiers (e.g., DepMap, BioGRID, PubMed).

## System Architecture

CHROMA operates across 5 deterministic and advisory layers:

### Layer 1: Evidence (Data Normalization)
Extracts raw assay data and normalizes it using **Fractional Rank Scoring** ($r_i/N_i$) via `pandas`. This resolves statistical skews between massive transcriptomic lists and smaller, high-dropout proteomic datasets without deleting missing proteins.

### Layer 2: Convergence (Integration)
The computational heart of the framework.
* **Robust Rank Aggregation (RRA):** Identifies genes consistently ranked better than chance.
* **Network Propagation:** Random Walk with Restart (RWR) over weighted molecular networks.
* **Co-Essentiality Proximity:** Maps targets against canonical sterol machinery and lysosomal/mitophagy structural anchors using DepMap data.

### Layer 3: Prioritization
Scores candidates using a Positive-Unlabeled (PU) classifier. Progress is validated strictly via Leave-one-pathway-out (LOPO) testing to ensure the model discovers novel biology rather than interpolating known gene clusters.

### Layer 4: Annotation (LLM Adversarial Debate)
A multi-agent AI layer featuring a **Triage Filter** to optimize API costs and prevent hallucination. 
* Highly canonical targets bypass the debate. 
* Novel and orphan targets trigger a 3-round debate between a **Proposer** (arguing biological relevance) and a **Skeptic** (querying databases for artifact/essentiality risks), moderated by an **Orchestrator**. 

### Layer 5: Triage & Output
Outputs highly structured candidate dossiers routed into distinct biological bins (e.g., *Convergent-Novel*, *Orphan/Essential*, *Convergent-Canonical*).
