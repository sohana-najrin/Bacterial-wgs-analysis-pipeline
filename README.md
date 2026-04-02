# 🧬 Bacterial Whole Genome Sequencing (BacWGS) Pipeline

## 📌 Overview
This repository provides an end-to-end pipeline for bacterial genome analysis including:

- Quality control
- Read trimming
- Genome assembly
- Genome annotation
- Taxonomic classification
- Phylogenetic analysis

---

## ⚙️ Tools Used
- FastQC
- fastp
- SPAdes
- QUAST
- Prokka
- Kraken2
- IQ-TREE / FastTree

---

## 🔄 Workflow

1. Quality Control (FastQC)
2. Trimming (fastp)
3. Genome Assembly (SPAdes)
4. Assembly Evaluation (QUAST)
5. Annotation (Prokka)
6. Taxonomic Classification
7. Phylogenetic Tree Construction

---

## ▶️ Usage

```bash
bash workflows/bacwgs_pipeline.sh sample1 
