# RegSeqDB

A structured, queryable database linking promoter sequence variation to gene expression and transcription factor (TF) and RNA polymerase (RNAP) binding in *E. coli*.

## Overview

RegSeqDB serves as the primary database for data related to a [preprint by Röschinger et al.](https://arxiv.org/abs/2505.08764) and the subsequent analysis performed by the Galagan Lab. It provides easy access to promoter sequences, RNA/DNA barcode counts across experimental conditions, and TF and RNAP binding sites and energies — to support reliable and repeatable downstream analyses of gene regulation in *E. coli*.

## Features

- **Single Promoter Search** — Query a promoter and transcription factor pair under any of 39 experimental conditions to retrieve expression and binding energy data with interactive Plotly visualizations.
- **Compare Between Conditions** — Compare expression and TF binding across two experimental conditions using log fold change regression plots.
- **Genomic Locus Viewer** — Visualize TSS, RNAP, and TF binding site positions for a given promoter.
- **CSV Export** — Download full search results for downstream analysis.

## Data

RegSeqDB contains:
- 96K+ promoter sequence variants
- RNA/DNA barcode counts across 39 experimental conditions
- TF binding site predictions and affinity scores from [BoltzNet](https://www.nature.com/articles/s41467-025-58862-8)
- RNAP binding site positions and energies from [Promoter Calculator](https://www.nature.com/articles/s41467-022-32829-5)

## Experimental Background

Researchers in the Phillips Lab sought to identify condition-dependent regulatory elements within *E. coli* promoter sequences by performing Reg-Seq — a method that introduces random mutations into promoters and measures the effects on downstream gene expression under varying experimental conditions. The Galagan Lab then applied [BoltzNet](https://www.nature.com/articles/s41467-025-58862-8) and [Promoter Calculator](https://www.nature.com/articles/s41467-022-32829-5) to predict TF and RNAP binding sites and their corresponding binding energies across all sequence variants, with the goal of investigating their value in predicting gene expression. Current analyses are exploring the use of these biophysical quantities to build predictive models of transcription. Considering that future analyses are likely to grow in complexity as additional models for transcription factor binding are integrated into BoltzNet, a structured and easily expandable database is imperative in supporting these long-term research objectives.

## References

1. Röschinger et al. [Condition-dependent regulatory elements in E. coli promoters identified by Reg-Seq.](https://arxiv.org/abs/2505.08764) Preprint.
2. BoltzNet — [Galagan Lab biophysical model for transcription factor binding prediction.](https://www.nature.com/articles/s41467-025-58862-8)
3. LaFleur, Hossain & Salis. [Automated model-predictive design of synthetic promoters to control transcriptional profiles in bacteria.](https://www.nature.com/articles/s41467-022-32829-5) *Nature Communications*, 2022.

## Team

Developed by Bradley Jenner, Zokhira Mukhammadyunusova, Mohammad Gharandouq, and Leah Morzenti under the supervision of Dr. James Galagan at Boston University, as part of BF768, Spring 2026, instructed by Dr. Gary Benson.

## Tech Stack

- **Backend** — Python, Flask, MariaDB
- **Frontend** — HTML, CSS, JavaScript, jQuery, Plotly.js
- **Hosted** — Boston University BioEd server
