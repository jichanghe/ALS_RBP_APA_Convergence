# ALS RBP Alternative Polyadenylation Convergence

Analysis code and processed data supporting the manuscript:

**Partial Convergence of Alternative Polyadenylation Changes Across ALS-Associated RNA-Binding Proteins**

Author: **Changhe Ji**

## Overview

This repository contains R scripts and processed APA datasets used to analyze alternative polyadenylation (APA) changes across five ALS-associated RNA-binding protein perturbation datasets:

- TDP-43 knockdown (KD)
- TDP-43 K263E
- FUS knockout (KO)
- FUS P525L
- MATR3 knockout (KO)

The repository includes individual-dataset analyses, cross-dataset comparisons, Gene Ontology analyses, and quantitative QC/revision utilities.

## Repository structure

- `01_TDP43_KD/` — TDP-43 knockdown APA and GO analyses
- `02_TDP43_K263E/` — TDP-43 K263E APA and GO analyses
- `03_FUS_KO/` — FUS knockout APA analyses
- `04_FUS_P525L/` — FUS P525L APA analyses
- `05_MATR3_KO/` — MATR3 knockout APA analyses
- `06_Cross_dataset/` — cross-dataset integration and GO analyses
- `07_QC_revision/` — quantitative QC and revision utilities
- `data/` — processed APA data used by cross-dataset analyses
- `results/` — analysis outputs

## Processed data

The repository includes five processed differential APA datasets:

- `data/TDP43_KD_differential_APA.xlsx`
- `data/TDP43_K263E_differential_APA.xlsx`
- `data/FUS_KO_differential_APA.xlsx`
- `data/FUS_P525L_differential_APA.xlsx`
- `data/MATR3_KO_differential_APA.xlsx`

Processed mechanistic APA tables used for revision/QC analyses are provided under `data/QC_revision/`.

Raw FASTQ and BAM files are not redistributed. Raw RNA-sequencing datasets are publicly available from the original studies, with accession information provided in Table 1 and Supplementary Table 3 of the associated manuscript.

## Running the analyses

Run scripts from the repository root.

Example cross-dataset analysis:

`Rscript 06_Cross_dataset/run_5model_APA_integration_PUBLICATION_with_PanelF.R`

Example recurrent GO analysis:

`Rscript 06_Cross_dataset/run_5model_APA_integration_PUBLICATION_recurrent_GO.R`

QC audit:

`Rscript 07_QC_revision/audit_revision_statistics.R`

Reviewer-ready quantitative summary:

`Rscript 07_QC_revision/make_reviewer_ready_summary.R`

Analysis outputs are written under `results/`.

## Reproducibility notes

The TDP-43 depletion dataset represents **TDP-43 knockdown (KD)**. Some historical script filenames and internal variable names retain the label `TDP43_KO` because those names were used during the original analysis; they do not indicate a TDP-43 knockout experiment.

The processed differential APA tables included in this repository support the cross-dataset analyses. Some individual-dataset scripts depend on upstream QAPA or quantification outputs generated from the original public RNA-sequencing datasets; these large intermediate files are not duplicated here.

`07_QC_revision/collect_revision_statistics.R` is retained as a legacy discovery/QC utility. It recursively scans a complete upstream QAPA analysis tree and can be run as:

`Rscript 07_QC_revision/collect_revision_statistics.R /path/to/QAPA_root`

The archived `results/QC_revision/01_revision_quantitative_summary.csv` was generated from the complete original analysis tree. Some upstream PAU matrices recorded in its provenance fields are not included in this compact repository.

The R environment used to prepare this archive is recorded in `sessionInfo.txt`.

## Citation

If you use this repository, please cite the associated manuscript:

**Ji C. Partial Convergence of Alternative Polyadenylation Changes Across ALS-Associated RNA-Binding Proteins.**

A permanent Zenodo DOI will be added after archival of the versioned GitHub release.

## License

This repository is distributed under the MIT License. See `LICENSE`.
