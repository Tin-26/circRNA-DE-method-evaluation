# Assessing the Impact of Parental Linear Gene Normalization on circRNA Differential Expression Analysis

This repository contains code and supplementary material associated with the article:

> **Assessing the impact of parental linear gene normalization in the performance of statistical models for circular RNA differential expression analysis**

The goal of this repository is to evaluate how different normalization strategies influence the performance of statistical models used for circRNA differential expression analysis, including the evaluation of using the linear gene information.

---

## Article

- [Preprint](https://www.biorxiv.org/content/10.64898/2026.03.06.710045v1)
- Published version: *to be added*

---

## Code availability

This repository contains the analysis code used for preprocessing, simulation, differential expression benchmarking, and visualization in the manuscript.

Main scripts:

- `Analysis_Scripts/Pre-processing/Run_pipeline.py`  
  Runs the preprocessing workflow for FASTQ quality control, circRNA identification, and count-matrix generation.

- `Analysis_Scripts/Simulation/Simulate_datasets.R`  
  Generates semi-parametric simulated circRNA count datasets using SPsimSeq.

- `Analysis_Scripts/Simulation/countsimQC.R`  
  Performs quality control of simulated datasets.

- `Analysis_Scripts/DE_benchmark.R`  
  Main benchmarking script for limma-voom, edgeR, DESeq2, and ciriDE-style analyses.

- `Analysis_Scripts/Data_exploration/`  
  PCA, exploratory statistics, and filtering/outlier checks.

The code is released under the CC0 1.0 Universal license; see `LICENSE`.

---

## Software requirements

For the binaries of the circRNA identifier tools, please refer to the original developers’ GitHub pages:

- **CIRI3**: [https://github.com/gyjames/CIRI3](https://github.com/gyjames/CIRI3)
- **CIRCexplorer2**: [https://github.com/YangLab/CIRCexplorer2/tree/master](https://github.com/YangLab/CIRCexplorer2/tree/master)
- **CLEAR**: [https://github.com/YangLab/CLEAR](https://github.com/YangLab/CLEAR)
- **CircTools2**: [https://github.com/jakobilab/circtools](https://github.com/jakobilab/circtools)

R packages:

- edgeR
- limma
- DESeq2
- tidyverse
- SPsimSeq
- SimSeq

Python tools/scripts require Python 3 and external command-line tools used by the preprocessing pipeline.

External tools:

- STAR
- featureCounts/Subread
- fastp
- FastQC

---
## Data availability

Raw sequencing files are not stored in this repository because of file size constraints. Public datasets can be retrieved from NCBI/SRA using the accession identifiers below.

| Dataset label | Source/accession | Used for | Availability |
|---|---|---|---|
| BC | PRJNA553624 | real-data benchmarking and simulation source data | Public NCBI/SRA |
| HCC-tissue | PRJNA716508 | real-data benchmarking and simulation source data | Public NCBI/SRA |
| HCC-PBMC | PRJNA754685 | real-data benchmarking and simulation source data | Public NCBI/SRA |
| SCLC | PRJNA1237743 | real-data benchmarking and simulation source data | Public NCBI/SRA |
| EBC1 | internal cohort - PRJNA1429817  | real-data benchmarking and simulation source data | Public NCBI/SRA |
| EBC2 | internal cohort - PRJNA1429817  | real-data benchmarking and simulation source data | Public NCBI/SRA |

## Files not included in this repository

The following files are not included directly in Git due to size and/or privacy constraints:

- raw FASTQ files
- BAM/SAM alignment files
- large count matrices
- simulated replicate datasets
- full benchmark result tables

Where possible, derived non-identifying files are provided through the external data archive listed above.

---

## Reference genome

All analyses were performed using the **Homo sapiens primary assembly** reference genome from **Ensembl release 112**.

- Files available at:  
  https://ftp.ensembl.org/pub/release-112/

---

## Reproducibility

Due to data size constraints, raw sequencing data are not included in this repository. 
The used datasets are available in the NCBI database. The needed accession identifiers can be found in the published article.
Scripts required to reproduce simulations, normalization strategies, and benchmarking analyses are provided.

Details on data sources and execution order are described in the manuscript and supplementary materials.

---

## Citation

If you use this code or reuse parts of the analysis, please cite:

*Author et al.* (Year).  
Assessing the impact of parental linear gene normalization in the performance of statistical models for circular RNA differential expression analysis.  
*Journal*. DOI

---

## Contact

For questions regarding this repository, please contact the first author:

- Erda Qorri
- Email: erda.qorri@brc.hu

---

## License

This repository is made available under the **Creative Commons CC0 1.0 Universal (CC0 1.0)** license.

To the extent possible under law, the author(s) have waived all copyright and related or neighboring rights to the contents of this repository.

See the `LICENSE` file for the full legal text.
