# Assessing the Impact of Parental Linear Gene Normalization on circRNA Differential Expression Analysis

This repository contains code and supplementary material associated with the article:

> **Assessing the impact of parental linear gene normalization in the performance of statistical models for circular RNA differential expression analysis**

The goal of this repository is to evaluate how different normalization strategies influence the performance of statistical models used for circRNA differential expression analysis, including the evaluation of using the linear gene information.

---

## Article

- Preprint / Published version: *to be added*

---

## Repository contents

- `Analysis_Scripts/` – data preprocessing, simulations, quality control, and execution of differential expression tools  
  - `Analysis_Scripts/Pre-processing/` – primary pipeline used to generate all circRNA and linear RNA count matrices analyzed in the manuscript
  - `Analysis_Scripts/DE-scripts/` - functions to run DE tools
- `supplementary/` – supplementary figures, tables, and additional analyses


---

## Analysis scope

The repository includes scripts for:

- Pre-processing circRNA and linear RNA count data
- Semi-parametric simulation of circRNA expression datasets
- Quality control of simulated datasets
- Running and benchmarking differential expression tools under different normalization strategies

This code is intended to support the analyses presented in the manuscript and is not designed as a general-purpose circRNA analysis pipeline.

---

## Reference genome

All analyses were performed using the **Homo sapiens primary assembly** reference genome from **Ensembl release 112**.

- Files available at:  
  https://ftp.ensembl.org/pub/release-112/

---

## Reproducibility

Due to data size constraints, raw sequencing data are not included in this repository. 
The used datasets are available in the NCBI database
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

For questions regarding this repository, please contact the corresponding author:

- X Y
- Email: 
---
## License

License information to be added.
