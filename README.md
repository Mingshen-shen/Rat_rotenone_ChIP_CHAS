## Deconvolution analysis using CHAS
This subproject deconvolutes the bulk H3K27ac data using CHAS, a deconvolution tool. So that we can determine which DARs were specific to individual cell types. Furthermore, this approach investigated cell type–specific gene regulatory changes and dysregulated pathways across the substantia nigra (SN) using GO analysis. 

### Aims
The aims of this project are:
- To deconvolute a bulk H3K27ac ChIP-seq datasets using CHAS.
- To identify cell type–specific differentially acetylated regions (DARs) across the substantia nigra (SN).
- To investigate dysregulated gene regulatory pathways and biological processes through GO enrichment analysis.

## Data Access
The raw and processed data used in this project are stored on the King’s College London (KCL) High Performance Computing (HPC) system. Due to institutional data governance policies, these datasets are not publicly available.

Researchers who require access to the data for collaborative purposes may contact the project contributors listed below.

## Contributors
- **Bomin Lee** — Responsible for preprocessing and the analysis of Cell Type-Specific Regulatory Landscapes.
- **Mingshen Shen** — Responsible for preprocessing and Deconvolution analysis using CHAS.
- **Maria Tsalenchuk** — Provided the raw ChIP-seq datasets used in this project.
- **Paulina Urbanaviciute**, **Sarah J. Marzi** — Supervision and project guidance.

## Reference
This project makes use of a deconvolution algorithm called CHAS from the following GitHub repository:

🔗 [Kitty Murphy et al.](https://github.com/Marzi-lab/CHAS)

The raw data used for CHAS deconvolution were obtained from the following GitHub repository:

🔗 [Maria Tsalenchuk](https://github.com/Marzi-lab/rotenone_rat_ChIP)


