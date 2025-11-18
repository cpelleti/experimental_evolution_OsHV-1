# Experimentally mimicking long term infections of Magallana gigas with the OsHV-1 virus reveals evolution through positive selection

## Introduction
This project was developed to investigate the genomic evolution and selection dynamics of OsHV-1 µVar during experimental infections in oyster populations. The accompanying R scripts implement a complete workflow from statistical analysis of mortality and viral load across generations, to de novo genome assembly, variant calling, and detailed intra-individual Single Nucleotide Variant (iSNV) analyses. The code enables reproducible reconstruction of viral consensus genomes, identification of genomic variability, analysis of substitution patterns, principal component analyses for population structure, detection of selection signatures, and visualization of genotype frequencies using Muller plots. This workflow is designed for large-scale genomic datasets and allows researchers to track viral evolution across host populations and generations. All raw sequencing data are publicly available in NCBI SRA under Bioproject PRJNA1216400, This framework can be applied to study evolutionary dynamics, viral adaptation, and population genomics in other pathogen-host systems.

## Description
This repository contains scripts that have been used to performed analysis of experimental evolution datas from Pelletier et al, 2025.\
Raw reads are available at http://www.ncbi.nlm.nih.gov/bioproject/PRJNA1216400. NCBI:BioProject: PRJNA1216400.

## Prerequisites
To run the analyses and reproduce the results presented in this project, the following prerequisites are required:

**Software requirements for R:**\
All analyses were performed in R version 4.2.2. The following R packages are needed to run the included scripts: 

→ <ins>Available</ins> on CRAN for current version of R 
These packages can be installed by running the command `install.packages(package_name)` 

1. tidyverse (version 2.0.0): Load multiple tidyverse packages in a single step, including dplyr (version 1.1.4), ggplot2 (3.4.0), tibble (version 3.3.0), stringr (version 1.5.1), and tidyr (version 1.3.1)
2. factoextra (version 1.0.7): Extract and visualize the output of multivariate data analyses 
3. NbClust (version 3.0.1): Determine the best number of clusters 
4. MullerPlot (version 0.13): Visualize emergence and loss of genotypes 
5. ggpubr (version 0.6.1): Additional plotting functionality for ggplot2 
6. gggenes (version 0.5.1): Drawing gene arrow maps in ggplot2 
7. rstatix (version 0.7.2): Framework for performing common statistical tests, computing effect sizes
8. dplyr (version 1.1.4): Data manipulation, providing tools to efficiently filter, transform, summarize, and organize data
9. knitr (version 1.50): Enables dynamic report generation by integrating R code, results, and narrative text into documents
10. stats (base R, version 4.5.1): Provides fundamental statistical functions, including probability distributions, hypothesis testing, linear and nonlinear modeling 


→ <ins>Not available</ins> on CRAN for current version of R 

1. poolSeq 
```
install.packages("remotes")
 # Install PoolSeq from GitHub
remotes::install_github("ThomasTaus/poolSeq") 
 ```
 Or follow the installation instructions provided on the GitHub repository: https://github.com/ThomasTaus/poolSeq/. Ensure all dependencies are installed first: data.table (>= 1.9.4), foreach (>= 1.4.2), stringi (>= 0.4-1), matrixStats (>= 0.14.2), Rcpp (>= 1.0.3) 
 
 2. ACER v1.0 
 ```
 # Install ACER from GitHub
remotes::install_github("MartaPelizzola/ACER")
```
Or follow the installation instructions provided on the GitHub repository: https://github.com/MartaPelizzola/ACER


**Software requirements for bash environment:**\
Fastp v0.20.1 for read quality filtering and trimming.\
Bowtie2 v2.4.1 for read mapping.\
metaSPAdes v3.15.4 for de novo genome assembly.\
Freebayes v1.3.5 for variant calling.\
MAFFT v1.4.0 for sequence alignment.\
vt v2005.11.10 for variant normalization and deduplication.

**Software requirements for analyses:**\
Geneious Prime v2022.2.2 for genome assembly and alignment visualization.


**Datasets:**\
Raw sequencing reads for 126 OsHV-1 µVar libraries, publicly available in NCBI SRA under Bioproject PRJNA1216400 (accession numbers SAMN46433918 to SAMN46434013).\
Reference genome of Magallana gigas (GCF_902806645.1) for host read filtering.

**Hardware Requirements:**\
A standard desktop or laptop computer with at least 16 GB RAM for genome assembly and variant calling.\
Multi-core CPU recommended for parallelized analyses in metaSPAdes, Bowtie2, and Freebayes.

**Usage Notes:**\
Scripts are designed to be run sequentially, starting from raw data preprocessing through to statistical analyses, variant calling, and visualization.

## Contents
| Folder/File   | Description |
| ------------- | ------------- |
| test_datas/   | Contains example datasets for testing the scripts and reproducing the analyses.  |
| 01_variant_calling.sh  | Shell script for variant calling pipeline, including read mapping and Freebayes variant detection.  |
| 02_stats_mortalities_viralLoad_Fig3.R   | R script for statistical analysis of mortality and viral load across oyster generations and populations, and generation of Figure 3.  |
| 03_stats_ts_tv_Fig5.R  | R script to analyze transition/transversion ratios, substitution patterns, and generate Figure 5.  |
| 04_Genomic_selection_detection_Fig6-7-8.R   | R script for detecting selection signatures across viral generations and populations, generating Figures 6, 7, and 8.  |
| 05_allelic_trajectories_Fig9.R  | R script to compute allelic frequency trajectories and generate Figure 9. |
| 06_muller_plot_Fig10.R  | R script to visualize genotype emergence and lineage dynamics using Muller plots (Figure 10). |

**File formats:**\
&nbsp;&nbsp;&nbsp;&nbsp;.R – R scripts, editable in RStudio or any text editor.\
&nbsp;&nbsp;&nbsp;&nbsp;.sh.txt – Shell script for Linux/macOS terminal execution.\
&nbsp;&nbsp;&nbsp;&nbsp;.csv or .tsv (in test_datas/) – data tables for input/output, compatible with R, Excel, or any text editor.

This structure allows sequential execution of the workflow from variant calling to statistical analyses and visualization of results.

## Usage
This repository provides scripts and pipelines for the statistical and genomic analysis of OsHV-1 in oyster populations, including mortality, viral load, viral genome assembly, variant calling, selection detection, and lineage tracking. Below are step-by-step instructions for using the code.

**1. Statistical analysis of mortality and viral load**

Input: Mortality and viral load data across oyster populations and generations.\
Execution: script 02_stats_mortalities_viralLoad_Fig3.R\
Output: Graphical representations of mortality and viral load; ANOVA tables showing differences among populations and generations.

**2. Assembly of viral NR-genomes**

Input: Raw Illumina paired-end reads from OsHV-1 infected oyster tissues.\
Execution: Quality filter and trim reads using Fastp ; Remove host sequences by mapping reads to the M. gigas reference genome using Bowtie2; Assemble OsHV-1 reads de novo with metaSPAdes ; Identify OsHV-1 contigs with BLASTn and reorder in Geneious to reconstruct NR-genomes.\
Output: Non-redundant viral genomes (NR-genomes) for each viral suspension.

**3. Variant calling**

Input: Filtered reads and the artificial OsHV-1 major NR-consensus._
Execution: script 01_variant_calling.sh\
Post-processing: Normalize, decompose, sort, and deduplicate variants using vt, then retain variants with ≥10% frequency.\
Output: VCF files containing intra-individual Single Nucleotide Variants (iSNVs).

**4. Substitution analysis and clustering**

Input: VCF files containing filtered intra-individual Single Nucleotide Variants (iSNVs).\
Execution: script 03_stats_ts_tv_Fig5.R\
Analyze transitions/transversions with BCFtools stats ; Average substitution counts per generation and population ; Visualize using ggplot2 ; Perform PCA on iSNV frequencies to explore sample clustering using factoextra in R.\
Output: Plots showing substitution patterns, ts/tv ratios, and sample clustering across oyster populations.

**5. Selection detection**

Input: VCF files containing filtered intra-individual Single Nucleotide Variants (iSNVs).\
Execution: scripts 04_Genomic_selection_detection_Fig6-7-8.R ; 05_CMH.R and 06_allelic_trajectories_Fig9.R\
Compare allele frequencies between initial and evolved generations ; Apply Cochran-Mantel-Haenszel (CMH) tests using ACER in R ; Adjust p-values using FDR and visualize significant iSNVs with Manhattan plots.\
Output: iSNVs under selection and their allele frequency trajectories across generations.

**6. Genotype proportion visualization**

Input: VCF files containing filtered intra-individual Single Nucleotide Variants (iSNVs).\
Execution: script 07_muller_plot_Fig10.R\
Assign iSNVs to their origin based on G0 viral inoculum ; Group variants by geographic origin and visualize relative abundances using MullerPlot in R.\
Output: Muller plots showing lineage emergence and dynamics throughout the experimental evolution.

## License
This work is licensed under the MIT License.

## Contact
For help with this bioinformatics project, please reach out to Camille Pelletier (camille.pelletier@ifremer.fr), Maude Jacquot (maude.jacquot@ifremer.fr), Germain Chevignon (germain.chevignon@ifremer.fr), or Benjamin Morga (benjamin.morga@ifremer.fr). They can assist with any questions, troubleshooting, or clarifications. Feel free to contact them directly or through the appropriate communication channels.

## Authors and acknowledgment
Authors : Camille Pelletier, Nicole Faury, Mickaël Mege, Lionel Dégremont, Maelle Hattinguais, Jeremie Vidal-Dupiol, Germain Chevignon, Maude Jacquot, Benjamin Morga.

We thank the staff of the Ifremer station at La Tremblade (ASIM), Frédéric Girardin and his team (Plateforme des Mollusques Marins de La Tremblade PMMLT), and Virginie François and his team (Plateforme des Mollusques Marins de Bouin PMMB). We also thank the SEBIMER team for maintaining bioinformatics tools.\
C.P.  was financially supported by a grant from the Ifremer Scientific Board and the Nouvelle-Aquitaine region. The present study was supported by scientific direction of Ifremer in context of the GT Huître project and by the FEAMP Gestinnov (FFEA470020FA1000007) GESTINNOV project and by DGAL (French General Directorate for Food) through the National Reference Laboratory for Mollusc Diseases and by the European-Union Reference Laboratory for Mollusc Diseases, Ifremer, La Tremblade. The authors acknowledge the Pôle de Calcul et de Données Marines (PCDM;  https://wwz.ifremer.fr/en/Research-Technology/Research-Infrastructures/Digital-infrastructures/Computation-Centre) for providing DATARMOR computing and storage resources.\
C.P. and co-authors thank J.M.E. for the discussions and exchanges. We would like to thank Bruno Petton and his team PHYTNESS for sampling the virus in the Brittany area, as well as Johan Vieira and the CAPENA team for sampling OsHV-1 in Arcachon Bay. 