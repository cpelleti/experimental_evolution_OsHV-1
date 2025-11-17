# Experimentally mimicking long term infections of Magallana gigas with the OsHV-1 virus reveals evolution through positive selection

## Introduction
This project was developed to investigate the genomic evolution and selection dynamics of OsHV-1 µVar during experimental infections in oyster populations. The accompanying R scripts implement a complete workflow from statistical analysis of mortality and viral load across generations, to de novo genome assembly, variant calling, and detailed intra-individual Single Nucleotide Variant (iSNV) analyses. The code enables reproducible reconstruction of viral consensus genomes, identification of genomic variability, analysis of substitution patterns, principal component analyses for population structure, detection of selection signatures, and visualization of genotype frequencies using Muller plots. This workflow is designed for large-scale genomic datasets and allows researchers to track viral evolution across host populations and generations. All raw sequencing data are publicly available in NCBI SRA under Bioproject PRJNA1216400, This framework can be applied to study evolutionary dynamics, viral adaptation, and population genomics in other pathogen-host systems.

## Description
This repository contains scripts that have been used to performed analysis of experimental evolution datas from Pelletier et al, 2025.
Raw reads are available at http://www.ncbi.nlm.nih.gov/bioproject/PRJNA1216400. NCBI:BioProject: PRJNA1216400.

## Prerequisites
To run the analyses and reproduce the results presented in this project, the following prerequisites are required:

Software:
R v4.2.2 (or later) with the following packages: ggplot2 (v3.5.2 or later), rstatix (v0.7.2), factoextra (v1.0.7), ACER (v1.0), MullerPlot (v0.13), stats (base R), BCFtools (v1.17).
Geneious Prime v2022.2.2 for genome assembly and alignment visualization.
Fastp v0.20.1 for read quality filtering and trimming.
Bowtie2 v2.4.1 for read mapping.
metaSPAdes v3.15.4 for de novo genome assembly.
Freebayes v1.3.5 for variant calling.
MAFFT v1.4.0 for sequence alignment.
vt v2005.11.10 for variant normalization and deduplication.

Datasets:
Raw sequencing reads for 126 OsHV-1 µVar libraries, publicly available in NCBI SRA under Bioproject PRJNA1216400 (accession numbers SAMN46433918 to SAMN46434013).
Reference genome of Magallana gigas (GCF_902806645.1) for host read filtering.

Hardware Requirements:
A standard desktop or laptop computer with at least 16 GB RAM for genome assembly and variant calling.
Multi-core CPU recommended for parallelized analyses in metaSPAdes, Bowtie2, and Freebayes.

Usage Notes:
Scripts are designed to be run sequentially, starting from raw data preprocessing through to statistical analyses, variant calling, and visualization.

## Contents
| Folder/File   | Description |
| ------------- | ------------- |
| Content Cell  | Content Cell  |
| Content Cell  | Content Cell  |

File formats:
.R – R scripts, editable in RStudio or any text editor.
.sh.txt – Shell script for Linux/macOS terminal execution.
.csv or .tsv (in test_datas/) – data tables for input/output, compatible with R, Excel, or any text editor.

This structure allows sequential execution of the workflow from variant calling to statistical analyses and visualization of results.

## Contact
For help with this bioinformatics project, please reach out to Camille Pelletier (camille.pelletier@ifremer.fr), Maude Jacquot (maude.jacquot@ifremer.fr), Germain Chevignon (germain.chevignon@ifremer.fr), or Benjamin Morga (benjamin.morga@ifremer.fr). They can assist with any questions, troubleshooting, or clarifications. Feel free to contact them directly or through the appropriate communication channels.

## Authors and acknowledgment
Authors : Camille Pelletier, Nicole Faury, Mickaël Mege, Lionel Dégremont, Maelle Hattinguais, Jeremie Vidal-Dupiol, Germain Chevignon, Maude Jacquot, Benjamin Morga.

We thank the staff of the Ifremer station at La Tremblade (ASIM), Frédéric Girardin and his team (Plateforme des Mollusques Marins de La Tremblade PMMLT), and Virginie François and his team (Plateforme des Mollusques Marins de Bouin PMMB). We also thank the SEBIMER team for maintaining bioinformatics tools.
C.P.  was financially supported by a grant from the Ifremer Scientific Board and the Nouvelle-Aquitaine region. The present study was supported by scientific direction of Ifremer in context of the GT Huître project and by the FEAMP Gestinnov (FFEA470020FA1000007) GESTINNOV project and by DGAL (French General Directorate for Food) through the National Reference Laboratory for Mollusc Diseases and by the European-Union Reference Laboratory for Mollusc Diseases, Ifremer, La Tremblade. The authors acknowledge the Pôle de Calcul et de Données Marines (PCDM;  https://wwz.ifremer.fr/en/Research-Technology/Research-Infrastructures/Digital-infrastructures/Computation-Centre) for providing DATARMOR computing and storage resources.
C.P. and co-authors thank J.M.E. for the discussions and exchanges. We would like to thank Bruno Petton and his team PHYTNESS for sampling the virus in the Brittany area, as well as Johan Vieira and the CAPENA team for sampling OsHV-1 in Arcachon Bay. 