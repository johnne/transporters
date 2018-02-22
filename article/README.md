# Ecological significance of microbial transporters in the Baltic Sea
John Sundh, Carlo Berg, Bärbel Müller-Karulis, Christofer M.G Karlsson, Anders F. Andersson, Johannes Alneberg, Jarone Pinhassi, Åke Hagström, Ulla Li Zweifel

[TOC]

## Overview
Here are scripts and notebooks used in the analysis of transporters for the article
[insert ref here](http://).

## Material and Methods

### Software
Install the conda environment.
```
conda env create -f env/environment.yml
conda activate transport-article
```

### Datasets
A metagenomic reference assembly was generated as described in
[BARM and BalticMicrobeDB: a reference metagenome and interface to meta-omic data for the Baltic Sea](https://www.nature.com/sdata/) (*in revision*).
Metatranscriptomic samples were mapped onto the reference assembly to quantify genes.

```
mkdir data
cd data
```

For the purpose of these analyses the following annotations and abundance values were used):

#### BARM reference annotations
Download the following:
* [BARM archive](https://drive.google.com/open?id=0B_prCMxfYyv7ZTRJSjJNNkl6ZGM)

Extract and compress the specific files needed:

```
tar -zxvf barm_files.tar.gz barm_files/annotations/all.TIGRFAM.standardized.tsv
gzip barm_files/annotations/all.TIGRFAM.standardized.tsv
tar -xzvf barm_files.tar.gz barm_files/annotations/taxonomy_per_gene.tsv
gzip barm_files.tar.gz barm_files/annotations/taxonomy_per_gene.tsv
```


#### Metagenomic mappings
```
mkdir mg
```

Download the following:
* [LMO gene TPM values](https://drive.google.com/open?id=0B_prCMxfYyv7Z1RXNHRFeFhRams)
* [LMO gene RAW counts](https://drive.google.com/open?id=0B_prCMxfYyv7LXA2TXlBMXYzcUU)

```
mv all_genes.tpm.tsv.gz mg/
mv all_genes.raw_counts.tsv.gz mg/
```

#### Metatranscriptomic mappings
##### Preprocessing
The 25 metatranscriptomic samples were preprocessed as follows:

1. Adapter trimming using [Cutadapt](https://github.com/marcelm/cutadapt) (v. 1.8.0) with the ILLUMINA universal adapter and default settings
2. Quality trimming using [Sickle](https://github.com/najoshi/sickle) (v. 1.210) with default settings
3. Read sorting (rRNA/non_rRNA) using [SortMeRNA](https://github.com/biocore/sortmerna) (v. 2.0) with default settings and all RNA databases

```
mkdir mt
```
Download the following:
* [LMO gene TPM values](https://)
* [LMO gene RAW counts](https://)

**TODO: Add info on how to obtain the RNA data**

### Notebooks

All tables, figures and results presented in the article were obtained using
the notebooks in this repository. They were executed in the following order:

1. [process_data.ipynb](process_data.ipynb)
2. [transporter_stats_taxprofiles.ipynb](transporter_stats_taxprofiles.ipynb)
3. [environmental_plot.ipynb](environmental_plot.ipynb)