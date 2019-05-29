# Ecological significance of microbial transporters in the Baltic Sea

[TOC]

## Overview
Here are scripts and notebooks used in the analysis of transporters for the article
[insert ref here](http://).

## Material and Methods

### Datasets
A metagenomic reference assembly was generated as described in
[BARM and BalticMicrobeDB, a reference metagenome and interface to meta-omic data for the Baltic Sea](https://www.nature.com/articles/sdata2018146).
#### Metagenomes
The metagenomic and sequence data was preprocessed
and mapped onto the reference assembly to quantify genes as outlined
in that paper.

#### Metatranscriptomes
##### Preprocessing
The 25 metatranscriptomic samples were preprocessed as follows:

1. Adapter trimming using [Cutadapt](https://github.com/marcelm/cutadapt) (v. 1.8.0) with the ILLUMINA universal adapter and default settings
2. Quality trimming using [Sickle](https://github.com/najoshi/sickle) (v. 1.210) with default settings
3. Read sorting (rRNA/non_rRNA) using [SortMeRNA](https://github.com/biocore/sortmerna) (v. 2.0) with default settings and all RNA databases

### Notebooks

All tables, figures and results presented in the article were obtained using
the notebooks in this repository. They were executed in the following order:

1. [01.process_data.ipynb](process_data.ipynb)
2. [02.cluster_samples.ipynb](02.cluster_samples.ipynb)
3. [03.transporter_stats_taxprofiles.ipynb](03.transporter_stats_taxprofiles.ipynb)
4. [04.diversity_and_example_categories.ipynb](diversity_and_example_categories.ipynb)
5. [05.process_and_run_deseq2.ipynb](05.process_and_run_deseq2.ipynb)
6. [06.plot_deseq2.ipynb](06.plot_deseq2.ipynb)
7. [07.environmental_correlations.ipynb](07.environmental_correlations.ipynb)
8. [08.environmental_plot.ipynb](08.environmental_plot.ipynb)