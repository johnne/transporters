# Ecological significance of microbial transporters in the Baltic Sea
John Sundh, Carlo Berg, Bärbel Müller-Karulis, Christofer M.G Karlsson, Anders F. Andersson, Johannes Alneberg, Jarone Pinhassi, Åke Hagström, Ulla Li Zweifel

## Overview
Here are scripts and notebooks used in the analysis of transporters for the article
[insert ref here](http://).

## Material and Methods
### Datasets
A metagenomic reference assembly was generated as described in
[BARM and BalticMicrobeDB: a reference metagenome and interface to meta-omic data for the Baltic Sea](https://www.nature.com/sdata/) (*in revision*).
Metatranscriptomic samples were mapped onto the reference assembly to quantify genes.

```
mkdir data
cd data
```

For the purpose of these analyses the following annotations and abundance values were used):

**Annotations**

* [BARM archive](https://drive.google.com/open?id=0B_prCMxfYyv7ZTRJSjJNNkl6ZGM)

```
tar -zxvf barm_files.tar.gz barm_files/annotations/all.TIGRFAM.standardized.tsv
```


**Metagenomic abundances**
```
mkdir mg
```

Download the following:
* [LMO gene TPM values](https://drive.google.com/open?id=0B_prCMxfYyv7Z1RXNHRFeFhRams)
* [LMO gene RAW counts](https://drive.google.com/open?id=0B_prCMxfYyv7LXA2TXlBMXYzcUU)

```
gunzip -c all_genes.tpm.tsv.gz > mg/all_genes.tpm.tsv
rm all_genes.tpm.tsv.gz

gunzip -c all_genes.raw_counts.tsv.gz > mg/all_genes.raw_counts.tsv
rm all_genes.raw_counts.tsv.gz
```

**Metatranscriptomic abundances**
```
mkdir mt
```

* [LMO gene TPM values](https://)
* [LMO gene RAW counts](https://)

**TODO: Add info on how to obtain the RNA data**