# transporters
This repository contains scripts and data files for the analysis of transporters in genomes and metagenomes.

## Defining transporters
Here, the Pfam, Tigrfam and COG databases were screened for protein families related to transport in cells. 
The first step is data collection and identification using regular expressions. The second step merges protein families based on the annotation of
all reviewed entries in [Uniprot](http://www.uniprot.org). That is, if several protein families are annotated on the same protein they are merged into a group.
Some very broad families that linked to a lot of other families (e.g. the [PF00005](http://pfam.xfam.org/family/PF00005) ABC-transporter family) 
were ignored (see the --edgecount flag for [merge_annotations.py](scripts/merge_annotations.py)) to avoid too large groups.

This first merging of families into "Transport groups" was refined using operon predictions from [OperonDB](http://operondb.cbcb.umd.edu/cgi-bin/operondb/operons.cgi).
If transporter families were found to occupy the same predicted operon, the were further merged.

The workflow for these steps are detailed below:

### 1. Identifying protein families
Pfam v. 28.0 (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam28.0/) and TIGRFAM v. 15.0 (ftp://ftp.jcvi.org/pub/data/TIGRFAMs/)
databases were downloaded and extracted. 

    wget -O data/Pfam-A.28.0.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam28.0/Pfam-A.hmm.gz
    wget -O data/TIGRFAM.15 ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz
    tar xvf data/TIGRFAM.15.hmm.tar.gz
    cat *.HMM | gzip -c > data/TIGRFAM.15.hmm.gz
    rm *.HMM

The COG database was accessed directly. The protein family descriptions were then matched using regular expression:

    regexp="[Tt]ransport|[Ee]fflux|[Uu]ptake|[Ss]ymport|[Aa]ntiport|[Ii]mport|[Pp]ermease|[Pp]hosphotransferase|[Ee]xport|[Ss]olute[- ]binding|[Ss]ecretory|[Ss]ecretion|[Tt]on[Bb]|[Ss]u[Dd]"
    python scripts/print_hmm_db.py data/Pfam-A.28.0.hmm.gz | egrep "$regexp" > data/Pfam-A.28.0_regexp_match.tab
    python scripts/print_hmm_db.py data/TIGRFAM.15.hmm.gz | egrep "$regexp" > data/TIGRFAM.15_regexp_match.tab
    python scripts/print_cog_db.py | egrep "$regexp" > data/COG_regexp_match.tab


### 2. Merging protein families
#### 2.1 Uniprot cross-reference
A [cross-reference table](http://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes) was downloaded from [Uniprot](http://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes) with TIGRFAM, PFAM and eggNOG annotations 
for 549,832 reviewed proteins and saved as [data/uniprot.2015_11.cross_ref.tab](data/uniprot.2015_11.cross_ref.tab). All entries matching the protein families
identified as transporters were then matched and stored:

    trans_fams=`cut -f1 data/COG_regexp_match.tab data/Pfam-A.28.0_regexp_match.tab data/TIGRFAM.15_regexp_match.tab | cut -f1 | tr '\n' '|' | sed 's/|$//g'`
    egrep "$trans_fams" data/uniprot.2015_11.cross_ref.tab > data/uniprot.2015_11.cross_ref.regexp_match.tab

Protein families were then merged based on entries in the Uniprot database. By default, outgoing edges from a single protein family are limited to 3.
Families with more edges than this threshold are saved to [data/transporter.families.filtered](data/transporter.families.filtered):
    
    cogfile="data/COG_regexp_match.tab"
    pfamfile="data/Pfam-A.28.0_regexp_match.tab"
    tigrfile="data/TIGRFAM.15_regexp_match.tab"
    python scripts/merge_annotations.py -i data/uniprot.2015_11.cross_ref.regexp_match.tab -f <(cut -f1 $cogfile $pfamfile $tigrfile) > data/transporters.merged.tab 2> data/transporter.families.filtered

Transport clusters, protein families and descriptions were then collated:

    python scripts/print_merged_to_multiline.py -i data/transporters.merged.tab -d <(cat $cogfile $pfamfile $tigrfile) > data/transporters.merged.multiline.tab

#### 2.2 Operon predictions
Protein family merging was further refined using gene operon predictions downloaded from OperonDB (ftp://ftp.cbcb.umd.edu/pub/data/operondb/).

    wget -O data/operon_predictions.tgz ftp://ftp.cbcb.umd.edu/pub/data/operondb/operon_predictions.tgz
    cat prediction/*.operons | gzip -c > data/operons.gz
    rm -r prediction data/operon_predictions.tgz

Predictions were parsed using [gene_operons_to_fams.py](scripts/gene_operons_to_fams.py). 
For this, **a mapping table of UniProt accessions to NCBI GI numbers are needed**. Create a file containing the ids to map:

    cut -f1 data/uniprot.2015_11.cross_ref.regexp_match.tab > data/uniprot.2015_11.cross_ref.regexp_match.ids

Then go to http://www.uniprot.org/uploadlists/, upload the [data/uniprot.2015_11.cross_ref.regexp_match.ids](data/uniprot.2015_11.cross_ref.regexp_match.ids) list and choose to map From 'UniProktKB AC/ID' To 'GI number'. Download the mapping table and store
it as [data/uniprot.2015_11.cross_ref.regexp_match.ids.to.gi.tab](data/uniprot.2015_11.cross_ref.regexp_match.ids.to.gi.tab).

Next, parse the operon predictions.

    python scripts/gene_operons_to_fams.py -g data/uniprot.2015_11.cross_ref.regexp_match.ids.to.gi.tab.txt -f data/uniprot.2015_11.cross_ref.regexp_match.tab -o data/operons.gz > data/families.merged.operons.tab

The merging table from gene operons was combined with the table from Uniprot annotations:

    python scripts/merge_annotations.py -i <(cat data/uniprot.2015_11.cross_ref.regexp_match.tab data/families.merged.operons.tab) -f <(cut -f1 $cogfile $pfamfile $tigrfile) > data/transporters.operons.merged.tab 2> data/transporter.operons.families.filtered
