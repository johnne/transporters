# transporters
Analysis of transporters in genomes and metagenomes.

## Identifying protein families
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


## Merging protein families
### Uniprot cross-reference
A cross-reference table was downloaded from [Uniprot](http://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes) with TIGRFAM, PFAM and eggNOG annotations 
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

### Operon predictions
Protein family merging was further refined using gene operon predictions downloaded from OperonDB (ftp://ftp.cbcb.umd.edu/pub/data/operondb/operon_predictions.tgz). All predictions were concatenated
and parsed using [gene_operons_to_fams.py](scripts/gene_operons_to_fams.py):

    python scripts/gene_operons_to_fams.py -g data/uniprot.2015_11.cross_ref.regexp_match.ids.to.gi.tab.txt -f data/uniprot.2015_11.cross_ref.regexp_match.tab -o data/all.operons > data/families.merged.operons.tab

The merging table from gene operons was combined with the table from Uniprot annotations:

    python scripts/merge_annotations.py -i <(cat data/uniprot.2015_11.cross_ref.regexp_match.tab data/families.merged.operons.tab) -f <(cut -f1 $cogfile $pfamfile $tigrfile) > data/transporters.operons.merged.tab 2> data/transporter.operons.families.filtered
