# transporters
Analysis of transporters in genomes and metagenomes.

## Identifying protein families
Pfam v. 27.0 (ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/) and TIGRFAM v. 15.0 (ftp://ftp.jcvi.org/pub/data/TIGRFAMs/)
databases were downloaded and extracted. The COG database was accessed directly. The protein family descriptions were then matched using regular expression:


    regexp="[Tt]ransport|[Ee]fflux|[Uu]ptake|[Ss]ymport|[Aa]ntiport|[Ii]mport|[Pp]ermease|[Pp]hosphotransferase|[Ee]xport|[Ss]olute[- ]binding|[Ss]ecretory|[Ss]ecretion|[Tt]on[Bb]|[Ss]u[Dd]"
    python scripts/print_hmm_db.py Pfam-A.27.0.hmm | egrep $regexp > data/Pfam-A.27.0_regexp_match.tab
    python scripts/print_hmm_db.py TIGRFAM.15.hmm | egrep $regexp > data/TIGRFAM.15_regexp_match.tab
    python scripts/print_cog_db.py | egrep $regexp > data/COG_regexp_match.tab


## Merging protein families
A cross-reference table was downloaded from [Uniprot](http://www.uniprot.org/uniprot/?query=*&fil=reviewed%3Ayes) with TIGRFAM, PFAM and eggNOG annotations 
for 549,832 reviewed proteins and saved as [data/uniprot.2015_11.cross_ref.tab](data/uniprot.2015_11.cross_ref.tab). All entries matching the protein families
identified as transporters were then matched and stored:

    trans_fams=`cut -f1 data/COG_regexp_match.tab data/Pfam-A.27.0_regexp_match.tab data/TIGRFAM.15_regexp_match.tab | cut -f1 | tr '\n' '|' | sed 's/|$//g'`
    egrep "$trans_fams" data/uniprot.2015_11.cross_ref.tab > data/uniprot.2015_11.cross_ref.regexp_match.tab

Protein families were then merged based on entries in the Uniprot database:
    
    cogfile="data/COG_regexp_match.tab"
    pfamfile="data/Pfam-A.27.0_regexp_match.tab"
    tigrfile="data/TIGRFAM.15_regexp_match.tab"
    python scripts/merge_annotations.py -i data/uniprot.2015_11.cross_ref.regexp_match.tab -f <(cut -f1 $cogfile $pfamfile $tigrfile) > data/transporters.merged.tab 2> data/transporter.families.filtered.tab

Transport clusters, protein families and descriptions were then collated:

    python scripts/print_merged_to_multiline.py -i data/transporters.merged.tab -d <(cat $cogfile $pfamfile $tigrfile) > data/transporters.merged.multiline.tab

### Merging protein families by operon predictions
**TODO**
