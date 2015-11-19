# transporters
Analysis of transporters in genomes and metagenomes.

## Identifying protein families
Pfam [v. 27.0](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz) and Tigrfam [v. 15](ftp://ftp.jcvi.org/pub/data/TIGRFAMs/) 
databases were downloaded and extracted. The COG database was accessed directly. The protein family descriptions were then matched using regular expression:

    regexp="[Tt]ransport|[Ee]fflux|[Uu]ptake|[Ss]ymport|[Aa]ntiport|[Ii]mport|[Pp]ermease|[Pp]hosphotransferase|[Ee]xport|[Ss]olute[- ]binding|[Ss]ecretory|[Ss]ecretion|[Tt]on[Bb]|[Ss]u[Dd]"
    python scripts/print_hmm_db.py Pfam-A.27.0.hmm | egrep $regexp > data/Pfam-A.27.0_regexp_match.tab
    python scripts/print_hmm_db.py TIGRFAM.15.hmm | egrep $regexp > data/TIGRFAM.15_regexp_match.tab
    python scripts/print_cog_db.py | egrep $regexp > data/COG_regexp_match.tab

## Merging protein families
