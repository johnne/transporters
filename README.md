# transporters
Analysis of transporters in genomes and metagenomes.

## Identifying protein families
    regexp="[Tt]ransport|[Ee]fflux|[Uu]ptake|[Ss]ymport|[Aa]ntiport|[Ii]mport|[Pp]ermease|[Pp]hosphotransferase|[Ee]xport|[Ss]olute[- ]binding|[Ss]ecretory|[Ss]ecretion|[Tt]on[Bb]|[Ss]u[Dd]"
    python scripts/print_hmm_db.py Pfam-A.27.0.hmm | egrep $regexp > data/Pfam-A.27.0_regexp_match.tab
    python scripts/print_hmm_db.py TIGRFAM.15.hmm | egrep $regexp > data/TIGRFAM.15_regexp_match.tab
    python scripts/print_cog_db.py | egrep $regexp > data/COG_regexp_match.tab
