#!/usr/bin/env python

import networkx as nx, pandas as pd, sys, csv
from argparse import ArgumentParser

def rowsplit(s): return s.rstrip(";").split(";")

def tab_to_dataframe(infile,families):
    df = pd.DataFrame(columns=['Node1','Node2'])
    j = 0
    with open(infile, 'r') as fh:
        for i,row in enumerate(csv.reader(fh,delimiter="\t")):
            if i==0: continue
            [gene,pf,tigr,cog] = row
            pfams = rowsplit(pf)
            tigrfams = rowsplit(tigr)
            cogs = rowsplit(cog)
            store = pfams+tigrfams+cogs
            if families: store = list(set(store).intersection(set(families)))
            elif len(store)==1:
                tmp = pd.DataFrame(data=[store[0],""],index=["Node1","Node2"],columns=[j]).T
                df = pd.concat([df,tmp])
                j+=1
            for i,fam in enumerate(store):
                for fam2 in store:
                    if fam==fam2: continue
                    tmp = pd.DataFrame(data=[fam,fam2],index=["Node1","Node2"],columns=[j]).T
                    df = pd.concat([df,tmp])
                    j+=1
    remaining = list(set(families).difference(set(df.Node1)))
    for fam in remaining:
        tmp = pd.DataFrame(data=[fam,""],index=["Node1","Node2"],columns=[j]).T
        df = pd.concat([df,tmp])
        j+=1
    return df

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Infile cross-reference table with columns ['gene_id','PFAMs','TIGRFAMs','COGs']")
    parser.add_argument("-f", "--families", nargs="*",
            help="Only create links for these families")
    parser.add_argument("-o", "--outfile", 
            help="Write table to outfile. Defaults to stdout")

    args = parser.parse_args()
    
    linkdf = tab_to_dataframe(args.infile, args.families)
    
    if args.outfile: out = args.outfile
    else: out = sys.stdout
    linkdf.to_csv(out,sep="\t")

if __name__ == '__main__':
    main()
