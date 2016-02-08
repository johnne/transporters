#!/usr/bin/env python

import sys,csv, pandas as pd
from argparse import ArgumentParser

def CalcTransCov(trans2fams, orf2ann, orfcov):
    orf2trans = {}
    r = {}
    for t in list(set(trans2fams.index)):
        families = trans2fams.loc[t,"Family"]
        if type(families)==str: families = [families]
        ## Get all ORFs corresponding to these families
        t_orfs = set(orf2ann[orf2ann.isin(families)].index)
        ## Assign transporter to orfs
        for orf in t_orfs: orf2trans[orf] = t
        ## Get number of families hit in the samples
        t_size = len(set(orf2ann[orf2ann.isin(families)]))
        ## Calculate coverage
        try: tcov = orfcov.loc[t_orfs].sum()
        except KeyError: 
            tcov = 0
            continue
        ## If coverage is 0, skip it
        if tcov.sum()==0: continue
        ## Divide by number of families hit
        r[t] = tcov/float(t_size)
    return r,orf2trans

def main():
    parser = ArgumentParser()
    parser.add_argument("-t", "--transporters", type=str, required=True,
            help="Transporter definitions")
    parser.add_argument("-c", "--genecov", type=str, required=True,
            help="Tab-delimited file with coverage for genes (rows) in samples (columns)")
    parser.add_argument("-a", "--annotations", type=str, required=True,
            help="Tab-delimited file with annotations for genes. Genes with multiple annotations (e.g. COG, PFAM, TIGRFAM) are given on multiple lines")

    args = parser.parse_args()

    orfann = pd.read_csv(args.annotations, header=0, sep="\t", index_col=0)
    orf2ann = pd.concat([orfann.COG,orfann.PFAM,orfann.TIGRFAM])
    orf2cov = pd.read_csv(args.genecov, header=0, sep="\t", index_col=0)
    trans2fams = pd.read_csv(args.transporters, header=0, sep="\t", index_col=0)

    transcov,orf2trans = CalcTransCov(trans2fams, orf2ann, orf2cov)
    tcov = pd.DataFrame(transcov).T
    tcov.fillna(0, inplace=True)
    tcov.to_csv(sys.stdout, sep='\t')
    
    o2t = pd.DataFrame([orf2trans]).T
    o2t.columns=["Transporter"]
    o2t.to_csv(sys.stderr, sep='\t')

if __name__ == '__main__': 
    main()
