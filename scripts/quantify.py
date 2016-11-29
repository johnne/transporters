#!/usr/bin/env python

import pandas as pd, numpy as np, logging, os, sys
from argparse import ArgumentParser
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def stats(tcov,cov,f):
    transporters_found = len(set(tcov.transporter))
    transporter_fractions = tcov.ix[:,2:].sum().div(cov.sum())*100
    tminsum = transporter_fractions.min()
    tmaxsum = transporter_fractions.max()
    tmeansum = transporter_fractions.mean()
    fn = os.path.basename(f)
    sys.stderr.write(fn+" "+str(transporters_found)+" transporters " + str(np.round(tminsum,2))+"-"+str(np.round(tmaxsum,2))+" mean:"+str(np.round(tmeansum,2))+"\n")

def get_rep(fams,dfsum):
    return list(dfsum.loc[fams].mean(axis=1).sort_values(ascending=False).index)[0]

def add_def(tmean,unclass,dfsum,famdf,cdef,i):
    for t in unclass:
        fams = list(set(cdef.loc[cdef.transporter==t,"family"]).intersection(set(dfsum.index)))
        if len(fams)==1: f = fams[0]
        else: f = get_rep(fams,dfsum)
        d = list(famdf.loc[famdf.family==fams[0],1].values)[0]
        if i > 0: tmean.loc[t,:i] = [d]*i
        else: tmean = pd.concat([pd.DataFrame(columns=["Name"],index=[t],data={"Name":d}),tmean],axis=1)
    return tmean

def write_reps(tmean, dfsum, cdef, f):
    with open(f, 'w') as fh:
        for t in tmean.index:
            fams = list(set(cdef.loc[cdef.transporter==t,"family"]).intersection(set(dfsum.index)))
            r = get_rep(fams,dfsum)
            fh.write(t+"\t"+r+"\n")

def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--definitions", required=True,
            help="Transport cluster definitions. Tab delimited with columns ['transport_id','family']")
    parser.add_argument("-a", "--annotations", required=True,
            help="Annotation file. Tab delimited with columns ['orf_id','family_id']")
    parser.add_argument("-q", "--quant", required=True,
            help="ORF abundance file. Rows are ORFs and columns are samples")
    parser.add_argument("-o", "--outfile", type=str,
            help="Write transporter means to outfile. Defaults to stdout")
    parser.add_argument("-c", "--classif", type=str,
            help="Classification file for transporters (optional)")
    parser.add_argument("-r", "--reps", type=str,
            help="Write family with highest mean across samples for each cluster to file")

    args = parser.parse_args()
    ## Definitions file
    logging.info("Reading definitions")
    cdef = pd.read_csv(args.definitions, header=None, sep="\t", names=["transporter","family"])
    logging.info("Read definitions for "+str(len(set(cdef.transporter)))+" transporters")
    ## Annotations file
    logging.info("Reading annotations for ORFs")
    ann = pd.read_csv(args.annotations, header=None, sep="\t", names=["orf","family"],usecols=[0,1])
    ann.family = ann.family.str.replace("PFAM","PF")
    logging.info("Read annotations for "+str(len(ann))+" ORFs")
    ## Coverage file
    logging.info("Reading coverage for ORFs")
    cov = pd.read_csv(args.quant, header=0, sep="\t", index_col=0)
    ## Merge annotations and coverage
    logging.info("Calculating sum for protein families")
    df = pd.merge(ann,cov,left_on="orf",right_index=True)
    ## Sum to protein family
    dfsum = df.groupby("family").sum()
    ## Merge protein family sum with transporter definitions
    logging.info("Merging with transporter definitions")
    tcov = pd.merge(cdef,dfsum,left_on="family",right_index=True)
    ## Calculate transporter means
    logging.info("Calculating transporter means")
    tmean = tcov.groupby("transporter").mean()
    ##Sort by cluster size
    order = cdef.loc[cdef.transporter.isin(tmean.index)].groupby("transporter").count().sort_values("family",ascending=False).index
    ## Classifications 
    if args.classif:
        logging.info("Adding classifications from "+args.classif)
        cclass = pd.read_csv(args.classif, header=0, sep="\t", index_col=0)
        tmean = pd.merge(cclass,tmean,left_index=True,right_index=True)
    
    if args.outfile: tmean.to_csv(args.outfile, sep="\t")
    else: tmean.to_csv(sys.stdout, sep="\t")
    ## Write representatives
    if args.reps: write_reps(tmean,dfsum,cdef,args.reps)
    ## Write stats
    stats(tcov,cov,args.annotations)

if __name__ == '__main__':
    main()
