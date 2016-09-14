#!/usr/bin/env python

import pandas as pd, gzip as gz, sys, csv, logging
from argparse import ArgumentParser

def ParserOrder(orderstring):
    orderitems = orderstring.split(",")
    ordernum = {}
    for i,item in enumerate(orderitems):
        try: 
            if item[0:4] == "TIGR": ordernum[i] = "TIGRFAM"
        except IndexError: continue
        try:
            if item[0:2] == "PF": ordernum[i] = "PFAM"
        except IndexError: continue
        try:
            if item[0:3] == "COG": ordernum[i] = "COG"
        except IndexError: continue
    return ordernum

def CalculateCoverage(transporters,orfcov,orfann,ordernum,method):
    reps = {}
    tcov = {}
    for t in transporters.index:
        fams = False
        for i in range(0,len(ordernum)):
            try: 
                fams = transporters.loc[t,ordernum[i]].split("|")
                break
            except AttributeError:  continue
        if not fams: continue
    
        ##Sum coverage per family defined for transporter
        t_orfann = orfann[orfann[1].isin(fams)] ## Get ORFs matching families
        if len(t_orfann) == 0: continue ## If no ORFs, skip transporter
    
        ## If error thrown here, no ORFs for this transporter group remain to calculate coverage
        try: t_orfann_cov = pd.concat([t_orfann,orfcov.loc[t_orfann.index]],axis=1)
        except KeyError: continue 
        
        t_fam_cov = t_orfann_cov.groupby(1).sum() ## Groupby family and sum
        if len(t_fam_cov)==1: rep = t_fam_cov.index[0] ## If only 1 family, use it as representative
        else: 
            ## Otherwise, use method and sort
            if method == 'mean': rep = t_fam_cov.mean(axis=1).sort_values(ascending=False).index[0]
            elif method == 'max': rep = t_fam_cov.max(axis=1).sort_values(ascending=False).index[0]
        reps[t] = rep
        tcov[t] = t_fam_cov.loc[rep] ## Store coverage
    
        ## Remove ORFs used to calculate coverage here, to avoid using ORFs multiple times
        used_orfs=list(orfann[orfann[1].isin(fams)].index)
        orfcov.drop(used_orfs,inplace=True,errors="ignore")
    tcov = pd.DataFrame(tcov).T
    return tcov,reps

def AddDescriptions(tcov,descs,reps):
    '''Add description column to transport coverage table'''
    cols = ["Description"]+list(tcov.columns)
    reps_df = pd.DataFrame([reps], index=["Family"]).T
    transporter_desc = pd.merge(reps_df,descs,left_on="Family", right_index=True)
    tcov = pd.merge(tcov,transporter_desc,left_index=True,right_index=True)[cols]
    return tcov    

def main():
    parser = ArgumentParser()
    parser.add_argument("-t", "--transporters", type=str, required=True,
            help="Transporter definitions, tab separated. Format: <Transporter> <PFAMs> <TIGRFAMs> <COGs> <Other db>")
    parser.add_argument("-c", "--cov", type=str, required=True,
            help="Normalized coverage for ORFs in all samples")
    parser.add_argument("-a", "--annotations", required=True,
            help="Annotations for ORFs, one annotation per line")
    parser.add_argument("-o", "--outfile", 
            help="Write transporter coverage to outfile")
    parser.add_argument("--descriptions", type=str,
            help="(Optional) Specify tab-delimited file with descriptions for families. Will be applied to transporters")
    parser.add_argument("--order", type=str, default="TIGR,COG,PF",
            help="Priority order for protein family databases. Defaults to TIGRFAM,COG,PFAM")
    parser.add_argument("--method", type=str, default="mean",
            help="Specify method (mean,max to use for selecting representative protein family for transporter groups.\
            Defaults to 'mean' which uses the family with the highest mean across samples.")
    parser.add_argument("-v", "--verbose", action="store_true",
            help="increase output verbosity")

    args = parser.parse_args()
    
    if args.verbose: logging.basicConfig(level=logging.INFO)
    
    ordernum = ParserOrder(args.order)
    
    orfcov = pd.read_csv(args.cov, index_col=0, sep="\t", header=0)
    orfcov = orfcov.div(orfcov.sum())*100 ## Normalize
    logging.info("Read coverage for "+str(len(orfcov))+" ORFs in "+str(len(orfcov.columns))+" samples.")

    orfann = pd.read_csv(args.annotations, index_col=0, sep="\t", header=None)
    logging.info("Read "+str(len(set(orfann[1])))+" annotations for "+str(len(set(orfann.index)))+" ORFs.")      

    transporters = pd.read_csv(args.transporters, index_col=0, sep="\t", header=None, names=["PFAM","TIGRFAM","COG","OTHER"])
    transporters.sort_values(by=ordernum[0], inplace=True)
    logging.info("Read definitions for "+str(len(transporters))+" transporters.")
    
    logging.info("Calculating coverage for transporters...")
    tcov,reps = CalculateCoverage(transporters,orfcov,orfann,ordernum,args.method)
    logging.info(str(len(tcov)) + " transporters remaining after calculating coverage")

    if args.descriptions:
        logging.info("Adding descriptions...")
        try: 
            descs = pd.read_csv(args.descriptions, header=None, sep="\t", index_col=0, names=["Description"])
            tcov = AddDescriptions(tcov,descs,reps)
        except IOError: pass  

    ## Write to file
    if args.outfile:
        tcov.to_csv(args.outfile, sep = '\t')
        logging.info("Wrote to "+args.outfile+".")
    else:
        tcov.to_csv(sys.stdout, sep = '\t')
        logging.info("Wrote to stdout.")

if __name__ == '__main__':
    main()
