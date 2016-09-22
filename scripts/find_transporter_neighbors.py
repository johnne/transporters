#!/usr/bin/env python

import pandas as pd, sys, csv, logging
from argparse import ArgumentParser

def Filter(bed,orf2trans):
    '''Filters contigs'''
    ## Remove contigs without transporters
    contigs = list(set(bed[bed.ORF.isin(orf2trans.index)].index))
    bed_f = bed.loc[contigs]
    ## Remove contigs with only 1 ORF
    bed_f = bed_f.loc[bed_f.groupby(level=0).count().ix[:,0]>1]
    return bed_f

def Merge(bed_f, orf2trans, orfann):
    ## Merge dataframes, first bed with all ORF annotations
    x = pd.merge(bed_f,orfann,right_index=True,left_on="ORF")
    ## Next, merge with ORF-> transporter map
    y = pd.merge(x,orf2trans, left_on="ORF",right_index=True, how="outer")
    return y

def FindNeighbors(merged,distance,logging):
    ## Add index
    merged["Index"] = range(0,len(merged))
    merged["Contig"] = merged.index
    merged.index = merged.Index

    ## Get index for transporter ORFs
    transporter_index = merged[merged[1]==merged[1]].index
    
    t_neighbors = {}

    ## For progress bar
    part = len(transporter_index)/10
    progress = 0

    distance = 3
    for j,index in enumerate(transporter_index, start=1):
        if j%part == 0: 
            progress+=10
            logging.info(str(progress)+"%")
        contig = merged.loc[index,"Contig"]
        t = merged.loc[index,1] ## Current transporter
        try: t_neighbors[t]
        except KeyError: t_neighbors[t] = {}
        ## Generate index range
        r = range(index-distance,index+distance+1)
        ## Line index range up with matching contigs
        m = pd.DataFrame(merged.loc[r,"Contig"]==contig)
        r_m = m[m["Contig"]==True].index ## Matching range
        ## Iterate the index range
        for i in r_m:
            if i==index: continue
            t_n = merged.loc[i,1] ## Transporter annotation
            if t_n==t_n: 
                try: t_neighbors[t][t_n]+=1
                except KeyError: t_neighbors[t][t_n] = 1
            a_n = merged.loc[i,"Family"]
            if type(a_n)==float: continue
            for item in a_n:
                try: t_neighbors[t][item]+=1
                except KeyError: t_neighbors[t][item] = 1
    ## Create dataframe
    df = pd.DataFrame(t_neighbors).T
    df.fillna(0,inplace=True)
    return df

def main():
    parser = ArgumentParser('''Finds neighbors to transporters on contigs and reports the annotated function''')
    parser.add_argument("-b", "--bed", required = True,
            help="Bed file for assembly. Defines the contigs and predicted ORFs. First 4 columns need to be <Contig> <Start> <End> <ORF id>.")
    parser.add_argument("-t", "--orf2trans", required=True,
            help="Table with ORF-> transporter mapping (see --orftable flag in calculate_transporter_abundance.py). \
                    Format is: <ORF> <transporter> [protein family].")
#    parser.add_argument("-t", "--transporters", type=str, required=True,
#            help="Transporter definitions, tab separated. Format: <Transporter> <PFAMs> <TIGRFAMs> <COGs> <Other db>")
    parser.add_argument("-a", "--annotations", required=True,
            help="Annotations for ORFs, one annotation per line")
    parser.add_argument("-d", "--distance", default=3,
            help="Distance in ORFs from transporter to analyze. Defaults to 3.")
    parser.add_argument("-o", "--outfile", type=str,
            help="Write data frame with protein family/transporter neighboring counts for each transporter")
    parser.add_argument("-v", "--verbose", action="store_true",
            help="increase output verbosity")

    args = parser.parse_args()    

    if args.verbose: logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.INFO)

    ## Read ORF->transporter map
    orf2trans = pd.read_csv(args.orf2trans, header=None, index_col=0, sep="\t")
    logging.info("Read "+str(len(orf2trans))+" ORF->transporter mappings")

    ## Read bed file
    bed = pd.read_csv(args.bed, sep="\t", header=None, index_col=0, names=["Contig","Start","End","ORF"], usecols=[0,1,2,3])
    logging.info("Read "+str(len(bed))+" ORF definitions on "+str(len(set(bed.index)))+" contigs")

    ## Filter bed file
    logging.info("Filtering to remove contigs without transporters/with too few ORFs")
    bed_f = Filter(bed,orf2trans)
    logging.info(str(len(set(bed_f.index)))+" contigs remaining after filtering")
    if len(bed_f.index)==0: sys.exit()

    ## Read ORF annotations
    orfann = pd.read_csv(args.annotations, index_col=0, sep="\t", header=None, names=["ORF","Family"], usecols=[0,1])
    logging.info("Read annotations for ORFs")

    ## Aggregate annotations into a list per ORF
    logging.info("Aggregating annotations...")
    orfann = orfann.groupby(level=0).agg(lambda col: list(col))
        
    ## Merge tables
    merged = Merge(bed_f, orf2trans, orfann)
    
    ## Find neighbors
    logging.info("Finding neighbors...")
    df = FindNeighbors(merged,args.distance,logging)
    
    ## Write
    if args.outfile: df.to_csv(args.outfile, sep="\t")
    else: df.to_csv(sys.stdout, sep="\t")

if __name__ == '__main__':
    main()
