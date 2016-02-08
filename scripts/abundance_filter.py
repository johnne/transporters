#!/usr/bin/env python

from argparse import ArgumentParser
import csv, sys, pandas as pd

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True,
            help="Transporter (rows) coverage in samples (columns)")
    parser.add_argument("-f", "--filter", type=float, default=1,
            help="Lowest maximum abundance across samples to include transporter")
    parser.add_argument("-s", "--sort", action="store_true",
            help="Also sort transporters by sum across all samples")

    args = parser.parse_args()

    tcov = pd.read_csv(args.infile, header=0, sep="\t", index_col=0)
    ## Normalize to percentages
    tcovn = tcov.div(tcov.sum())*100
    
    t_in = tcovn[tcovn.max(axis=1)>=args.filter].index
    t_out =  tcovn[tcovn.max(axis=1)<args.filter].index

    tcovn_in = tcovn.loc[t_in]

    if args.sort: tcovn_in = tcovn_in.loc[tcovn_in.sum(axis=1).sort(ascending=False, inplace=False).index]

    tcovn_in.to_csv(sys.stdout, sep="\t")    

if __name__ == '__main__':
    main()
    
