#!/usr/bin/env python

import csv, sys
from argparse import ArgumentParser

def readtrans(transfile):
    hin = open(transfile)
    hincsv = csv.reader(hin, delimiter = '\t')
    t2fams = {}
    for row in hincsv:
        t = row[0]
        t2fams[t] = []
        for item in row[1:]:
            for f in item.split("|"): 
                if f == "": continue
                t2fams[t].append(f)
    return t2fams

def refine(r, t2fams):
    hin = open(r)
    hincsv = csv.reader(hin, delimiter = '\t')
    trans = {}
    for row in hincsv:
        t = row[0]
        trans[t] = []
        transporters = row[-1].split("|")
        for transporter in transporters: trans[t] += t2fams[transporter]
    return trans

def write(trans):
    hout = sys.stdout
    houtcsv = csv.writer(hout,delimiter = '\t')
    for t,fams in trans.iteritems():
        cogs = []
        tigrs = []
        pfams = []
        other = []
        for f in fams:
            try:
                if f[0:4]=="TIGR": tigrs.append(f)
                elif f[0:3] == "COG": cogs.append(f)
                elif f[0:2] == "PF": pfams.append(f)
                else: other.append(f)
            except IndexError: continue
        houtcsv.writerow([t,"|".join(pfams),"|".join(tigrs),"|".join(cogs),"|".join(other)])
    hout.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("-t", "--trans", required=True,
            help="Transport definitions (mapping transporters to families)")
    parser.add_argument("-r", "--refine", required=True,
            help="Refined transporter mergings")
    args = parser.parse_args()

    t2fams = readtrans(args.trans)

    trans = refine(args.refine,t2fams)

    write(trans)

if __name__ == '__main__':
    main()
