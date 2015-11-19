#!/usr/bin/env python

from argparse import ArgumentParser
import csv, sys

def read_gi2uni(f):
    m = {}
    hin = open(f)
    hincsv = csv.reader(hin, delimiter = '\t')
    for i,row in enumerate(hincsv):
        if i==0: continue
        m[row[1]] = row[0]
    hin.close()
    return m

def read_uni2fam(f):
    m = {}
    hin = open(f)
    hincsv = csv.reader(hin, delimiter = '\t')
    for row in hincsv:
        uni = row[0]
        m[uni] = []
        for item in row:
            for fam in (item.rstrip(";")).split(";"):
                if fam[0:4] == "TIGR" or fam[0:3] == "COG" or fam[0:2] == "PF": m[uni].append(fam)
    hin.close()
    return m            

def read_tclust2fam(f):
    m = {}
    hin = open(f)
    hincsv = csv.reader(hin, delimiter = '\t')
    for i, row in enumerate(hincsv):
        if i==0: continue
        m[row[1]] = row[0]
    hin.close()
    return m

def match(gi, gi2uni, uni2fams):
    try: uni = gi2uni[gi]
    except KeyError: return False

    fams = uni2fams[uni]
    fams = list(set(fams))
    return fams

def read_operons(f, gi2uni, uni2fams, conf):
    hin = open(f)
    hincsv = csv.reader(hin, delimiter = ' ')
    for row in hincsv:
        gi1 = row[0]
        gi2 = row[1]
        c = float(row[2])
        if c<conf: continue ## Filter operons by confidence
        try:
            gi2uni[gi1]
            gi2uni[gi2]
        except KeyError: continue

        fams1 = match(gi1,gi2uni,uni2fams)
        fams2 = match(gi2,gi2uni,uni2fams)
        if fams1 and fams2:
            print "|".join(fams1),"|".join(fams2)

def main():
    parser = ArgumentParser()
    parser.add_argument("-g", "--uniprottogi", type=str, required=True,
            help="Uniprot mapping file of UniprotKB to GI accessions")
    parser.add_argument("-f", "--uniprottofams", type=str, required=True,
            help="Uniprot to protein family annotations, cross-reference table")
    parser.add_argument("-o", "--operons", type=str, required=True,
            help="Operon database output file (see http://operondb.cbcb.umd.edu/cgi-bin/operondb/operons.cgi)")
    parser.add_argument("-c", "--confidence", type=float, default=50.0,
            help="Minimum confidence from operon predictions")

    args = parser.parse_args()

    gi2uni = read_gi2uni(args.uniprottogi)
    uni2fams = read_uni2fam(args.uniprottofams)
    read_operons(args.operons, gi2uni, uni2fams, args.confidence)

if __name__=='__main__':
    main()
