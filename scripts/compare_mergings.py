#!/usr/bin/env python

from argparse import ArgumentParser
import csv, sys

def read_merged(f):
    m = {}
    t = {}
    hin = open(f)
    hincsv = csv.reader(hin, delimiter = '\t')
    for i, row in enumerate(hincsv):
        if i==0:continue
        m[row[1]] = {"T": row[0], "D": row[2]}
        try: t[row[0]].append(row[1])
        except KeyError: t[row[0]] = [row[1]]
    hin.close()
    return m,t

def compare(m1,t1,m2,t2):
    hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '\t')

    filtered = len(set(m1.keys()).difference(set(m2.keys())))

    for f, d in m1.iteritems():
        out = [f,d["D"],"","",""]
        if not f in m2.keys():
            out[2] = "FILTERED"
            houtcsv.writerow(out)
            continue

        ## If family is present in both, check overlap and mergings
        t_1 = d["T"]
        mergefams1 = t1[t_1] ## List of families merged to current family in m1
        
        t_2 = m2[f]["T"]
        mergefams2 = t2[t_2] ## Families merged to current family in m2
        
        ## Get lengths
        l1 = len(set(mergefams1))
        l2 = len(set(mergefams2))

        ## Calculate overlap and difference
        overlap = len(set(mergefams1).intersection(set(mergefams2)))
        if overlap == 0: overlap = "NONE|NONE"
        else: overlap = str(overlap)+"/"+str(l1)+"|"+str(overlap)+"/"+str(l2)
        out[3] = overlap

        diff1 = set(mergefams1).difference(set(mergefams2))
        diff2 = set(mergefams2).difference(set(mergefams1))

        if len(diff1)==len(diff2)==0: out[4] = "IDENTICAL"
        else: out[4] = "|".join(list(diff2))

        if l2>l1: 
            out[2] = "MERGED"
            houtcsv.writerow(out)
            continue
        elif len(mergefams2)==len(mergefams1):
            out[2] = "UNTOUCHED"
            houtcsv.writerow(out)
            continue
        else:
            out[2] = "IN SMALLER"
            houtcsv.writerow(out)
            continue
    hout.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("--m1", type=str, required=True,
            help="First merging file")
    parser.add_argument("--m2", type=str, required=True,
            help="Second merging file")
    args = parser.parse_args()

    m1,t1 = read_merged(args.m1)
    m2,t2 = read_merged(args.m2)

    compare(m1,t1,m2,t2)

if __name__ == '__main__': 
    main()
