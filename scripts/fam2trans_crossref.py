#!/usr/bin/env python

import csv,sys
from argparse import ArgumentParser

def readtrans(transfile):
    hin = open(transfile)
    hincsv = csv.reader(hin, delimiter = '\t')
    f2t = {}
    for row in hincsv:
        t = row[0]
        for item in row[1:]:
            for f in item.split("|"): 
                if f == "": continue
                f2t[f] = t
    return f2t

def write(f2t, f):
    hin = open(f)
    hincsv = csv.reader(hin, delimiter = '\t')
    houtcsv = csv.writer(sys.stdout, delimiter = '\t')
    i = 0
    for row in hincsv:
        l = []
        for item in row[1:]:
            for fam in item.split(";"):
                try: l.append(f2t[fam])
                except KeyError: continue
        if len(l)==0: continue
        i+=1
        newrow = [i]
        newrow.append(";".join(list(set(l))))
        houtcsv.writerow(newrow)
    hin.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("-t", "--trans", type=str, required=True,
            help="Transporter mergings")
    parser.add_argument("-f", "--fams", type=str, required=True,
            help="Family cross-ref")
    args = parser.parse_args()

    f2t = readtrans(args.trans)

    write(f2t, args.fams)

if __name__ == '__main__':
    main()
