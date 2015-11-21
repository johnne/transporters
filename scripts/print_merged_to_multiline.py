#!/usr/bin/env python

from argparse import ArgumentParser
import csv,sys

def readdescs(descriptions):
    d = {}
    hin = open(descriptions)
    hincsv = csv.reader(hin, delimiter = '\t')
    for row in hincsv:
        d[row[0]] = row[1]
    hin.close()
    return d

def readmerged(infile):
    tgroups = {}
    hin = open(infile)
    hincsv = csv.reader(hin, delimiter = '\t')
    for i,row in enumerate(hincsv):
        if i==0:continue
        fams = []
        for item in row[1:]:
            if item=="": continue
            fams+=item.split("|")
        tgroups[row[0]] = fams
    hin.close()
    return tgroups

def write(tgroups, descs):
    hout = sys.stdout
    houtcsv = csv.writer(hout, delimiter = '\t')
    houtcsv.writerow(["TransporterGroup","Family","Description"])
    for t,fams in tgroups.iteritems():
        for fam in fams:
            houtcsv.writerow([t,fam,descs[fam]])
    hout.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True,
            help="Infile with annotations merged per transporter cluster")
    parser.add_argument("-d", "--descriptions", type=str, required=True,
            help="File with descriptions for protein families")

    args = parser.parse_args()

    descs = readdescs(args.descriptions)

    tgroups = readmerged(args.infile)

    write(tgroups, descs)

if __name__ == '__main__':
    main()
