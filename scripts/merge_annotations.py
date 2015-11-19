#!/usr/bin/env python
"""Merges protein families based on their localization on entries in a protein database"""

import sys, csv
from argparse import ArgumentParser

def readlimits(famfile):
    if not famfile: return []
    limits = []
    hin = open(famfile)
    for line in hin: limits.append(line.rstrip())
    hin.close()
    return limits

def readannots(infile, limits = []):
    a = {}
    hin = open(infile, 'r')
    hincsv = csv.reader(hin, delimiter = '\t')
    for i,row in enumerate(hincsv):
        if i==0: continue
        fams = []
        for item in row[1:]:
            item = item.rstrip(";")
            fams+=item.split(";")
        fams = list(set(fams))
        try: fams.remove("")
        except ValueError: pass

        for fam in fams:
            if len(limits)>0 and not fam in limits: continue
            for fam2 in fams:
                if len(limits)>0 and not fam2 in limits: continue            
                if fam==fam2: continue
                try: a[fam].append(fam2)
                except KeyError: a[fam] = [fam2]
    hin.close()

    for key,value in a.iteritems(): a[key] = list(set(value))
    return a

def mergeannots(a):
    am = {}
    parsed = {}
    i = 0
    for fam,famlist in a.iteritems():
        try: 
            parsed[fam]
            continue
        except KeyError: pass
        famlist.append(fam)
        print "Start:",len(famlist),
        curlist = [x for x in famlist]
        newfams = []
        while True:
            for f in curlist:
                f_tmp = a[f]
                f_diff = set(f_tmp).difference(set(famlist))
                newfams+=list(f_diff.difference(set(newfams)))
            famlist+=newfams
            if len(newfams)==0: break
            curlist = [x for x in newfams]
            newfams = []
        print "End:",len(famlist)
        am["T"+str(i)] = famlist
        i+=1
        for f in famlist:
            parsed[f] = ""
        
    return am

def write(am):
    hout = sys.stdout
    houtcsv = csv.writer(hout,delimiter = '\t')
    houtcsv.writerow(["TCluster","PFAM","TIGRFAM","COG"])

    for tc,fams in am.iteritems():
        cogs = []
        tigrs = []
        pfams = []
        for f in fams:
            try:
                if f[0:4]=="TIGR": tigrs.append(f)
                elif f[0:3] == "COG": cogs.append(f)
                elif f[0:2] == "PF": pfams.append(f)
            except IndexError: continue
        houtcsv.writerow([tc,"|".join(pfams),"|".join(tigrs),"|".join(cogs)])
    hout.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True,
            help="Infile with cross-reference between protein families (required)")
    parser.add_argument("-f", "--famfile", type=str, required=False,
            help="Limit merging to these families (optional)")
    args = parser.parse_args()

    if not args.infile: sys.exit(parser.print_help())

    limits = readlimits(args.famfile)

    a = readannots(args.infile, limits)
    am = mergeannots(a)
    write(am)
    
if __name__ == '__main__': 
    main()
