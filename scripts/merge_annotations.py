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

def mergeannots(a, edgecount):
    am = {}
    parsed = {}
    filtered = []
    i = 0
    for fam,famlist in a.iteritems():
        try: 
            parsed[fam]
            continue
        except KeyError: pass
        famlist.append(fam)
        curlist = [x for x in famlist]
        newfams = []
        while True:
            for f in curlist:
                f_tmp = a[f]
                edges = len(f_tmp)
                ## Set limit on outgoing edges
                if edges > edgecount:
                    filtered.append(f)
                    continue 
                f_diff = set(f_tmp).difference(set(famlist))
                newfams+=list(f_diff.difference(set(newfams)))
            famlist+=newfams
            if len(newfams)==0: break
            curlist = [x for x in newfams]
            newfams = []
        filtered_famlist = list(set(famlist).difference(set(filtered)))
        if len(filtered_famlist)>0:
            am["T"+str(i)] = filtered_famlist
            i+=1
        for f in famlist:
            parsed[f] = ""
    return (am,list(set(filtered)))

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
    parser.add_argument("-e", "--edgecount", type=int, default=3,
            help="Maximum number of outgoing connections from a single protein family. Defaults to 3. \
                    Choose a lower number to limit the size of merged transporter clusters.")
    args = parser.parse_args()

    if not args.infile: sys.exit(parser.print_help())

    limits = readlimits(args.famfile)
    a = readannots(args.infile, limits)
    (am,filtered_fams) = mergeannots(a, args.edgecount)
    write(am)
    
    ## Write protein families with more outgoing edges than the limit
    hout = sys.stderr
    for f in filtered_fams: hout.write(f+"\n")
    hout.close()

if __name__ == '__main__': 
    main()
