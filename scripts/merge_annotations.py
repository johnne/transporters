#!/usr/bin/env python
"""Merges protein families based on their localization on entries in a protein database"""

import sys, csv
from argparse import ArgumentParser

def readlimits(famfile):
    '''Reads protein families to include (to the exclusion of all others)'''
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
        
        ## Handle case if there is only one family for the row
        if len(fams)==1:
            if len(limits)>0 and not fams[0] in limits: continue ## Don't add if we're filtering
            a[fams[0]] = []
            continue

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
    am = {} ## Holds the merged annotations, keys are transporters, values are lists of families
    parsed = {} ## Holds all families that have been parsed already as keys
    filtered = [] ## Holds families with too high edgecount
    i = 0
    ## Iterate protein family and the list of protein families it matches to
    for fam,famlist in a.iteritems():
        ## If famlist is empty, add the family as its own transporter
        if len(famlist)==0:
            am["T"+str(i)] = [fam]
            parsed[fam] = ""
            i+=1
            continue

        ## If family has been parsed, continue
        if fam in parsed.keys(): continue
        ## If this is the first time parsing this family, append it to its own list
        famlist.append(fam)
        ## Make a copy of the family list
        curlist = [x for x in famlist]
        ## Initialize the list where we'll save all new families
        newfams = []
        while True:
            for f in curlist:
                ## Iterate the current list and get list of matching families
                fams_tmp = a[f]
                edges = len(fams_tmp)
                ## Check outgoing edges from family 'f',
                ## if greater than the edgecount threshold, add it to the 'filtered' list and continue
                if edges > edgecount:
                    filtered.append(f)
                    continue 
                ## If not, add family 'f' to fams_tmp
                fams_tmp.append(f)
                ## Then add the new families detected for 'f' to 'newfams'
                f_diff = set(fams_tmp).difference(set(famlist))
                newfams+=list(f_diff.difference(set(newfams)))
            ## Add all newly found families to 'famlist'
            famlist+=newfams
            ## If there are no new families, break the loop
            if len(newfams)==0: break
            ## Update the current list so that we will only loop the new families found
            curlist = [x for x in newfams]
            ## Reset newfams
            newfams = []
        ## Remove filtered families from famlist
        filtered_famlist = list(set(famlist).difference(set(filtered)))
        ## Assign transporters to the am dictionary
        if len(filtered_famlist)>0:
            am["T"+str(i)] = filtered_famlist
            i+=1
        ## Add all parsed families to parsed
        for f in famlist: parsed[f] = ""
    return (am,list(set(filtered)))

def write(am):
    hout = sys.stdout
    houtcsv = csv.writer(hout,delimiter = '\t')
    houtcsv.writerow(["TransportGroup","PFAM","TIGRFAM","COG"])
    for tg,fams in am.iteritems():
        cogs = []
        tigrs = []
        pfams = []
        for f in fams:
            try:
                if f[0:4]=="TIGR": tigrs.append(f)
                elif f[0:3] == "COG": cogs.append(f)
                elif f[0:2] == "PF": pfams.append(f)
            except IndexError: continue
        houtcsv.writerow([tg,"|".join(pfams),"|".join(tigrs),"|".join(cogs)])
    hout.close()

def write_filtered(filtered_fams, a):
    edgecounts = {}
    for f in filtered_fams: edgecounts[f] = len(a[f])+1
    ## Sort by edgecount
    import operator
    edgecounts = {}
    sorted_edges = sorted(edgecounts.items(), key=operator.itemgetter(1), reverse=True)
    hout = sys.stderr
    for item in sorted_edges: hout.write(item[0]+"\t"+str(item[1])+"\n")
    hout.close()


def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str, required=True,
            help="Infile with cross-reference between protein families (required)")
    parser.add_argument("-f", "--famfile", type=str, required=False,
            help="Limit merging to these families (optional)")
    parser.add_argument("-e", "--edgecount", type=int, default=6,
            help="Maximum number of outgoing connections from a single protein family. Defaults to 6. \
                    Choose a lower number to limit the size of merged transporter clusters.")
    args = parser.parse_args()

    if not args.infile: sys.exit(parser.print_help())

    limits = readlimits(args.famfile)
    a = readannots(args.infile, limits)
    (am,filtered_fams) = mergeannots(a, args.edgecount)
    write(am)
    
    ## Write protein families with more outgoing edges than the limit
    write_filtered(filtered_fams, a)

if __name__ == '__main__': 
    main()
