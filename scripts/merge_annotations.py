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

def readcorr(corrfile):
    import pandas as pd
    corr =pd.read_csv(corrfile, sep="\t", index_col=0, header=0)
    return corr.to_dict()

def readannots(infile, frac, limits = [], corr = {}, corrmin=0):
    a = {}
    famcount = {}
    hin = open(infile, 'r')
    hincsv = csv.reader(hin, delimiter = '\t')
    for i,row in enumerate(hincsv):
        fams = []
        for item in row[1:]:
            item = item.rstrip(";")
            fams+=item.split(";")
        fams = list(set(fams))
        try: fams.remove("")
        except ValueError: pass

        for fam in fams:
            try: famcount[fam]+=1
            except KeyError: famcount[fam] = 1

        ## Handle case if there is only one family for the row
        if len(fams)==1:
            if len(limits)>0 and not fams[0] in limits: continue ## Don't add if we're filtering
            if not fams[0] in a.keys(): a[fams[0]] = [] ## Only set empty list if first time seen
            continue

        for fam in fams:
            if len(limits)>0 and not fam in limits: continue
            for fam2 in fams:
                if len(limits)>0 and not fam2 in limits: continue            
                if fam==fam2: continue
                try: a[fam].append(fam2)
                except KeyError: a[fam] = [fam2]        
    hin.close()

    for key,value in a.iteritems():
        fc1 = famcount[key] ## Number of occurrences of family key
        ## Calculate fraction of times this family is seen together with other families in the list
        u = list(set(value)) ## Create unique list of families
        keep = []
        for fam in u:
            if corr:
                try: c = corr[key][fam]
                except KeyError: c = 1
                if c<corrmin: continue
            fc2 = famcount[fam] ## Number of occurrences of current family
            c = value.count(fam) ## Number of times current family is seen with family key

            ## Calculate fraction based on the least abundant family
            if fc1>fc2: f = c/float(fc2)
            else: f = c/float(fc1)
            if f>=frac: keep.append(fam)
        a[key] = list(set(keep))

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

def write(am, filtered):
    hout = sys.stdout
    houtcsv = csv.writer(hout,delimiter = '\t')
    for tg,fams in am.iteritems():
        cogs = []
        tigrs = []
        pfams = []
        other = []
        for f in fams:
            if f in filtered: continue
            try:
                if f[0:4]=="TIGR": tigrs.append(f)
                elif f[0:3] == "COG": cogs.append(f)
                elif f[0:2] == "PF": pfams.append(f)
                else: other.append(f)
            except IndexError: continue
        houtcsv.writerow([tg,"|".join(pfams),"|".join(tigrs),"|".join(cogs),"|".join(other)])

def write_filtered(outfile, filtered_fams, a):
    edgecounts = {}
    for f in filtered_fams: edgecounts[f] = len(a[f])+1
    ## Sort by edgecount
    import operator
    sorted_edges = sorted(edgecounts.items(), key=operator.itemgetter(1), reverse=True)
    hout = open(outfile, 'w')
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
    parser.add_argument("--minfrac", type=float, default = 0.5,
            help="Minimum fraction of times that a link between two families has to occur. \
                    A higher value will decrease the number of transporters merged. Defaults to 0.5.")
    parser.add_argument("-c", "--corr", type=str, required=False,
            help="Read correlation matrix for families")
    parser.add_argument("--corrmin", type=float, default=0.5,
            help="Minimum correlation coefficient between families to merge. See the -c flag")
    parser.add_argument("--filterout", type=str,
            help="Write filtered families, with number of outgoing edges found, to specified outfile")
    args = parser.parse_args()

    if not args.infile: sys.exit(parser.print_help())
    
    ## Read family limits
    limits = readlimits(args.famfile)
    
    ## Read correlations if specified
    if args.corr: corr = readcorr(args.corr)
    else: corr = False

    ## Read annotations
    a = readannots(args.infile, args.minfrac, limits, corr, args.corrmin)

    ## Merge and also return families filtered by edgecount
    (am,filtered_fams) = mergeannots(a, args.edgecount)
    
    ## Write mergings
    write(am, filtered_fams)
    
    ## Write protein families with more outgoing edges than the limit
    if args.filterout: write_filtered(args.filterout, filtered_fams, a)

if __name__ == '__main__': 
    main()
