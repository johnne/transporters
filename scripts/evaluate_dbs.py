#!/usr/bin/env python

from argparse import ArgumentParser
import sys,csv

def CheckUniqueMap(fams,t,d):
    for fam in fams:
        if fam=='':continue
        try: d[fam].append(t)
        except KeyError: d[fam] = [t]
    return d

def IdentFams(item,pfams,tigrs,cogs):
    '''Identifies the type of protein family database'''
    fams = item.split("|")
    for fam in fams: 
        if fam[0:2] == "PF": return [fams,tigrs,cogs]
        elif fam[0:4] == "TIGR": return [pfams,fams,cogs]
        elif fam[0:3] == "COG": return [pfams,tigrs,fams]

def Readdb(infile):
    '''Reads the transporter definitions'''
    hin = open(infile, 'r')
    hincsv = csv.reader(hin, delimiter = '\t')
    db = {}
    fu = {} ## Store unique family -> transport mapping
    for row in hincsv:
        t = row[0]
        pfams = []
        tigrs = []
        cogs = []
        for item in row[1:]:
            if item == '': continue
            [pfams,tigrs,cogs] = IdentFams(item,pfams,tigrs,cogs)
        fu = CheckUniqueMap(pfams+tigrs+cogs,t,fu)
        db[t] = {"PFAM": len(pfams), "TIGRFAM": len(tigrs), "COG": len(cogs)}
    return (db,fu)

def Evaldb(db,fu):
    '''Calculates and prints simple statistics of the protein family database definitions'''
    ## 1. Calculate number of transporters unique to each database
    ## 2. Calculate total number of transporters defined by each database
    ## 3. Calculate average number of families per transporter for each database

    ## Store definitions per database
    pf = []
    ti = []
    co = []

    ## Store unique definitions
    pf_u = []
    ti_u = []
    co_u = []

    ## Store the sum of families defined
    pf_sum = 0
    ti_sum = 0
    co_sum = 0
    
    for t,counts in db.iteritems():
        pf_sum+=counts["PFAM"]
        ti_sum+=counts["TIGRFAM"]
        co_sum+=counts["COG"]
        
        if counts["PFAM"] > 0: pf.append(t)
        if counts["TIGRFAM"] > 0: ti.append(t)
        if counts["COG"] > 0: co.append(t)

        if counts["PFAM"] == 0 and counts["TIGRFAM"] == 0: co_u.append(t)
        if counts["PFAM"] == 0 and counts["COG"] == 0: ti_u.append(t)
        if counts["TIGRFAM"] == 0 and counts["COG"] == 0: pf_u.append(t)
    pf_av = pf_sum/float(len(pf))
    ti_av = ti_sum/float(len(ti))
    co_av = co_sum/float(len(co))

    print "Unique/Total/Average per transporter for PFAM:",str(len(pf_u))+"/"+str(len(pf))+"/"+str(round(pf_av,2))
    print "Unique/Total/Average per transporter for TIGRFAM:",str(len(ti_u))+"/"+str(len(ti))+"/"+str(round(ti_av,2))
    print "Unique/Total/Average per transporter for COG:",str(len(co_u))+"/"+str(len(co))+"/"+str(round(co_av,2))

    
    print
    ## 4. Check unique family -> transporter mapping
    unique = True
    for fam,ts in fu.iteritems():
        if len(ts)>1: 
            print "Non-unique family -> transporter map:",fam,"|".join(ts)
            unique = False
    if unique: print "All family -> transporter maps are unique"

def main():
    parser = ArgumentParser()
    parser.add_argument("-d", "--db", type=str, required=True,
            help="Transporter definitions, tab separated. Format: <Transporter> <Other db> <PFAMs> <TIGRFAMs> <COGs>")
    args = parser.parse_args()

    (db,fu) = Readdb(args.db)
    Evaldb(db,fu)

if __name__ == '__main__': main()
