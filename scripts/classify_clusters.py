#!/usr/bin/env python

import pandas as pd, sys, logging
from argparse import ArgumentParser

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def read_map(f,godf):
    m = {}
    with open(f, 'r') as fh:
        for line in fh:
            line = line.rstrip()
            fam = line.split(" ")[0].split(":")[-1]
            term = "GO:"+line.split(":")[-1]
            try: godf.loc[term]
            except KeyError: continue
            try: m[fam].append(term)
            except KeyError: m[fam] = [term]
    return m

def count_terms(t,m,fams):
    terms = {}
    df = pd.DataFrame(columns=["term","tigr","pf","cog"])
    for fam in fams:
        if not fam in m: continue
        for term in m[fam]: 
            tigr = pfam = cog = 0
            if fam[0:2]=="PF": pfam = 1
            elif fam[0:2]=="TI": tigr=1
            elif fam[0:2]=="CO": cog = 1
            tmp = pd.DataFrame(columns=df.columns,data={"term": term, "tigr": tigr, "pf": pfam, "cog": cog}, index=[t])
            df = pd.concat([df,tmp])
            try: terms[term]+=1
            except KeyError: terms[term] = 1
    #df = pd.DataFrame([terms],index=[t]).T
    #df.sort_values(t,ascending=False, inplace=True)
    s = df.groupby("term").sum()
    term_sum = s.sum(axis=1)
    df = pd.concat([s,term_sum],axis=1)
    df.columns = ["tigr","pf","cog",t]
    return df.sort_values(t,ascending=False)
    #print(df.groupby("term").sum())
    #return df

def classify(tdf,godf,m):
    c = {}
    ambiguous = unambiguous = unclassified = 0
    for i,t in enumerate(tdf.index.unique()):
        c[t] = {0: "", 1: "", 2: "", 3: "", 4: "", 5: "", 6: "", 7: "", 8: ""}
        fams = tdf.loc[t,'Family']
        if type(fams)==str: fams = [fams]
        terms = count_terms(t,m,fams)
        terms = pd.merge(terms,godf.loc[terms.index],left_index=True,right_index=True)
        #continue
        ## Get the count of the term with the highest count
        max_count = terms[t].max()
        ## Get all terms with that count
        term_ids = list(set(terms.loc[terms[t]==max_count].index))
        if term_ids == []:
            for key in c[t].keys(): c[t][key] = "UNCLASSIFIED"
            unclassified +=1
            continue
        ## If more than one term_id with the max count, then first choose terms from TIGRFAMs, then COG, then PFAM
        if len(term_ids)>1:
            term_ids_ti = list(set(terms.loc[terms["tigr"]==max_count].index))
            if term_ids_ti==[]:
                term_ids_co = list(set(terms.loc[terms["cog"]==max_count].index))
                if term_ids_co==[]: 
                    term_ids_pf = list(set(terms.loc[terms["pf"]==max_count].index))
                    if term_ids_pf == []: term_ids = term_ids
                    else: term_ids = term_ids_pf
                else: term_ids = term_ids_co
            else: term_ids = term_ids_ti
    
        term_names = "/".join(sorted(list(set(terms.loc[term_ids,"name"]))))
        parents = [p.split("|") for p in terms.loc[term_ids,"parent_names"]]
        parentdf = pd.DataFrame(parents)
        ambig = False
        for col in parentdf.columns:
            s = parentdf.loc[:,col]
            s = s[s.notnull()]
            col_names = list(set(s))
            if len(col_names)>1: ambig = True
            c[t][col] = "/".join(col_names)
        if ambig: ambiguous+=1
        else: unambiguous+=1
        c[t][6] = term_names
        for col in list(range(1,9)):
            if c[t][col]=="": 
                if c[t][col-1][:13]=="Unclassified.": c[t][col] = c[t][col-1]
                else: c[t][col] = "Unclassified."+c[t][col-1]
    logging.info(str(len(set(tdf.index)))+" clusters.")
    logging.info(str(unambiguous)+" unambiguous classifications")
    logging.info(str(ambiguous)+ " ambiguous classifications")
    logging.info(str(unclassified) + " unclassified clusters")
    return pd.DataFrame(c).T

def main():

    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
            help="Transporter cluster table")
    parser.add_argument("-g", "--gotable", required=True,
            help="Gene ontology table of terms and hierarchies")
    parser.add_argument("-m", "--mapfile", required=True,
            help="Gene ontology mapping between TIGRFAMs, COGs and PFAMs")

    args = parser.parse_args()

    tdf = pd.read_csv(args.infile, header=None, sep="\t", index_col=0, names=['Family'])
    godf = pd.read_csv(args.gotable, header=0, sep="\t", index_col=0)

    m = read_map(args.mapfile,godf)

    tdf_c = classify(tdf,godf,m)
    tdf_c = tdf_c.loc[tdf.index.unique()]

    tdf_c.columns = ["Category"+str(x) for x in tdf_c.columns]
    tdf_c.to_csv(sys.stdout, sep="\t")

if __name__ == '__main__':
    main()
