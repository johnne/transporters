#!/usr/bin/env python

import sys

def print_hmm(hmm):
    print hmm['ACC'].split(".")[0]+"\t"+hmm['DESC']
    hmm = {"ACC": "", "DESC": "", "//": ""}
    return hmm    

def main():
    try: infile = sys.argv[1]
    except IndexError: sys.exit("Usage: print_hmm_db.py <infile.hmm>")
    
    if ".gz" in infile:
        import gzip as gz
        hin = gz.open(infile)
    else: hin = open(infile)
    hmm = {"ACC": "", "DESC": "", "//": ""}
    for line in hin:
        line = line.rstrip()
        items = line.rsplit()
        try: hmm[items[0]]
        except KeyError: continue
        if items[0]=="//": hmm = print_hmm(hmm)
        else: hmm[items[0]] = " ".join(items[1:])
    hin.close()

if __name__ == '__main__':
    main()

