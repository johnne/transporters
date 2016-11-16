#!/usr/bin/env python

from urllib import request
import sys

def main():
    r = request.urlopen("ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab")
    for line in r.readlines(): 
        line = line.decode('utf-8','ignore')
        if line[0]=="#": continue
        line = line.rstrip()
        items = line.split("\t")
        print(items[0]+"\t"+items[-1])

if __name__ == '__main__':
    main()
