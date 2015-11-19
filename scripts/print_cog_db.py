#!/usr/bin/env python

import urllib2, sys

def main():
    r = urllib2.urlopen("ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog")
    for line in r.readlines():
        try: 
            if line[0]=="[": 
                line = line.rstrip()
                items = line.rsplit()
                print items[1]+"\t"+" ".join(items[2:])
        except IndexError: continue

if __name__ == '__main__':
    main()
