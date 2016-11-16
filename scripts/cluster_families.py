#!/usr/bin/env python

import networkx as nx, pandas as pd, sys, csv, logging
from argparse import ArgumentParser

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

def trim_graph(g,maxlink,minoc,oc):
    trimmed_nodes = []
    trimmed_edges = []
    for e in g.edges():
        (n1,n2) = e
        if oc[n1][n2]<minoc: trimmed_edges.append((n1,n2))
    g.remove_edges_from(trimmed_edges)

    for n in g.nodes():
        edges = g[n].keys()
        if len(edges)>maxlink: 
            g.remove_node(n) 
            trimmed_nodes.append(n)
    g.remove_edges_from(trimmed_edges)
    if "" in trimmed_nodes: trimmed_nodes.remove("")
    return [g,list(set(trimmed_nodes)),list(set(trimmed_edges))]

def cluster(g):
    clustered = []
    clusters = {}
    i = 1
    for n in g.nodes():
        c = [n]
        if n in clustered: continue
        edges = list(nx.dfs_edges(g,n))
        for e in edges:
            n1,n2 = e
            clustered+=[n1,n2]
            c+=[n1,n2]
        c = list(set(c))
        clusters[i] = {'num': len(c), 'fams': c[:]}
        i+=1
    return clusters

def count_occurrences(df):
    occur = {}
    for fam in list(set(df.Node1)):
        r = df.loc[df.Node1==fam]
        if list(set(r.Node2))[0]=='': 
            occur[fam] = {}
            continue
        tmp = df.loc[df.Node1==fam].groupby("Node2").count().ix[:,0]   
        occur[fam] = tmp.to_dict()
    return occur

def write(cdf):
    for i in cdf.index:
        clust = "T"+str(i)
        fams = cdf.loc[i,"fams"]
        for fam in fams: sys.stdout.write(clust+"\t"+fam+"\n")

def main():
    parser = ArgumentParser()
    parser.add_argument("-i", "--infile", type=str,
            help="Table with one row per link between families, columns: ['Family1','Family2'].\
                    If not specified the program reads from stdin.")
    parser.add_argument("--maxlink", default=6, type=int,
            help="Maximum number of allowed outgoing links (edges) for a single protein family. Defaults to 6")
    parser.add_argument("--minoc", default=10, type=int,
            help="Minimum number of occurrences for linking two families. Defaults to 10")
    parser.add_argument("--trimmed_out", type=str,
            help="Write trimmed families to file")

    args = parser.parse_args()
    
    if args.infile: linkdf = pd.read_csv(args.infile, sep="\t", header=0)
    else: linkdf = pd.read_csv(sys.stdin, sep="\t", header=0)
    linkdf.fillna("",inplace=True)
    oc = count_occurrences(linkdf)
    
    ## Create graph from data frame
    g = nx.from_pandas_dataframe(linkdf,source="Node1",target="Node2")
    g.remove_node('')
    ## Trim nodes by outgoing edges and edges by occurrence
    [gt,trimmed_nodes,trimmed_edges] = trim_graph(g,args.maxlink, args.minoc, oc)
    if args.trimmed_out: pd.DataFrame(trimmed_nodes).to_csv(args.trimmed_out,sep="\t",index=False,header=False)
    logging.info("Removed "+str(len(trimmed_edges))+" links due to low occurrence")    
    logging.info("Removed "+str(len(trimmed_nodes))+" families with too many links ("+str(len(gt.nodes()))+" remaining)")
    
    ## Create clusters for graph
    clusters = cluster(gt)
    logging.info(str(len(clusters))+" clusters created")    
    cdf = pd.DataFrame(clusters).T
    ## Sort by number of families in cluster
    cdf.sort_values("num",ascending=False,inplace=True)
    cdf.index = list(range(1,len(cdf)+1))
    
    ## Write clusters sorted by size
    write(cdf)
if __name__ == '__main__':
    main()

