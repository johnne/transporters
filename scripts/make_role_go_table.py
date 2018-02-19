#!/usr/bin/env python

import sys
import pandas as pd
from argparse import ArgumentParser


def parse_roles(role_names):
    res = {}
    for i in role_names.index:
        r = role_names.loc[i]
        role_id = r["role_id"]
        name = r["role_name"]
        role_type = r["role_type"].rstrip(":")
        if not role_id in res.keys():
            res[role_id] = {}
        res[role_id][role_type] = name
    return pd.DataFrame(res).T


def main():
    parser = ArgumentParser()
    parser.add_argument("--tigr2name",
                        help="TIGRFAM to TIGRFAM name map")
    parser.add_argument("--tigr2role",
                        help="TIGRFAM to TIGRFAM role map")
    parser.add_argument("--role_names",
                        help="TIGRFAM Role names")
    parser.add_argument("--tigr2go",
                        help="TIGRFAM to Gene ontology map")
    parser.add_argument("--go2name",
                        help="Gene Ontology name map")
    parser.add_argument("--transporters",
                        help="Transport cluster to protein family map")
    args = parser.parse_args()
    tigr2name = pd.read_table(args.tigr2name, header=None, names=["TIGRFAM","TIGRFAM_NAME"])
    tigr2role = pd.read_table(args.tigr2role, header=None, names=["TIGRFAM", "ROLE_ID"])
    role_names = pd.read_table(args.role_names, header=None, names=["role_id", "role_type", "role_name"],
                               usecols=[1, 2, 3])
    role_names_df = parse_roles(role_names)

    tigr2go = pd.read_table(args.tigr2go, header=None, names=["TIGRFAM","go_id"], usecols=[0,1])
    go_names = pd.read_table(args.go2name, header=None, names=["go_id", "go_name"])

    trans_df = pd.read_table(args.transporters, header=None, names=["transporter","FAM"])

    # Add TIGR roles
    df_roles = pd.merge(tigr2role, role_names_df, left_on="ROLE_ID", right_index=True)

    # Add GO terms
    df_go = pd.merge(tigr2go,go_names, left_on="go_id", right_on="go_id")

    # Merge the roles and go frames
    df = pd.merge(df_roles, df_go, left_on="TIGRFAM", right_on="TIGRFAM", how="outer")
    df = pd.merge(df,tigr2name, left_on="TIGRFAM", right_on="TIGRFAM", how="left")
    # Add transporter info
    df = pd.merge(trans_df,df, left_on="FAM", right_on="TIGRFAM")
    df = df[["transporter","TIGRFAM","mainrole","sub1role","go_id","go_name","TIGRFAM_NAME"]]

    fh = sys.stdout
    df.to_csv(fh, sep="\t", index=False)
if __name__ == '__main__':
    main()