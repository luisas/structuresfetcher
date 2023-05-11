#!/usr/bin/env python3

import pandas as pd
import sys
import csv


hits = sys.argv[1]


hits_df = pd.read_csv(hits, sep = "\t", header = None)
hits_df["name"] = hits_df[0].str.replace("/", "_")
hits_df =hits_df.rename(columns = {8:"start_cut", 9: "end_cut"})
hits_df["line"] = hits_df.apply(lambda row : "extract_structure.py "+row["name"]+" "+str(row["start_cut"])+" "+str(row["end_cut"])+" "+(row["name"]), axis=1)
hits_df["line"].to_csv("cut_structures_tmp.sh", sep=" ", header=None, index=False, quoting=csv.QUOTE_NONE, escapechar=",")