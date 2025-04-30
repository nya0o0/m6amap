#!/usr/bin/env python

import pandas as pd
import glob


def load_m6add_data(datatype):
    """
    Load and preprocess m6ADD dataset files.
    """
    files = glob.glob(f"tests/data/GSE*_{datatype}.txt")
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        dfs.append(df)
    df = pd.concat(dfs)
    df["chr"] = df["chr"].astype(str).str.lower().str.replace("chr", "")
    df["genome_pos"] = ((df["start"] + df["end"]) / 2).astype(int)
    return df

def fuzzy_merge(m6anet_df, m6add_df, window=3):
    """
    Perform fuzzy matching between m6anet and m6add dataframes based on genomic position.
    """
    matches = []
    for _, row in m6anet_df.iterrows():
        chr_match = m6add_df["chr"] == str(row["chromosome"])
        pos_match = (m6add_df["genome_pos"] >= row["genome_pos"] - window) & \
                    (m6add_df["genome_pos"] <= row["genome_pos"] + window)
        match_df = m6add_df[chr_match & pos_match]
        for _, match_row in match_df.iterrows():
            combined = row.to_dict()
            combined.update({f"m6add_{col}": match_row[col] for col in m6add_df.columns})
            matches.append(combined)
    return pd.DataFrame(matches)

