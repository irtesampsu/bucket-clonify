#!/usr/bin/env python3
"""
Compare every bucketing run’s assignments to the baseline,
computing ARI per sample and saving a summary table.
"""
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from pathlib import Path
import re
import argparse

argparser = argparse.ArgumentParser(description='Compare bucketing assignments.')
argparser.add_argument('--collection', type=str, default=None, help='MongoDB collection name')
args = argparser.parse_args()
if args.collection:
    print(f"Using collection: {args.collection}")
    OUT_DIR = Path("out")/args.collection
else:
    OUT_DIR = Path("out")



# OUT_DIR = Path("out")
SUMMARY_CSV = Path(f"ari_summary_{args.collection}.csv")

def load_assignments(path):
    df = pd.read_csv(path)
    df.columns = ["seq_id","clone_id"]
    return df.set_index("seq_id")["clone_id"]

def parse_out_path(p):
    """
    Given a path like out/faiss_k7_n32/sample1_assignments.csv
    return (method, params, sample):
      ("faiss", "k7_n32", "sample1")
    """
    parts = p.parts
    # parts[-2] == "<method>_<params>"
    m = parts[-2]
    sample = p.name.replace("_assignments.csv","")
    if "_" in m:
        method, params = m.split("_", 1)
    else:
        method, params = m, ""
    return method, params, sample

def main():
    # find baseline files
    baseline_files = list(OUT_DIR.glob("baseline*/*_assignments.csv"))
    if not baseline_files:
        raise RuntimeError("No baseline assignments found under out/baseline*/")
    # map sample -> baseline Series
    baseline = {}
    for p in baseline_files:
        _, _, sample = parse_out_path(p)
        baseline[sample] = load_assignments(p)

    rows = []
    # iterate all other assignment files
    for p in OUT_DIR.glob("*/*_assignments.csv"):
        method, params, sample = parse_out_path(p)
        if method.startswith("baseline") :# or method == "bktree":
            continue
        # if method != "bktree":
        #     continue
        if sample not in baseline:
            print(f"⚠️  no baseline for sample {sample}, skipping {p}")
            continue

        y_true = baseline[sample]
        y_pred = load_assignments(p)
        # intersect seq_ids
        idx = y_true.index.intersection(y_pred.index)
        # print non-matching seq_ids
        # if len(idx) != len(y_true) or len(idx) != len(y_pred):
        # print(f"⚠️  {p} has {len(y_true)-len(idx)} seq_ids only in baseline")
        # print(f"⚠️  {p} has {len(y_pred)-len(idx)} seq_ids only in predictions")
        # print(f"⚠️  {p} has {len(idx)} seq_ids in common")
        # print(f"⚠️  {p} has {len(y_true)} seq_ids in total in baseline")
        # print(f"⚠️  {p} has {len(y_pred)} seq_ids in total in predictions")
        # print(f"⚠️  {y_pred.index.difference(idx)}")
        # print(f"⚠️  {y_pred.index.has_duplicates}")
            # print(idx)
        # for i in y_pred.index:
        #     # if i not in y_pred.index:
        #     #     print(f"⚠️  {i} not in {p}")
        #     if i not in y_true.index:
        #         print(f"⚠️  {i} not in baseline")
        
        # print("Calculating ARI...")
        ari = adjusted_rand_score(y_true.loc[idx], y_pred.loc[idx])
        nmi = normalized_mutual_info_score(y_true.loc[idx], y_pred.loc[idx])
        
        rows.append({
            "method": method,
            "params": params,
            "sample": sample,
            "ari": ari,
            "nmi": nmi
        })

    # write summary
    df = pd.DataFrame(rows)
    df.to_csv(SUMMARY_CSV, index=False)
    print(f"Wrote ARI summary to {SUMMARY_CSV}")

if __name__ == "__main__":
    main()
