#!/usr/bin/env python3
import subprocess
import time
import csv
import os
from pathlib import Path

# === CONFIG ===
PYTHON = "python -W ignore"
SCRIPT  = "src/clonify3.py"
COMMON_ARGS = [
    "--db", "clonify_db",
    "--ip", "localhost",
    "--port", "27017",
    "--split_by", "gene",
    "--threads", "20",
    "--nt",
    "--no_update",
]

# define your experiments
EXPERIMENTS = []

# 1) baseline (no bucketing)
EXPERIMENTS.append({
    "name": "baseline",
    "args": []
})

# 2) faiss: vary kmer
for k in [5,7,9]:
    EXPERIMENTS.append({
        "name": f"faiss_k-{k}",
        "args": ["-b", "faiss", "--kmer", str(k)]
    })

# 3) minhash: vary kmer × nperm
for k in [5,7,9]:
    for nperm in [16,32,64]:
        EXPERIMENTS.append({
            "name": f"minhash_k-{k}_n-{nperm}",
            "args": ["-b", "minhash", "--kmer", str(k), "--nperm", str(nperm)]
        })

# 4) bktree: vary threshold
for th in [3,5,7]:
    EXPERIMENTS.append({
        "name": f"bktree_th-{th}",
        "args": ["-b", "bktree", "--threshold", str(th)]
    })

# output summary file
summary_path = Path("results_summary.csv")
with summary_path.open("w", newline="") as csvf:
    writer = csv.writer(csvf)
    writer.writerow(["experiment","out_dir","log_file","seconds","returncode"])

    for exp in EXPERIMENTS:
        name    = exp["name"]
        out_dir = Path("out")/name
        log_f   = Path("logs")/f"{name}.log"
        out_dir.mkdir(parents=True, exist_ok=True)
        log_f.parent.mkdir(exist_ok=True)

        cmd = (
            f"{PYTHON} {SCRIPT} "
            + " ".join(COMMON_ARGS) + " "
            + "--out " + str(out_dir) + " "
            + " ".join(exp["args"])
        )

        print(f"[RUNNING] {name} --> log: {log_f}")
        t0 = time.perf_counter()
        # spawn shell so redirection works
        ret = subprocess.run(cmd + f" > {log_f} 2>&1", shell=True)
        elapsed = time.perf_counter() - t0

        writer.writerow([name, str(out_dir), str(log_f), f"{elapsed:.1f}", ret.returncode])
        csvf.flush()

        if ret.returncode != 0:
            print(f"ERROR:  {name} exited {ret.returncode}, see {log_f}")

print(f"All done. Summary written to {summary_path}")
