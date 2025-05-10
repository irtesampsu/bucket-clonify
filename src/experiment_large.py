#!/usr/bin/env python3
import subprocess
import time
import csv
from pathlib import Path

# === CONFIG ===
PYTHON = "python -W ignore"
SCRIPT  = "src/clonify3.py"
COLLECTION = "SRR8283831"
COMMON_ARGS = [
    "--db", "clonify_db",
    "--ip", "localhost",
    "--port", "27017",
    "--split_by", "gene",
    "--threads", "20",
    "--coll", COLLECTION,
    "--nt",
    "--no_update",
    "--verbose"
]

# only one bucketing method
BUCKET_ARGS = ["-b", "bktree", "--threshold", "7"]

# sample sizes from 10000 to 70000 in steps of 10000
SAMPLE_SIZES = list(range(10000, 70000, 10000))

# build experiments
EXPERIMENTS = []
for n in SAMPLE_SIZES:
    EXPERIMENTS.append({
        "name": f"baseline_n{n}",
        "args": ["--sample-size", str(n)]
    })
    EXPERIMENTS.append({
        "name": f"bktree_th7_n{n}",
        "args": ["--sample-size", str(n)] + BUCKET_ARGS
    })

summary_path = Path("large_results_summary.csv")
with summary_path.open("w", newline="") as csvf:
    writer = csv.writer(csvf)
    writer.writerow(["experiment","out_dir","log_file","seconds","returncode"])

    for exp in EXPERIMENTS:
        name    = exp["name"]
        out_dir = Path("out")/ COLLECTION/name
        log_f   = Path("logs")/COLLECTION/f"{name}.log"
        out_dir.mkdir(parents=True, exist_ok=True)
        log_f.parent.mkdir(exist_ok=True)

        cmd = (
            f"{PYTHON} {SCRIPT} "
            + " ".join(COMMON_ARGS)
            + " --out " + str(out_dir) + " "
            + " ".join(exp["args"])
        )

        print(f"[RUNNING] {name} --> log: {log_f}")
        t0 = time.perf_counter()
        ret = subprocess.run(cmd + f" > {log_f} 2>&1", shell=True)
        elapsed = time.perf_counter() - t0

        writer.writerow([name, str(out_dir), str(log_f), f"{elapsed:.1f}", ret.returncode])
        csvf.flush()

        if ret.returncode != 0:
            print(f"ERROR: {name} exited {ret.returncode}, see {log_f}")

print(f"All done. Summary written to {summary_path}")
