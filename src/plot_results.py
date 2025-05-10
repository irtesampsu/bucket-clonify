import math
import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

# === CONFIG ===
argparser = argparse.ArgumentParser(description='Plot ARI results from CSV.')
argparser.add_argument('--collection', type=str, default=None, help='MongoDB collection name')
args = argparser.parse_args()
if args.collection:
    print(f"Using collection: {args.collection}")
    COLLECTION_NAME = args.collection
    out_dir = f"out/{COLLECTION_NAME}/"
    os.makedirs(out_dir, exist_ok=True)
    candidates = [f"ari_summary_{COLLECTION_NAME}.csv"]
else:
    out_dir = "out/"
    candidates = ['ari_summary.csv']    
# Locate the ARI summary CSV

csv_path = next((c for c in candidates if os.path.exists(c)), None)
if csv_path is None:
    raise FileNotFoundError("Could not find 'ari_summary.csv' in the repo root or 'metrics/'")

# Load data
df = pd.read_csv(csv_path)

# For consistent ordering, create a combined label
df['config'] = df.apply(
    lambda row: f"{row.method}_{row.params}" if row.params else row.method,
    axis=1
)

def plot_metrics(metric):
    ylabel_dict = {"ari": f"Adjusted Rand Index (ARI)", "nmi": f"Normalized Mutual Information (NMI)"}
    # Plot one figure per sample
    for sample, group in df.groupby('sample'):
        plt.figure(figsize=(8, 5))
        # sort by config for reproducibility
        group = group.sort_values('config')
        plt.plot(group['config'], group[metric], marker='o')
        plt.xticks(rotation=45, ha='right')
        plt.title(f"{metric} vs. Configuration for Sample: {sample}")
        plt.xlabel("Bucketing Configuration")
        plt.ylabel(ylabel_dict[metric])
        plt.tight_layout()
        # Save each figure
        fig_dir = os.path.join(out_dir, "figure")
        os.makedirs(fig_dir, exist_ok=True)
        out_path = f"{fig_dir}/fig_{metric}_{sample}.png"
        plt.savefig(out_path)
        plt.close()
        print(f"Saved figure for {sample}: {out_path}")

    # Additionally, one figure per method showing ARI trend over parameters
    for method, group in df.groupby('method'):
        plt.figure(figsize=(8, 5))
        for sample, samp_grp in group.groupby('sample'):
            # sort by params (lexicographically)
            samp_grp = samp_grp.sort_values('params')
            
            plt.plot(samp_grp['params'], samp_grp[metric], marker='o', label=sample)
        plt.xticks(rotation=90, ha='right')
        plt.title(f"{metric} Parameter Sweep for Method: {method}")
        plt.xlabel("Parameter Settings")
        plt.ylabel(ylabel_dict[metric])
        plt.legend(title="Sample")
        plt.tight_layout()
        out_path = f"{fig_dir}/fig_{metric}_{method}_sweep.png"

        plt.savefig(out_path)
        plt.close()
        print(f"Saved parameter-sweep figure for {method}: {out_path}")


    # Prepare subplot grid
    samples = sorted(df['sample'].unique())
    n_samples = len(samples)
    cols = 4  # you can adjust columns per row
    rows = math.ceil(n_samples / cols)

    if n_samples <= 1:
        print("No samples found in the CSV.")
        # exit(1)
        return

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), squeeze=False)

    # Plot each sample in its own subplot
    for idx, sample in enumerate(samples):
        r = idx // cols
        c = idx % cols
        ax = axes[r][c]
        grp = df[df['sample'] == sample].sort_values('config')
        ax.plot(grp['config'], grp[metric], marker='o')
        ax.set_title(sample)
        xticks = grp['config'].tolist()
        xticks = [x.replace("_", " ").replace("-", "=").title() for x in xticks]
        ax.set_xticklabels(xticks, rotation=90, ha='right')
        ax.set_ylabel(ylabel_dict[metric])
        ax.set_xlabel('Config')

    # Hide any unused subplots
    for idx in range(n_samples, rows * cols):
        r = idx // cols
        c = idx % cols
        axes[r][c].axis('off')

    plt.tight_layout()
    plt.savefig(f'fig_{metric}_all_samples.png')
    plt.close()

# plt.show()


plot_metrics("ari")
plot_metrics("nmi")
print("Plots saved successfully.")