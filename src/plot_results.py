import math
import pandas as pd
import matplotlib.pyplot as plt
import os

# Locate the ARI summary CSV
candidates = ['ari_summary.csv']
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

# Plot one figure per sample
for sample, group in df.groupby('sample'):
    plt.figure(figsize=(8, 5))
    # sort by config for reproducibility
    group = group.sort_values('config')
    plt.plot(group['config'], group['ari'], marker='o')
    plt.xticks(rotation=45, ha='right')
    plt.title(f"ARI vs. Configuration for Sample: {sample}")
    plt.xlabel("Bucketing Configuration")
    plt.ylabel("Adjusted Rand Index (ARI)")
    plt.tight_layout()
    # Save each figure
    out_path = f"out/figure/fig_ari_{sample}.png"
    plt.savefig(out_path)
    plt.close()
    print(f"Saved figure for {sample}: {out_path}")

# Additionally, one figure per method showing ARI trend over parameters
for method, group in df.groupby('method'):
    plt.figure(figsize=(8, 5))
    for sample, samp_grp in group.groupby('sample'):
        # sort by params (lexicographically)
        samp_grp = samp_grp.sort_values('params')
        plt.plot(samp_grp['params'], samp_grp['ari'], marker='o', label=sample)
    plt.xticks(rotation=45, ha='right')
    plt.title(f"ARI Parameter Sweep for Method: {method}")
    plt.xlabel("Parameter Settings")
    plt.ylabel("Adjusted Rand Index (ARI)")
    plt.legend(title="Sample")
    plt.tight_layout()
    out_path = f"out/figure/fig_ari_{method}_sweep.png"
    plt.savefig(out_path)
    plt.close()
    print(f"Saved parameter-sweep figure for {method}: {out_path}")


# Prepare subplot grid
samples = sorted(df['sample'].unique())
n_samples = len(samples)
cols = 4  # you can adjust columns per row
rows = math.ceil(n_samples / cols)

fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), squeeze=False)

# Plot each sample in its own subplot
for idx, sample in enumerate(samples):
    r = idx // cols
    c = idx % cols
    ax = axes[r][c]
    grp = df[df['sample'] == sample].sort_values('config')
    ax.plot(grp['config'], grp['ari'], marker='o')
    ax.set_title(sample)
    ax.set_xticklabels(grp['config'], rotation=30, ha='right')
    ax.set_ylabel('ARI')
    ax.set_xlabel('Config')

# Hide any unused subplots
for idx in range(n_samples, rows * cols):
    r = idx // cols
    c = idx % cols
    axes[r][c].axis('off')

plt.tight_layout()
plt.savefig('fig_ari_all_samples.png')
# plt.show()