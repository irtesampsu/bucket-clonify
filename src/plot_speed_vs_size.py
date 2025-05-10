import pandas as pd
import matplotlib.pyplot as plt

# ← change this to your actual CSV path if needed
FILENAME = 'large_results_summary.csv'

# Load your experiment summary
df = pd.read_csv(FILENAME)

# Parse out sample size (n) and method
df['sample_size'] = df['experiment'].str.extract(r'n(\d+)').astype(int)
df['method']      = df['experiment'].apply(
    lambda x: 'baseline' if x.startswith('baseline') else 'bktree'
)

# Ensure timing is numeric
df['seconds'] = df['seconds'].astype(float)

# Plot both curves on one chart
plt.figure()
for method, grp in df.groupby('method'):
    plt.plot(grp['sample_size'], grp['seconds'], marker='o', label=method)

plt.xlabel('Sample Size (N)')
plt.ylabel('Running Time (seconds)')
plt.title('Clonify3: baseline vs. bktree (th=7)')
plt.legend()
plt.tight_layout()
plt.savefig('speed_vs_size.png')
