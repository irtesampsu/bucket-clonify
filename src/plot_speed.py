import pandas as pd
import matplotlib.pyplot as plt
import os

# Locate the speed summary CSV
candidates = ['results_summary.csv', 'speed_summary.csv', 'metrics/results_summary.csv']
csv_path = next((c for c in candidates if os.path.exists(c)), None)
if csv_path is None:
    raise FileNotFoundError(
        "Could not find a speed summary CSV ('results_summary.csv' or 'speed_summary.csv')"
    )

# Load data
df = pd.read_csv(csv_path)

# Ensure columns are correct
# Expected columns: experiment, out_dir, log_file, seconds, returncode
if 'experiment' not in df.columns or 'seconds' not in df.columns:
    raise ValueError("CSV must contain 'experiment' and 'seconds' columns")

# Sort by seconds for better visibility
df_sorted = df.sort_values('seconds')

# Create bar chart of runtime per experiment
plt.figure(figsize=(10, 6))
plt.bar(df_sorted['experiment'], df_sorted['seconds'])
plt.xticks(rotation=45, ha='right')
plt.ylabel('Runtime (seconds)')
plt.xlabel('Experiment')
plt.title('Clonify Bucketing Methods: Runtime per Experiment')
plt.tight_layout()

# Save and show
output_fig = 'fig_runtime_comparison.png'
plt.savefig(output_fig)
plt.show()

print(f"Saved runtime comparison plot to {output_fig}")