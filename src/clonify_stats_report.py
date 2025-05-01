
import pymongo
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

import matplotlib
matplotlib.use('Agg')


# MongoDB connection details
DB_NAME = "clonify_db"  # Replace with actual DB name
IP = "localhost"
PORT = 27017
JUNC_FIELD = "junc_nt"  # or "junc_nt"

# Connect to MongoDB
client = pymongo.MongoClient(IP, PORT)
db = client[DB_NAME]
collections = db.list_collection_names()

summary_stats = []
all_data = []

# Create output directory
plot_dir = "out/plots"
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)


for coll_name in collections:
    print(f"Processing collection: {coll_name}")
    coll = db[coll_name]
    cursor = coll.find({'chain': 'IGH'}, {
        '_id': 0, 'seq_id': 1, 'v_gene.full': 1, 'j_gene.full': 1,
        JUNC_FIELD: 1, 'var_muts_nt.muts': 1
    })

    data = []
    for doc in cursor:
        if JUNC_FIELD not in doc:
            continue
        junc_len = len(doc[JUNC_FIELD])
        mut_count = len(doc.get('var_muts_nt', {}).get('muts', []))
        v_gene = doc['v_gene']['full']
        j_gene = doc['j_gene']['full']
        data.append({
            'sample': coll_name,
            'junction_length': junc_len,
            'mutation_count': mut_count,
            'v_gene': v_gene,
            'j_gene': j_gene
        })

    if not data:
        continue

    df = pd.DataFrame(data)
    all_data.append(df)

    # Summary statistics
    summary_stats.append({
        'sample': coll_name,
        'num_seqs': len(df),
        'junction_len_mean': df['junction_length'].mean(),
        'junction_len_std': df['junction_length'].std(),
        'junction_len_median': df['junction_length'].median(),
        'mutation_mean': df['mutation_count'].mean(),
        'mutation_std': df['mutation_count'].std(),
        'unique_v_genes': df['v_gene'].nunique(),
        'unique_j_genes': df['j_gene'].nunique()
    })

    # Plot junction length distribution
    plt.figure()
    sns.histplot(df['junction_length'], kde=True, bins=20)
    plt.title(f"Junction Length Distribution: {coll_name}")
    plt.xlabel("Junction Length")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{coll_name}_junction_length.png")
    plt.close()

    # Plot mutation count distribution
    plt.figure()
    sns.histplot(df['mutation_count'], bins=20)
    plt.title(f"Mutation Count Distribution: {coll_name}")
    plt.xlabel("Mutation Count")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{coll_name}_mutation_count.png")
    plt.close()

        # V gene usage barplot
    plt.figure(figsize=(10, 4))
    v_counts = df['v_gene'].value_counts().head(20)  # top 20 V genes
    sns.barplot(x=v_counts.values, y=v_counts.index)
    plt.title(f"Top V Gene Usage: {coll_name}")
    plt.xlabel("Count")
    plt.ylabel("V Gene")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{coll_name}_v_gene_usage.png")
    plt.close()

    # J gene usage barplot
    plt.figure(figsize=(6, 4))
    j_counts = df['j_gene'].value_counts()
    sns.barplot(x=j_counts.values, y=j_counts.index)
    plt.title(f"J Gene Usage: {coll_name}")
    plt.xlabel("Count")
    plt.ylabel("J Gene")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{coll_name}_j_gene_usage.png")
    plt.close()

    # V-J gene combination heatmap
    pivot = df.groupby(['v_gene', 'j_gene']).size().unstack(fill_value=0)
    plt.figure(figsize=(12, 8))
    sns.heatmap(pivot, cmap="YlGnBu", annot=False)
    plt.title(f"V-J Gene Combination Heatmap: {coll_name}")
    plt.xlabel("J Gene")
    plt.ylabel("V Gene")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{coll_name}_vj_heatmap.png")
    plt.close()


# Combine all data and summary
full_df = pd.concat(all_data, ignore_index=True)
summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv("clonify_sample_summary_stats.csv", index=False)
print("Summary stats saved to clonify_sample_summary_stats.csv")
print("Plots saved in the 'out/plots/' directory.")
