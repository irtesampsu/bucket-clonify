import argparse
import csv
from pymongo import MongoClient

def config_db(serverip = "localhost",port=27017, db_name = "clonify_db", collection_name = "antibodies_7k"):
    # Parameters
    client = MongoClient(serverip, port)
    db = client[db_name]
    collection = db[collection_name]
    return client, collection


def get_nt_mutations(seq_align, germ_align):
    mutations = []
    for i, (s, g) in enumerate(zip(seq_align, germ_align)):
        if s != g and g != 'N':  # Ignore masked/unknown bases
            mutations.append({'loc': i + 1, 'mut': f'{g}>{s}'})
    return mutations



def load_airr_data(tsv_path, collection):
    """
    Load AIRR data from a TSV file into MongoDB.
    Args:
        tsv_path (str): Path to the TSV file.
        collection (pymongo.collection.Collection): MongoDB collection to insert data into.
    """
    row_processed = 0
    with open(tsv_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            doc = {
                'seq_id': row['sequence_id'],
                'v_gene': { 'full': row['v_call'] },
                'j_gene': { 'full': row['j_call'] },
                'junc_aa': row.get('junction_aa', ''),
                'junc_nt': row.get('junction', ''),
                'var_muts_nt': {
                    'muts': get_nt_mutations(row['sequence_alignment'], row['germline_alignment'])
                },
                'chain': row.get('locus', 'heavy')
            }
            collection.insert_one(doc)
            row_processed += 1
            if row_processed % 1000 == 0:
                print(f"Processed {row_processed} rows...")

def main():
    # take command line arguments
    parser = argparse.ArgumentParser(description='Load AIRR data into MongoDB.')
    parser.add_argument('--tsv_path', type=str, help='Path to the AIRR TSV file')
    parser.add_argument('--serverip', type=str, required=True, help='MongoDB server IP address')
    parser.add_argument('--port', type=int, default=27017, help='MongoDB server port') # default is 27017
    parser.add_argument('--db', type=str, default='clonify_db', help='MongoDB database name')
    parser.add_argument('--collection', type=str, default='antibodies_7k', help='MongoDB collection name')

    print("Parsing arguments...")
    args = parser.parse_args()
    tsv_path = args.tsv_path 
    serverip = args.serverip 
    port = args.port
    db_name = args.db
    collection_name = args.collection

    print(f"Loading AIRR data from {tsv_path} into MongoDB...")
    # Check if the TSV file exists
    try:
        with open(tsv_path, 'r') as f:
            pass
    except FileNotFoundError:
        print(f"File {tsv_path} not found.")
        return
    
    # Configure database
    client, collection = config_db(serverip, port, db_name, collection_name)
    print(f"Connected to MongoDB at {serverip}:{port}, database: {db_name}, collection: {collection_name}")

    # Load AIRR data into MongoDB
    load_airr_data(tsv_path, collection)

    # Close the database connection
    client.close()

print("Done loading AIRR data into MongoDB.")

if __name__ == "__main__":
    main()
