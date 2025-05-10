import os
import argparse
import csv
from pymongo import MongoClient

MAX_ROWS = 70000

def config_db(serverip = "localhost",port=27017, db_name = "clonify_db", collection_name = "antibodies_7k"):
    # Parameters
    client = MongoClient(serverip, port)
    db = client[db_name]
    collection = db[collection_name]
    # Check if the collection exists
    if collection_name in db.list_collection_names():
        print(f"Collection {collection_name} already exists.")
        if collection_name == "SRR8283831":
            # Drop the collection if it exists
            print(f"Dropping collection {collection_name}...")
            db.drop_collection(collection_name)
            # Uncomment the next line to drop the collection
        # db.drop_collection(collection_name)
        return client, collection
        print(f"Collection {collection_name} dropped.")

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
    skipped_rows = 0
    with open(tsv_path, 'r') as f:
        if tsv_path.endswith('.csv'):
            reader = csv.DictReader(f, delimiter=',')
        else:
            reader = csv.DictReader(f, delimiter='\t')
            
        for row in reader:
            v_gene_first = row['v_call'].split(',')[0].strip()
            j_gene_first = row['j_call'].split(',')[0].strip()
            # Skip rows with empty junction sequences
            if not row['junction_aa'] or not row['junction']:
               skipped_rows += 1
               continue
            # Skip rows with missing V or J genes
            if not v_gene_first or not j_gene_first:
                skipped_rows += 1
                continue
            doc = {
                'seq_id': row['sequence_id'] if 'sequence_id' in row else str(row_processed),
                'v_gene': {'full': v_gene_first},
                'j_gene': {'full': j_gene_first},
                'junc_aa': row.get('junction_aa', ''),
                'junc_nt': row.get('junction', ''),
                'var_muts_nt': {
                    'muts': get_nt_mutations(row['sequence_alignment'], row['germline_alignment'])
                },
                'chain': row.get('chain', 'IGH'),
            }
            collection.insert_one(doc)
            row_processed += 1
            if row_processed % 1000 == 0:
                print(f"Processed {row_processed} rows...")
            
            if row_processed >= 70000:
                break

    print(f"Finished processing. Total rows processed: {row_processed}, skipped rows: {skipped_rows}")
    
def main():
    # take command line arguments
    parser = argparse.ArgumentParser(description='Load AIRR data into MongoDB.')
    parser.add_argument('--tsv_folder', type=str, help='Path to the AIRR TSV file')
    parser.add_argument('--serverip', type=str, required=True, help='MongoDB server IP address')
    parser.add_argument('--port', type=int, default=27017, help='MongoDB server port') # default is 27017
    parser.add_argument('--db', type=str, default='clonify_db', help='MongoDB database name')
    # parser.add_argument('--collection', type=str, default='antibodies_7k', help='MongoDB collection name')

    print("Parsing arguments...")
    args = parser.parse_args()
    tsv_folder = args.tsv_folder
    serverip = args.serverip 
    port = args.port
    db_name = args.db
    # collection_name = args.collection

    print(f"Loading AIRR data from {tsv_folder} into MongoDB...")
    # Check if the TSV directory exists
    if not os.path.exists(tsv_folder):
        print(f"Directory {tsv_folder} does not exist.")
        return
    # Get all tsv files from tsv_folder
    tsv_files = [f for f in os.listdir(tsv_folder) if f.endswith('.tsv') or f.endswith('.csv')] 
    if not tsv_files:
        print(f"No TSV files found in {tsv_folder}.")
        return
    # Load each TSV file into MongoDB
    for tsv_file in tsv_files:
        collection_name = tsv_file.split('.')[0]
        client, collection = config_db(serverip, port, db_name, collection_name)
        if collection is None:
            print(f"Collection {collection_name} already exists. Skipping...")
            continue
        print(f"Connected to MongoDB at {serverip}:{port}, database: {db_name}, collection: {collection_name}")
        tsv_path = os.path.join(tsv_folder, tsv_file)
        load_airr_data(tsv_path, collection)
        print(f"Loaded {tsv_file} into MongoDB collection {collection_name}")
    
    # Close the database connection
    client.close()



if __name__ == "__main__":
    main()
