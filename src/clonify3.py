import os
import time
import math
import argparse
import itertools
from multiprocessing import Pool, cpu_count

from pymongo import MongoClient

import numpy as np

import fastcluster as fc
from scipy.cluster.hierarchy import fcluster

from Bio import pairwise2
from Levenshtein import distance

from itertools import zip_longest
import time
from functools import wraps
from contextlib import contextmanager
from faiss_bucketing import faiss_bucketing


@contextmanager
def time_block(name="Block"):
    start = time.time()
    yield
    end = time.time()
    print(f"[TIMER] {name} executed in {end - start:.4f} seconds.")

def time_function(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        with time_block(name=func.__name__):
            result = func(*args, **kwargs)
        return result
    return wrapper



parser = argparse.ArgumentParser()
parser.add_argument('-d', '--db', required=True)
parser.add_argument('-c', '--coll', default=None)
parser.add_argument('-i', '--ip', default='localhost')
parser.add_argument('-p', '--port', default=27017, type=int)
parser.add_argument('-o', '--out', dest='output', default='')
parser.add_argument('-s', '--split_by', default='none', choices=['none', 'fam', 'gene'])
parser.add_argument('-t', '--threads', type=int, default=None)
parser.add_argument('-x', '--dist', dest='distance_cutoff', default=0.35, type=float)
parser.add_argument('-z', '--no_split', action='store_true', default=False)
parser.add_argument('-n', '--nt', dest='is_aa', action='store_false', default=True)
parser.add_argument('-u', '--no_update', dest='update', action='store_false', default=True)
parser.add_argument('-k', '--chunksize', type=int, default=500)
parser.add_argument('-b', '--bucket', action='store_true', default=False)
args = parser.parse_args()


class Seq:
    def __init__(self, data, junc_query):
        self.id = data['seq_id']
        v_full = data['v_gene']['full']
        j_full = data['j_gene']['full']
        self.v_fam = v_full.split('-')[0][4:]
        self.v_gene = '-'.join(v_full.split('*')[0].split('-')[1:])
        self.v_all = v_full.split('*')[1]
        self.j_gene = j_full.split('*')[0][4:]
        self.j_all = j_full.split('*')[1]
        self.junc = data[junc_query]
        self.junc_len = len(self.junc)
        assert self.junc_len > 0, f"Junction length is zero for sequence ID {self.id}"
        self.muts = [f"{d['loc']}{d['mut']}" for d in data.get('var_muts_nt', {}).get('muts', [])]

    def v_gene_string(self):
        return f'v{self.v_fam}-{self.v_gene}'

    def v_fam_string(self):
        return f'v{self.v_fam}'


def clonify(pair):
    ichunk, jchunk = pair
    return [get_scores(i, jchunk) for i in ichunk]


def get_scores(i, jchunk):
    results = []
    for j in jchunk:
        if i.id == j.id:
            results.append(0.0)
            continue
        LD = get_LD(i, j)
        vPenalty = vCompare(i, j)
        jPenalty = jCompare(i, j)
        lenPenalty = abs(i.junc_len - j.junc_len) * 2
        editLength = min(i.junc_len, j.junc_len)
        mutBonus = sharedMuts(i, j)
        if mutBonus > (LD + vPenalty + jPenalty):
            mutBonus = LD + vPenalty + jPenalty - 0.001
        results.append((LD + vPenalty + jPenalty + lenPenalty - mutBonus) / editLength)
    return results


def get_LD(i, j):
    if i.junc_len == j.junc_len:
        identity = pairwise2.align.globalms(i.junc, j.junc, 1, 0, -50, -50, score_only=True, one_alignment_only=True)
        return i.junc_len - identity
    else:
        return distance(i.junc, j.junc)


def vCompare(i, j):
    return 10 if i.v_gene != j.v_gene else 0


def jCompare(i, j):
    return 8 if i.j_gene != j.j_gene else 0


def sharedMuts(i, j):
    if i.id == j.id:
        return 0.0
    return sum(0.35 for mut in i.muts if mut and mut in j.muts)


def get_seqs(database, collection):
    conn = MongoClient(args.ip, args.port)
    db = conn[database]
    c = db[collection]
    junc_query = 'junc_aa' if args.is_aa else 'junc_nt'
    results = c.find({'chain': 'IGH'}, {'_id': 0, 'seq_id': 1, 'v_gene': 1, 'j_gene': 1, junc_query: 1, 'var_muts_nt': 1})
    return [Seq(r, junc_query) for r in results if r and junc_query in r]


def get_collections():
    if args.coll:
        return [args.coll]
    conn = MongoClient(args.ip, args.port)
    db = conn[args.db]
    return sorted(db.list_collection_names())


def update_db(database, collection, clusters):
    conn = MongoClient(args.ip, args.port)
    db = conn[database]
    c = db[collection]
    clust_sizes = []
    clust_count = 0
    for clust_id, seqs in clusters.items():
        clust_size = len(seqs)
        seq_ids = [s.id for s in seqs]
        if clust_size > 1:
            if args.update:
                c.update_many({'seq_id': {'$in': seq_ids}}, {'$set': {'clonify': {'id': clust_id, 'size': clust_size}}})
            clust_count += 1
            clust_sizes.append(clust_size)
    clustered_seqs = sum(clust_sizes)
    avg_clust_size = clustered_seqs / clust_count if clust_count else 0
    max_clust_size = max(clust_sizes) if clust_sizes else 0
    return [clust_count, clustered_seqs, avg_clust_size, max_clust_size]


def count_cpus():
    return args.threads if args.threads else cpu_count()


def split_by_gene(seqs):
    split = {}
    for seq in seqs:
        split.setdefault(seq.v_gene_string(), []).append(seq)
    return split


def split_by_fam(seqs):
    split = {}
    for seq in seqs:
        split.setdefault(seq.v_fam_string(), []).append(seq)
    return split


def get_chunksize(input):
    if args.no_split or len(input) < args.chunksize:
        return len(input)
    s = float(len(input)) / cpu_count()
    return int(s) if int(math.ceil(s)) * (cpu_count() - 1) > len(input) else int(math.ceil(s))


def chunk_maker(n, iterable, fillvalue=None):
    args = [iter(iterable)] * n
    return [[e for e in t if e is not None] for t in zip_longest(*args, fillvalue=fillvalue)]


def grouper_nofill(n, iterable):
    it = iter(iterable)
    return iter(lambda: list(itertools.islice(it, n)), [])


def build_cluster_dict(count, vh):
    return {f"lineage_{vh}_{c}": [] for c in range(1, count)}


@time_function
def build_matrix(ichunks, chunksize, size, chunk_count):
    matrix = np.zeros((size, size))
    print(f'number of processes: {chunk_count}')
    print(f'matrix: {matrix.shape}')
    print(f'total calculations: {matrix.size}')
    p = Pool(processes=chunk_count)
    for x, seq in enumerate(grouper_nofill(chunk_count, ichunks)):
        result = p.imap(clonify, seq)
        for i, r in enumerate(result):
            matrix[x * chunksize:x * chunksize + len(r), i * chunksize:i * chunksize + len(r[0])] = r
    p.close()
    p.join()
    return matrix

@time_function
def squareform(matrix):
    cm = []
    for i, row in enumerate(matrix[:-1]):
        cm.extend(row[i + 1:])
    return cm


def make_clusters(input_seqs, vh):
    chunksize = get_chunksize(input_seqs)
    print(f'Chunksize is: {chunksize}')
    chunks = chunk_maker(chunksize, input_seqs)
    iter_chunks = itertools.product(chunks, repeat=2)
    distMatrix = build_matrix(iter_chunks, chunksize, len(input_seqs), len(chunks))
    print('condensing the distance matrix...')
    con_distMatrix = squareform(distMatrix)
    print('clustering...')

    with time_block(name="Clustering"):
        linkageMatrix = fc.linkage(np.asarray(con_distMatrix), method='average', preserve_input=False)
        flatCluster = fcluster(linkageMatrix, args.distance_cutoff, criterion='distance')
    
    print('building cluster dict...')
    clusters = build_cluster_dict(max(flatCluster) + 1, vh)
    print('assigning sequences to clusters...')
    for s, cl in enumerate(flatCluster):
        clusters[f'lineage_{vh}_{cl}'].append(input_seqs[s])
    return clusters


def analyze_collection(coll):
    startTime = time.time()
    print(f'\n\n========================================\nprocessing collection: {coll}\n========================================\n')
    print(f'Querying MongoDB (db: {args.db}, collection: {coll}) for the input sequences...')
    seqs = get_seqs(args.db, coll)
    print(f'...done. Retrieved {len(seqs)} sequences.\n')
    if len(seqs) < 2:
        return False
    split_seqs = {'v0': seqs}
    if args.split_by == 'gene':
        split_seqs = split_by_gene(seqs)
    elif args.split_by == 'fam':
        split_seqs = split_by_fam(seqs)
    

    print('Sorting sequences into clonal families...')
    clusters = {}
    for vh in sorted(split_seqs.keys()):
        if len(split_seqs[vh]) <= 1:
            continue
        print(f'\n--------\n{vh}\n--------')

        if args.bucket:
            buckets = faiss_bucketing(split_seqs[vh], k=5)
            bucket_id = 1
            for bucket in buckets:
                if len(bucket) > 1:
                    clusters.update(make_clusters(bucket, vh + str(bucket_id)))
                bucket_id += 1
        else:
            clusters.update(make_clusters(split_seqs[vh], vh))
        
    print('...done.\n')
    if args.update:
        print('Updating MongoDB...')
    else:
        print('Calculating cluster statistics...')
    stats = update_db(args.db, coll, clusters)
    print('...done.\n')
    if args.output:
        print('Writing clonal families to file...')
        write_output(args.output, coll, clusters)
        print('...done.\n')
    else:
        print('No output directory was provided. Lineage assignments are not being written to file.\n')
    print(f'Querying MongoDB took {round(time.time() - startTime, 2)} seconds.')
    print(f"{len(seqs)} sequences were segregated into {stats[0]} clonal families.")
    print(f'The average cluster size was {stats[2]:.2f}.')
    print(f'The largest cluster contains {stats[3]} sequences.')
    print(f'{stats[1]} sequences were assigned to clonal families ({100.0 * stats[1] / len(seqs):.2f}%).\n\n')


def write_output(out_dir, collection, data):
    out_file = os.path.join(out_dir, collection + '_clones.txt')
    with open(out_file, 'w') as f:
        for c, seqs in data.items():
            if len(seqs) < 2:
                continue
            f.write(f'#{c}\n')
            for seq in seqs:
                f.write(f'>{seq.id}\n{seq.junc}\n')
            f.write('\n')


def main():
    total_start_time = time.time()
    for c in get_collections():
        analyze_collection(c)
    print(f'Total execution time: {time.time() - total_start_time:.2f} seconds.')

if __name__ == '__main__':
    main()
