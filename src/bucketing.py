from datasketch import MinHash, MinHashLSH
from collections import defaultdict
import hashlib

def get_kmers(seq, k=5):
    return set(seq[i:i+k] for i in range(len(seq) - k + 1))

def compute_minhash(kmers, num_perm=128):
    m = MinHash(num_perm=num_perm)
    for kmer in kmers:
        m.update(kmer.encode('utf-8'))
    return m

def build_lsh_index(seqs, k=5, num_perm=128, threshold=0.8):
    lsh = MinHashLSH(threshold=threshold, num_perm=num_perm)
    mh_table = {}
    for seq in seqs:
        kmers = get_kmers(seq.junc, k)
        m = compute_minhash(kmers, num_perm)
        key = hashlib.md5(seq.id.encode()).hexdigest()
        lsh.insert(key, m)
        mh_table[seq.id] = (key, m, seq)
    return lsh, mh_table

def assign_to_best_bucket(seqs, lsh, mh_table):
    # group seqs by their best bucket (max intersection count)
    bucket_map = defaultdict(list)
    for seq in seqs:
        key, m, _ = mh_table[seq.id]
        candidates = lsh.query(m)
        if not candidates:
            bucket_map[seq.id].append(seq)  # singleton bucket
            continue
        best = max(candidates, key=lambda c: m.jaccard(mh_table[c][1]))
        bucket_map[best].append(seq)
    return list(bucket_map.values())
