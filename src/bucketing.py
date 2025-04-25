# kmer_filter.py

from datasketch import MinHash, MinHashLSH
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

def kmerize(sequence, k=7):
    """Generate a set of k-mers from a sequence."""
    return {sequence[i:i+k] for i in range(len(sequence) - k + 1)}

def compute_sketch(seq, k=7, num_perm=128):
    """Compute MinHash for a Seq object."""
    kmers = kmerize(seq.junc, k)
    m = MinHash(num_perm=num_perm)
    for kmer in kmers:
        m.update(kmer.encode("utf8"))
    return seq.id, m

def build_lsh_index(seqs, k=7, threshold=0.85, num_perm=128, threads=4):
    """
    Build an LSH index from Seq objects.
    Returns the LSH index and a dict of {id: MinHash}.
    """
    lsh = MinHashLSH(threshold=threshold, num_perm=num_perm)
    minhashes = {}

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for sid, mh in executor.map(lambda s: compute_sketch(s, k, num_perm), seqs):
            minhashes[sid] = mh
            lsh.insert(sid, mh)

    return lsh, minhashes

def get_candidate_pairs(seqs, lsh, minhashes, id_index):
    """
    Return a dict mapping sequence IDs to candidate similar sequence IDs.
    """
    candidates = defaultdict(list)
    for s in seqs:
        similar_ids = lsh.query(minhashes[s.id])
        candidates[s.id] = [id_index[nid] for nid in similar_ids if nid != s.id]
    return candidates

def index_seqs_by_id(seqs):
    """Returns a dict mapping ID to Seq object."""
    return {s.id: s for s in seqs}
