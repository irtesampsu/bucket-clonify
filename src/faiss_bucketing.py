import faiss
import numpy as np
from collections import defaultdict
from sklearn.feature_extraction.text import CountVectorizer

def get_kmers(seq, k=5):
    return ' '.join(seq[i:i+k] for i in range(len(seq) - k + 1))

def vectorize_kmers(seqs, k=5):
    corpus = [get_kmers(seq.junc, k) for seq in seqs]
    vectorizer = CountVectorizer(analyzer='word')
    X = vectorizer.fit_transform(corpus).astype(np.float32)
    X = X.toarray()
    faiss.normalize_L2(X)  # normalize for cosine similarity
    return X, vectorizer

def build_faiss_index(X, nlist=50, nprobe=10):
    d = X.shape[1]
    if X.shape[0] < 1000:
        index = faiss.IndexFlatIP(d)  # brute-force with cosine
        index.add(X)
        return index
    else:
        nlist = min(nlist, max(1, X.shape[0] // 40))
        quantizer = faiss.IndexFlatIP(d)
        index = faiss.IndexIVFFlat(quantizer, d, nlist, faiss.METRIC_INNER_PRODUCT)
        index.train(X)
        index.add(X)
        index.nprobe = min(nprobe, nlist)
        return index


def assign_to_buckets(index, X, seqs, top_k=1):
    D, I = index.search(X, top_k)
    buckets = defaultdict(list)
    for i, indices in enumerate(I):
        best = indices[0]
        buckets[best].append(seqs[i])
    return list(buckets.values())


def faiss_bucketing(seqs, k = 5):
    X, _ = vectorize_kmers(seqs, k)
    N = len(seqs)
    nlist = min(int(np.sqrt(N)), max(1, N // 40))
    nprobe = max(1, nlist // 10)

    index = build_faiss_index(X, nlist=nlist, nprobe=nprobe)
    buckets = assign_to_buckets(index, X, seqs, top_k=1)
    return buckets