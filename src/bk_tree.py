# -*- coding: utf-8 -*-
"""
BK-tree bucketing for Clonify
============================
Builds a Burkhard–Keller tree keyed by CDR3 sequences (seq.junc) and
returns disjoint buckets whose members are all within an edit-distance
threshold t of the bucket seed.

Usage (in clonify3.py):
    from bk_bucketing import bucket_by_bktree
    buckets = bucket_by_bktree(split_seqs[vh], t=3)
    for b in buckets:
        if len(b) > 1:
            clusters.update(make_clusters(b, vh))
"""

from Levenshtein import distance as edit_dist
from collections import defaultdict

class BKNode:
    __slots__ = ("seq", "children")
    def __init__(self, seq):
        self.seq = seq          # a Seq object
        self.children = {}      # key: edit distance, value: BKNode

class BKTree:
    def __init__(self, dist=edit_dist):
        self.root = None
        self.dist = dist

    def add(self, seq_obj):
        if self.root is None:
            self.root = BKNode(seq_obj)
            return
        node = self.root
        while True:
            d = self.dist(seq_obj.junc, node.seq.junc)
            if d in node.children:
                node = node.children[d]
            else:
                node.children[d] = BKNode(seq_obj)
                break

    def query(self, seq_str, t):
        """Return all Seq objects within distance ≤ t of seq_str."""
        if self.root is None:
            return []
        result, stack = [], [self.root]
        while stack:
            node = stack.pop()
            d = self.dist(seq_str, node.seq.junc)
            if d <= t:
                result.append(node.seq)
            lo, hi = d - t, d + t
            for child_d, child_node in node.children.items():
                if lo <= child_d <= hi:
                    stack.append(child_node)
        return result

# ----------  public convenience API  ---------- #

def bucket_by_bktree(seqs, t=3):
    """
    Build BK-tree and return list of buckets (each a list[Seq]).
    Greedy pass: pop an unvisited sequence, retrieve its neighbours,
    mark them visited, repeat.
    """
    tree = BKTree()
    for s in seqs:
        tree.add(s)

    visited = set()
    buckets = []
    for s in seqs:
        if s.id in visited:
            continue
        candidates = tree.query(s.junc, t)
        bucket = []
        for x in candidates:
            if x.id not in visited:
                # mark all candidates as visited
                visited.add(x.id)
                bucket.append(x)        
        buckets.append(bucket)
    return buckets
