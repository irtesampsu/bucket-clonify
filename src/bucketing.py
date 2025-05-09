from datasketch import MinHash, MinHashLSH
from collections import defaultdict
from clonify3 import vprint
import hashlib

def get_kmers(seq, k=5):
    return set(seq[i:i+k] for i in range(len(seq) - k + 1))

def compute_minhash(kmers, num_perm=128):
    m = MinHash(num_perm=num_perm)
    for kmer in kmers:
        m.update(kmer.encode('utf-8'))
    return m

def build_lsh_index(seqs, k=7, num_perm=128, threshold=0.8):
    lsh = MinHashLSH(threshold=threshold, num_perm=num_perm)
    mh_table = {}
    for seq in seqs:
        #TODO: What happens if we compute minhash on V-J or sth
        kmers = get_kmers(seq.junc, k)
        m = compute_minhash(kmers, num_perm)
        
        lsh.insert(seq, m) 
        mh_table[seq.id] = m
    return lsh, mh_table

def assign_to_best_bucket(seqs, lsh, mh_table):
    # group seqs by their best bucket (max intersection count)
    bucket_map = {}
    bucket_key = 0
    
    already_considered_seqs = set()
    single_entry_list = []
    for seq in seqs:
        # For now, if some seq has already been into some bucket, we do not consider it again
        if seq.id in already_considered_seqs:
            continue
        else:
            bucket_key += 1

        m = mh_table[seq.id]
        candidates = lsh.query(m)
        if len(candidates) <= 1:
            #TODO: Discuss: We can skip this case without adding
            # But we are getting rid of sequences which might not be a good for cluster comparison
            single_entry_list.append(candidates[0])
            already_considered_seqs.add(candidates[0].id)
            continue
        # we can fix the threshold when we build the index, so no filtering here
        selected_cand_list = []
        candidates_id_list = []
        for candidate in candidates:
            if candidate.id not in already_considered_seqs:
                selected_cand_list.append(candidate)
                candidates_id_list.append(candidate.id)
        
        already_considered_seqs.update(candidates_id_list)  # Use update instead of +=
        # print(candidates_id_list)
        # print(len(candidates))
        # print(candidates)
        if len(selected_cand_list) > 1:
            bucket_map[bucket_key] = selected_cand_list
        else:
            single_entry_list.append(selected_cand_list[0])
        # print(len(bucket_map[bucket_key]))
    if(len(single_entry_list) > 1):
        bucket_map[bucket_key+1] = single_entry_list
    elif(len(single_entry_list) == 1):
        # only one element, instead of skipping it, adding it to the first bucket
        bucket_map[1].append(single_entry_list[0])
    bucket_list = list(bucket_map.values())
    threshold = 10
    final_buckets = []
    seperate_bucket = []
    for bucket in bucket_list:
        if len(bucket) < threshold:
            seperate_bucket += bucket
        else:
            final_buckets.append(bucket)
    if len(seperate_bucket) < threshold:
        if len(final_buckets) > 0:
            final_buckets[-1] += seperate_bucket
        else:
            final_buckets.append(seperate_bucket)
    else:
        final_buckets.append(seperate_bucket)

    bucket_mem_count_list = [len(bucket) for bucket in final_buckets]
    vprint('Number of elements inside bucket: ', bucket_mem_count_list)
    return final_buckets
