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

def build_lsh_index(seqs, k=7, num_perm=128, threshold=0.5):
    lsh = MinHashLSH(threshold=threshold, num_perm=num_perm)
    mh_table = {}
    for seq in seqs:
        #TODO: What happens if we compute minhash on V-J or sth
        kmers = get_kmers(seq.junc, k)
        m = compute_minhash(kmers, num_perm)
        
        # key = hashlib.md5(seq.id.encode()).hexdigest()
        lsh.insert(seq, m) #TODO: why need key? why not use seq
        mh_table[seq.id] = m #TODO: do we need to store key and seq?
    return lsh, mh_table

def assign_to_best_bucket(seqs, lsh, mh_table):
    # group seqs by their best bucket (max intersection count)
    # bucket_map = defaultdict(list)
    bucket_map = {}
    bucket_key = 0
    # seq_to_bucket_dic = {}
    already_considered_seqs_list = []
    # seq_length = []
    single_entry_list = []
    for seq in seqs:
        # seq_length.append(len(seq.junc))
        # For now, if some seq has already been into some bucket, we do not cosider it again
        if seq.id in already_considered_seqs_list:
            continue
        else:
            bucket_key += 1

        m = mh_table[seq.id]
        candidates = lsh.query(m)
        if len(candidates) <= 1:
            #TODO: We can skip this case without adding
            # But we are getting rid of sequences which might not be a good for cluster comparison
            # TODO: to discuss
            # bucket_map[seq.id].append(seq)  # singleton bucket
            single_entry_list.append(candidates[0])
            already_considered_seqs_list.append(candidates[0].id)
            continue
        # we can fix the threshold when we build the index, so no filtering here
        # best = max(candidates, key=lambda c: m.jaccard(mh_table[c][1]))
        
        selected_cand_list = []
        candidates_id_list = []
        for candidate in candidates:
            if candidate.id not in already_considered_seqs_list:
                selected_cand_list.append(candidate)
                candidates_id_list.append(candidate.id)
        
        # candidates_id_list = [candidate.id for candidate in candidates]
        already_considered_seqs_list += candidates_id_list
        

        # for candidate in candidates: 
            # seq_to_bucket_dic[candidate] = bucket_key
        # print(candidates_id_list)
        # print(len(candidates))
        # print(candidates)
        if len(selected_cand_list) > 1:
            bucket_map[bucket_key] = selected_cand_list
        else:
            single_entry_list.append(selected_cand_list[0])
        # print(len(bucket_map[bucket_key]))
    # print(f'Num of seqs: {len(seqs)}')
    # print(f'Num of buckets: {len(bucket_map.keys())}')
    # bucket_mem_count_list = [len(bucket_map[key]) for key in bucket_map.keys()]
    # print('Number of elements inside bucket: ', bucket_mem_count_list)
    # print(seq_length)
        # bucket_map[best].append(seq)
    if(len(single_entry_list) > 1):
        # print(f'single entry {len(single_entry_list)}')
        bucket_map[bucket_key+1] = single_entry_list
    elif(len(single_entry_list) == 1):
        # only one element, instead of skipping it, adding it to the first bucket
        # if 1 not in bucket_map.keys():
        #     bucket_map[1] = []
        bucket_map[1].append(single_entry_list[0])
    bucket_mem_count_list = [len(bucket_map[key]) for key in bucket_map.keys()]
    print('Number of elements inside bucket: ', bucket_mem_count_list)
    return list(bucket_map.values())
