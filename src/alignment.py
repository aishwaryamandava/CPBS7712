
#!/usr/bin/env python
import sys
import argparse
import os

## fetch the query sequence indices and non-overlapping seeds
def fetch_query_index(query_path,seed_len):
    f=open(query_path)
    for line in f:
        if line[0] != '>':
            query_seq=line.strip()
    
    query_to_kmer_index = {}
    for i in range(0,len(line.strip())-(seed_len-1),seed_len):
        nonoverlap_kmer=query_seq[i:i+seed_len]
        if nonoverlap_kmer not in query_to_kmer_index:
            query_to_kmer_index[nonoverlap_kmer]=[i]
        else:
            query_to_kmer_index[nonoverlap_kmer].append(i)
    
    return query_seq, query_to_kmer_index

## Calculate edit distance between query sequence and contig
def edit_distance(query_seq,contig):

    edist_matrix = [[0 for m in range(len(query_seq)+1)]for n in range(len(contig)+1)]
    
    for j in range(1,len(query_seq)+1):
        for i in range(1,len(contig)+1):
            if query_seq[j-1] == contig[i-1]:
                edist_matrix[i][j] = edist_matrix[i-1][j-1]
            else:
                edist_matrix[i][j] = min(edist_matrix[i-1][j-1],edist_matrix[i][j-1],edist_matrix[i-1][j])+1

    return edist_matrix[len(contig)][len(query_seq)]

## First calculate look for seed alignments and then calculate edit distance for filtered contigs
def seed_extend_align(query_seq,seed_len,dfs_contig_list,query_to_kmer_index):   
    align_query={key:{'index':[],'edist':0,'contig_len':0,'edist/contig_len':0} for key in dfs_contig_list}
    for contig_name,contig_seq in dfs_contig_list.items():
        align_query[contig_name]['contig_len']=len(contig_seq)
        for base in range(0,len(contig_seq)-(seed_len-1)):
            if contig_seq[base:base+seed_len] in query_to_kmer_index:
                align_query[contig_name]['index'].append(query_to_kmer_index[contig_seq[base:base+seed_len]])

        if len(align_query[contig_name]['index'])>1 and all(align_query[contig_name]['index'][i] <= align_query[contig_name]['index'][i + 1] for i in range(len(align_query[contig_name]['index']) - 1)):
            align_query[contig_name]['edist']=edit_distance(query_seq,dfs_contig_list[contig_name])
            align_query[contig_name]['edist/contig_len']=edit_distance(query_seq,dfs_contig_list[contig_name])/align_query[contig_name]['contig_len']
        
    sorted_align_query = dict(sorted(
        ((key, value) for key, value in align_query.items() if 'index' in value and len(value.get('index', [])) > 1 and all(value.get('index', [])[i] <= value.get('index', [])[i+1] for i in range(len(value.get('index', [])) - 1))),
        #key=lambda x: (x[1]['edist/contig_len'])
        key=lambda x: (-x[1]['contig_len'],x[1]['edist'])
    ))

    longest_contig=dfs_contig_list[list(sorted_align_query.keys())[0]]

    return sorted_align_query,longest_contig

## Create a dictionary of the contig kmers and indices
def kmer_contig(k_len, longest_contig, kmer_id_coord):
    kmer_contig_dict = {}
    for i in range(0, len(longest_contig) - (k_len - 1)):
        kmer = longest_contig[i:i + k_len]
        if kmer in kmer_id_coord:
            if kmer not in kmer_contig_dict:
                kmer_contig_dict[kmer] = {}
            for id, coords in kmer_id_coord[kmer].items():
                kmer_contig_dict[kmer][id] = {'seq_coords': [(coords[0][0], coords[0][1])]}
                kmer_contig_dict[kmer][id]['contig_coords'] = [(i, i + k_len - 1)]

    return kmer_contig_dict

## Identify coordinate overlaps between reads and the longest contig
def aln_file(kmer_contig_dict,sorted_align_query):
    unique_ids=set()
    longest_contig_dict={}
    aln_list=[]

    for inner_dict in kmer_contig_dict.values():
        for read_id in inner_dict.keys():
            unique_ids.add(read_id)

    for ids in unique_ids:
        seq_coords = []
        contig_coords = []
        for k, v in kmer_contig_dict.items():
            if ids in v:
                seq_coords.extend(v[ids]['seq_coords'])
                contig_coords.extend(v[ids]['contig_coords'])
        longest_contig_dict[ids] = {'seq_coords': seq_coords,
                                'contig_coords': contig_coords}

        s_start = longest_contig_dict[ids]['seq_coords'][0][0]
        s_end = longest_contig_dict[ids]['seq_coords'][0][1]
        q_start = longest_contig_dict[ids]['contig_coords'][0][0]
        q_end = longest_contig_dict[ids]['contig_coords'][0][1]

        prev_s_start = s_start
        prev_q_start = q_start

        for seq_coord, contig_coord in zip(longest_contig_dict[ids]['seq_coords'], longest_contig_dict[ids]['contig_coords']):
            coord1, coord2 = seq_coord
            qcoord1, qcoord2 = contig_coord
            if coord1 - prev_s_start < 2 and qcoord1 - prev_q_start < 2:
                if coord1 < s_start:
                    s_start = coord1
                if coord2 > s_end:
                    s_end = coord2
                if qcoord1 < q_start:
                    q_start = qcoord1
                if qcoord2 > q_end:
                    q_end = qcoord2
                prev_s_start = coord1
                prev_q_start = qcoord1

        aln_list.append((ids[1:],list(sorted_align_query.keys())[0],s_start,s_end,q_start,q_end))

    return aln_list
        