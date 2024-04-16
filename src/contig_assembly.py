#!/usr/bin/env python

import sys
import argparse
import os
#from .contigs import *
import contigs
import alignment

## Parse the input arguments
parser=argparse.ArgumentParser()
parser.add_argument('-s','--SequencerReads',type=str,help='Path to sequencer reads file which is in FASTA format')
parser.add_argument('-q','--QuerySequence',type=str,help='Path to query sequence which is in FASTA format')
parser.add_argument('-k','--kmerLen',type=int,help='k-mer size')
parser.add_argument('-sl','--seedLen',type=int,help='seed length for creating hash index for alignment')
parser.add_argument('-o','--outputPath',type=str, help='Path for output .fasta and .aln files')

args=parser.parse_args()

## Generate dictionaries
id_kmer_dict,id_read_dict,id_presuffix_dict,nodes_unique,edges_all,suffixes,kmer_id_coord=contigs.generate_kmer_nodes(seq_path=args.SequencerReads,k_len=args.kmerLen)
edges_hash=contigs.build_debruijn_graph(id_presuffix_dict)
nodes_degree,starting_nodes,ending_nodes=contigs.calculate_node_degree(id_presuffix_dict)

## de Bruijn graph
edges_hash_pop=edges_hash.copy()
contig_list={}
contig_len=0
longest_contig=''
for i in starting_nodes:
    ep=contigs.eulerian_path_(edges_hash_pop,i)
    contig=contigs.concat_contig(ep,k_len=args.kmerLen)
    contig_list['starting_node_'+i]=contig
    if len(contig) > contig_len:
        contig_len=len(contig)
        longest_contig=contig

## Depth-first search
edges_hash=contigs.build_debruijn_graph(id_presuffix_dict)
dfs_contig_list={}
dfs_contig_len=0
dfs_longest_contig=''
for i in starting_nodes:
    e_paths=contigs.depthfirst_eulerian_path_(edges_hash,start=i,k_len=args.kmerLen)
    for dfs_idx,dfs_contig in enumerate(e_paths):
        dfs_contig_list['contig_start_node_'+i+'_path_'+str(dfs_idx)]=dfs_contig
        if len(dfs_contig) > dfs_contig_len:
            dfs_contig_len=len(dfs_contig)
            dfs_longest_contig=dfs_contig

## Create hash index
query_seq,query_to_kmer_index=alignment.fetch_query_index(args.QuerySequence,seed_len=args.seedLen)

## query contig seeds in the hash index database
sorted_align_query,longest_contig = alignment.seed_extend_align(query_seq,seed_len=args.seedLen,dfs_contig_list=dfs_contig_list,query_to_kmer_index=query_to_kmer_index)

## create contig kmers 
kmer_contig_dict=alignment.kmer_contig(k_len=args.seedLen,longest_contig=longest_contig,kmer_id_coord=kmer_id_coord)

## Identify coordinate overlaps between reads and the longest contig
aln_list=alignment.aln_file(kmer_contig_dict=kmer_contig_dict,sorted_align_query=sorted_align_query)

output_path_fasta=os.path.join(args.outputPath,"./ALLELES.fasta")
with open(output_path_fasta, 'w') as output_file:
    header_line=">longest_contig\n"
    output_file.write(header_line)
    contig_line=longest_contig+"\n"
    output_file.write(contig_line)

output_path_aln=os.path.join(args.outputPath,"./ALLELES.aln")
with open(output_path_aln, 'w') as output_file:
    header_line="sseqid\tqseqid\tsstart\tsend\tqstart\tqend\n"
    output_file.write(header_line)
    for contig_line in aln_list:
        contig_line_write = '\t'.join(map(str, contig_line)) + '\n'
        output_file.write(contig_line_write)

output_path_contigs=os.path.join(args.outputPath,"./CONTIGS.fasta")
with open(output_path_contigs, 'w') as output_file:
    #header_line="sseqid\tqseqid\tsstart\tsend\tqstart\tqend\n"
    for k,v in dfs_contig_list.items():
        header_line=">"+str(k)+"\n"
        output_file.write(header_line)
        contig_line=str(v)+"\n"
        output_file.write(contig_line)