#!/usr/bin/env python
import sys
import argparse
import os

#### Generate k-mers and (k-1)mer prefix and suffix
def generate_kmer_nodes(seq_path,k_len):

    id_kmer_dict = {}
    id_read_dict = {}
    id_presuffix_dict={}
    kmer_id_coord = {}
    nodes_unique=set()
    edges_all = []
    suffixes=set()
    starting_kmers = []

    f=open(seq_path)
    for line in f:
        if line[0] == '>':
            identifier = line.strip()
            id_kmer_dict[identifier]=[]
            id_read_dict[identifier] = []
            id_presuffix_dict[identifier]=[]
        else:
            id_read_dict[identifier] = line.strip()
            for i in range(0,len(line.strip())-(k_len-1)):
                if i==0:
                    starting_kmers.append(line.strip()[i:i+k_len-1])
                id_kmer_dict[identifier].append(line.strip()[i:i+k_len])
                id_presuffix_dict[identifier].append([line.strip()[i:i+k_len-1],line.strip()[i+1:i+k_len]])
                nodes_unique.add(line.strip()[i:i+k_len-1])
                nodes_unique.add(line.strip()[i+1:i+k_len])
                edges_all.append((line.strip()[i:i+k_len-1],line.strip()[i+1:i+k_len]))
                suffixes.add(line.strip()[i+1:i+k_len])
                if line.strip()[i:i+k_len] not in kmer_id_coord:
                    kmer_id_coord[line.strip()[i:i+k_len]]={}
                    kmer_id_coord[line.strip()[i:i+k_len]][identifier]=[(i,i+k_len-1)]
                else:
                    if identifier not in kmer_id_coord[line.strip()[i:i+k_len]]:
                       kmer_id_coord[line.strip()[i:i+k_len]][identifier]=[(i,i+k_len-1)]
                    else:
                       kmer_id_coord[line.strip()[i:i+k_len]][identifier].append((i,i+k_len-1))   
    edges_all_unique = set(edges_all)
    return id_kmer_dict,id_read_dict,id_presuffix_dict,nodes_unique,edges_all,suffixes,kmer_id_coord


#### de Bruijn graph represented as an adjacency list
def build_debruijn_graph(id_presuffix_dict):
    edges_hash={}
    for id, list_of_edges in id_presuffix_dict.items():
        for sub_list in list_of_edges:
            if sub_list[0] in edges_hash.keys():
                edges_hash[sub_list[0]].append(sub_list[1])
            else:
                edges_hash[sub_list[0]]=[sub_list[1]]
    return edges_hash

#### Calculate the in-degree and out-degree for each node in the de Bruijn graph
#### This helps in selecting the start and end nodes
def calculate_node_degree(id_presuffix_dict):
    nodes_degree={}
    for id, list_of_edges in id_presuffix_dict.items():
        for sub_list in list_of_edges:
            if sub_list[0] not in nodes_degree:
                nodes_degree[sub_list[0]]=[0,1]
            else:
                nodes_degree[sub_list[0]][1]+=1
            if sub_list[1] not in nodes_degree:
                nodes_degree[sub_list[1]]=[1,0]
            else:
                nodes_degree[sub_list[1]][0]+=1
    
    starting_nodes = [key for key, value in nodes_degree.items() if (value[1]-value[0]) == 1]
    ending_nodes = [key for key, value in nodes_degree.items() if (value[0]-value[1]) == 1]

    return nodes_degree,starting_nodes,ending_nodes

#### Approach1: traversing the graph. This is for the eulerian path. This only gives one eulerian path for each start node
def eulerian_path_(edges_hash_pop,start_node):
    eulerian_path=[]
    eulerian_path.append(start_node)
    next_node=start_node

    while True:
        try:
            if next_node not in edges_hash_pop:
                break

            next_node=edges_hash_pop[next_node].pop(0)

            if len(edges_hash_pop[next_node]) == 0:
                edges_hash_pop.pop(next_node)

            eulerian_path.append(next_node)

        except (KeyError,IndexError):
            break

    return eulerian_path 

#### Approach2: traversing the graph using depth first search. This is for the eulerian path.
def depthfirst_eulerian_path_(graph,start):
    paths = []
    stack = {start: [start]}  # Initialize stack as a dictionary
    #graph=edges_hash

    while stack:
        node, path = stack.popitem()
        if node not in graph or not graph[node]:
            paths.append(concat_contig(path))
        else:
            for next in graph[node]:
                if next not in path:
                    stack[next] = path + [next]
    return(paths)

#### Given a list of nodes that are visited, this function returns a contiguous sequence
def concat_contig(visited):
    contig=""
    for i,ep in enumerate(visited):
        if i == 0:
            contig+=ep
        else:
            contig+=ep[-1] 
    return contig
