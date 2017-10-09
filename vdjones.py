import tqdm
import os
import logging
import sys
import copy

from collections import defaultdict
import cPickle as pkl
import numpy as np

import parallel as par
import adj_list as aj
import output as out

def strs_to_int8(l):
    return [np.fromstring(x.strip(), np.int8) for x in l]


def gauss_interval(n, k):
    """
    Given n tasks, where task i takes time ~ (n - i), as with computing
    all pairwise distances between a set of points, determines the fairest
    splits to divide work evenly in k continuous ranges.

    Args:
        n: total number of tasks
        k: number of desired chunks

    Returns:
        list of k (start, end) lists, each of which corresponds to a 
        "fair" chunk of the total workload.
    """
    pairs = []
    j = 0
    for i in range(k-1):
        pair = [n-j,n-j]
        roots = np.roots([1., 1., -(j**2 + j) - (n**2 + n)/k])
        j = int(np.max(roots))
        pair[0] = n-j
        pairs.append(pair)
    pairs.append([0, n-j])
    return pairs[::-1]


def load_data(FASTA_FILE, V_SEG_LENGTH=300):
    """
    Import CDRH3 and V segment DNA from FASTA file
    """

    cdr3s = []
    all_v_segs = []
    v_seg_counts = defaultdict(int)

    all_lines = None
    with open(FASTA_FILE, 'r') as f:
        all_lines = f.readlines()

    for line in tqdm.tqdm(all_lines):
        if line[0] == '>':
            s = line.split(';')
            accession = s[0].split(' ')[0][1:]
            cdr3_dna = s[12]
            v_dna = s[13][-V_SEG_LENGTH:]
            
            # Filter out all clones whose V segments aren't at least V_SEG_LENGTH nucleotides long             
            if len(v_dna) == V_SEG_LENGTH:
                cdr3s.append((accession, cdr3_dna))
                all_v_segs.append((accession, v_dna))
                v_seg_counts[v_dna] += 1
                    
    s = set()
    for x in all_v_segs:
        s.add(x[1])
        
    uniq_v_segs = list(s)

    # {uniq v seg idx: uniq v seg}
    v_seg_dict = {}
    for idx, v_seg in enumerate(uniq_v_segs):
        v_seg_dict[v_seg] = idx

    # v_segs is a list of the unique v segments, represented in such a way that 
    # is easy to compute connectivity matrix on.
    v_segs_int8 = strs_to_int8(uniq_v_segs)

    return v_seg_counts, v_seg_dict, v_segs_int8, cdr3s, all_v_segs


def gen_adj_list(v_segs, start, stop):
    """
    Generate adjacency list for V segments in v_segs
    """
    num_segs = len(v_segs)
    logging.info("Starting!")
    adj_list = defaultdict(set)
    one_percent = (stop - start) / 100
    for i in range(start, stop):
        if ((i - start) % one_percent == 0):
            logging.info('Completed {}%'.format(100 * (i-start)/(stop-start)))
        v_seg = v_segs[i]
        for j in range(i+1, num_segs):
            if aj.isNeighbor2(v_seg, v_segs[j], 300):
                adj_list[i].add(j)
                adj_list[j].add(i)
    return adj_list


def merge_adj_lists(adj_lists):
    """
    Merge all the dicts resulting from the adjacency lists.
    """

    merged_list = defaultdict(set)
    for i, l in enumerate(adj_lists):
        for key in l:
            merged_list[key] |= l[key]
    return merged_list


def pick_from_cluster(clusters, all_v_segs, merged_list, v_seg_dict):
    cluster_reps = {}
    for cluster in clusters:
        most_connected = None
        most_connected_degree = -1
        # seq is an index into the original FASTA file
        for seq in clusters[cluster]:
            acc, v_seg = all_v_segs[seq]
            if len(v_seg) == 300:
                neighbors = merged_list[v_seg_dict[v_seg]]
                if len(neighbors) > most_connected_degree:
                    most_connected = v_seg
                    most_connected_degree = len(neighbors)
        
        if most_connected is not None:
            cluster_reps[cluster] = (most_connected, most_connected_degree)

    return cluster_reps


def find_dominant_clones(v_seg_counts, connectivity_order, merged_cluster_adj_list,
                         non_singletons_v_segs):
    # Find dominant clones by collapsing low connectivity, low frequency nodes into their more connected neighbors.

    surviving_nodes = {}

    v_seg_scores = copy.copy(v_seg_counts)

    for v in connectivity_order:
        survive = True
        max_nb = None
        max_nb_conn = -1
        
        my_conn = len(merged_cluster_adj_list[v])
        my_freq = v_seg_scores[non_singletons_v_segs[v]]
        
        for nb in merged_cluster_adj_list[v]:
            nb_conn = len(merged_cluster_adj_list[nb])
            nb_freq = v_seg_scores[non_singletons_v_segs[nb]]
            
            if nb_conn > my_conn and nb_freq >= my_freq and nb_conn > max_nb_conn:
                max_nb = nb
                max_nb_conn = nb_conn
                survive = False
                v_seg_scores[non_singletons_v_segs[nb]] += my_freq
                break
                
        if survive:
            surviving_nodes[v] = v_seg_scores[non_singletons_v_segs[v]]

    return surviving_nodes    


if __name__ == "__main__":

    PROJECT_ROOT = '/Users/rishi/dbio/mnt'
    FASTA_FILE = os.path.join(PROJECT_ROOT, 'humans-separated/cluster/091/091-annot.fa')

    logging.basicConfig(stream=sys.stdout,
                        format='%(asctime)s %(levelname)s %(process)d: ' +
                        '%(message)s',
                        level=logging.INFO)

    ## Load Data
    logging.info("Loading data")
    v_seg_counts, v_seg_dict, v_segs_int8, cdr3s, all_v_segs = load_data(FASTA_FILE)
    logging.info("Data loaded")

    ## Connectivity
    logging.info("Starting adjacency list generation")
    num_segs = len(v_segs_int8)
    intervals = gauss_interval(num_segs, 6)
    inputs = [(v_segs_int8, pair[0], pair[1]) for pair in intervals]
    adj_lists = par.submit_jobs(gen_adj_list, inputs, 6)
    logging.info("Merging adjacency lists")
    merged_list = merge_adj_lists(adj_lists)

    ## CDR3 Clustering
    logging.info("Clustering sequences by CDR3")
    clusters = defaultdict(list)
    for idx, cdr3 in tqdm.tqdm(enumerate(cdr3s)):
        clusters[cdr3[1]].append(idx)


    ## Pick from Clusters
    logging.info("Picking from clusters")
    cluster_reps = pick_from_cluster(clusters, all_v_segs, merged_list, v_seg_dict)
    cluster_rep_v_segs = [cluster_reps[cdr3][0] for cdr3 in cluster_reps]
    logging.info("Rendering non-redundant and non-singleton")
    cluster_rep_v_segs = list(set(cluster_rep_v_segs))
    non_singletons_v_segs = [x for x in cluster_rep_v_segs if len(merged_list[v_seg_dict[x]]) > 0]

    ## Connectivity of Cluster Reps
    logging.info("Starting adjacency list of cluster reps generation")
    cluster_v_segs_int8 = [np.fromstring(v_seg.strip(), np.int8) for v_seg in non_singletons_v_segs]
    num_cluster_segs = len(cluster_v_segs_int8)
    intervals = gauss_interval(num_cluster_segs, 4)
    cluster_inputs = [(cluster_v_segs_int8, pair[0], pair[1]) for pair in intervals]
    cluster_adj_lists = par.submit_jobs(gen_adj_list, cluster_inputs, 4)
    logging.info("Merging adjacency lists")
    merged_cluster_adj_list = merge_adj_lists(cluster_adj_lists)
    
    ## Find dominant clones
    # Sorted the cluster reps in increasing order of connectivity and then frequency
    logging.info("Finding dominant clones")
    connectivity_order = sorted(merged_cluster_adj_list.keys(), 
                                key=lambda x: (len(merged_cluster_adj_list[x]), 
                                               v_seg_counts[non_singletons_v_segs[x]]))
    surviving_nodes = find_dominant_clones(v_seg_counts, connectivity_order, 
                                           merged_cluster_adj_list,
                                           non_singletons_v_segs)

    ## Output them to a FASTA file
    logging.info("Outputting dominant clones as FASTA")
    sorted_survivors = sorted(surviving_nodes.keys(), 
                              key=lambda x: surviving_nodes[x], reverse=True)

    sorted_survivors_v_segs = [non_singletons_v_segs[v] for v in sorted_survivors]
    sorted_survivors_scores = [surviving_nodes[v] for v in sorted_survivors]

    out.dump_fasta(sorted_survivors_v_segs, 'dominant.clones', sorted_survivors_scores)
