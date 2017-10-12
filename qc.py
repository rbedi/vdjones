import matplotlib
matplotlib.use('agg')

import numpy as np
from matplotlib import pyplot as plt
from collections import defaultdict
import os
import copy
import sklearn.metrics
import tqdm
import logging
import sys
import argparse as ap

import output as out

def compute_quality(PUTATIVE_FILE, IMGT_CLASSIFIED, DEBUG_PATH=None, cutoff=None):
    all_genes = set()

    logging.info("Opening file")
    with open(IMGT_CLASSIFIED, 'r') as inputfile:
        lines = inputfile.readlines()

    logging.info("Parsing lines")
    imgt_counts = defaultdict(int)
    for line in tqdm.tqdm(lines):
        if line[0] == '>': 
            features = line.split(';')
            if len(features) > 1:
                v_assignment = features[1].split(' ')
                if len(v_assignment) == 3:
                    assignment = v_assignment[0]
                    align_len = int(v_assignment[1])
                    muts = int(v_assignment[2])
                    imgt_counts[assignment] += 1
                    if muts < 4 and align_len > 250:
                        all_genes.add(assignment)

    total_count = float(sum(imgt_counts.values()))
    imgt_freqs = {k : imgt_counts[k] / total_count for k in imgt_counts}

    logging.info("All genes present in source file: {}".format(all_genes))

    found_genes = defaultdict(int)
    ROC_data = []
    PRC_data = []
    WROC_data = []
    FP = 0.0
    TP = 0.0
    
    lines = open(PUTATIVE_FILE, 'r').readlines()
    if cutoff is not None:
        lines = lines[:cutoff]

    vseg_to_status = {}

    for idx, line in tqdm.tqdm(enumerate(lines)):
        if line[0] == '>': 
            IMGT_class, matchlen, dist = line.split(';')[1].split(' ')
            if int(dist) < 4:
                found_genes[IMGT_class] += 1
                if found_genes[IMGT_class] == 1:
                    vseg_to_status[lines[idx+1].strip()] = 'GERMLINE'
                    TP += 1
                else:
                    vseg_to_status[lines[idx+1].strip()] = 'ALTHIT'
            else:
                FP += 1
                vseg_to_status[lines[idx+1].strip()] = 'FP'
            
            FN = (len(all_genes) - len(found_genes)) * 1.0
            perc_explained = sum([imgt_freqs[k] for k in found_genes])
            ROC_data.append([TP, FP])
            WROC_data.append([perc_explained, FP])
            PRC_data.append([TP / (TP + FP), TP / (TP + FN)])

    if DEBUG_PATH is not None:
        logging.info('Dumping germline status Cytoscape attrs file')
        vseg_to_id_file = os.path.join(DEBUG_PATH, 'cluster_reps.txt')
        vseg_to_id_dict = out.recover_ordered_list(vseg_to_id_file)

        attrs = {}
        for vseg in vseg_to_id_dict:
            if vseg in vseg_to_status:
                attrs[vseg_to_id_dict[vseg]] = vseg_to_status[vseg]

        germline_attrs_file = os.path.join(DEBUG_PATH, 'for-cytoscape-germline-status.attr')
        out.dump_cyto_attr(attrs, 'IsGermline', germline_attrs_file)

    not_found_genes = all_genes - set(found_genes.keys())
    print "These genes were not found: {}".format(not_found_genes)

    ROC_raw = np.array(ROC_data, dtype='float32')
    PRC_raw = np.array(PRC_data, dtype='float32')
    WROC_raw = np.array(WROC_data, dtype='float32')
    
    ROC_raw[:, 0] /= len(all_genes)
    ROC_raw[:, 1] /= np.max(ROC_raw[:, 1])
    WROC_raw[:, 1] /= np.max(WROC_raw[:, 1])

    return ROC_raw, PRC_raw, WROC_raw

def evalCurve(data, output_path, title='ROC'):
    plt.figure()
    plt.title(title)
    for i, d in enumerate(data):
        raw, name, color = d
        plt.plot(raw[:, 1], raw[:, 0], color=color, label=name)
        AUC = sklearn.metrics.auc(raw[:, 1], raw[:, 0])
        print "{} AUC of {}: {}".format(title, name, AUC)
    plt.legend(loc='lower left')
    plt.savefig(output_path)

if __name__ == "__main__":

    p = ap.ArgumentParser(description='VDJones QC')
    p.add_argument(metavar='/path/to/output', type=str, dest='outpath',
                   help='Path to output')
    p.add_argument(metavar='091 Human', type=str, dest='run_name',
                   help='Title identifying run being QCed')

    args = p.parse_args()

    logging.basicConfig(stream=sys.stdout,
                        format='%(asctime)s %(levelname)s %(process)d: ' +
                        '%(message)s',
                        level=logging.INFO)

    RESULTS_DIR = args.outpath
    FULL_ALGO_GENES = os.path.join(RESULTS_DIR, 'qc/full-alg-v-segs.VDJ.H3.L3.CH1.dnaH3.fa')
    FREQ_ONLY_GENES = os.path.join(RESULTS_DIR, 'qc/most-freq-v-segs.VDJ.H3.L3.CH1.dnaH3.fa')
    IMGT_CLASSIFIED = os.path.join(RESULTS_DIR, 'input_seqs_imgt_annot.fa')
    DEBUG_PATH = RESULTS_DIR

    ROC_raw, PRC_raw, WROC_raw = compute_quality(FULL_ALGO_GENES, IMGT_CLASSIFIED, DEBUG_PATH)#, cutoff=750)
    ROC_baseline, PRC_baseline, WROC_baseline = compute_quality(FREQ_ONLY_GENES, IMGT_CLASSIFIED)#, cutoff=750)

    evalCurve([(ROC_raw, 'Full Alg.', 'green'),
               (ROC_baseline, 'Freq. Only', 'red')], 
               os.path.join(RESULTS_DIR,'qc/roc.png'),
              title='ROC ({})'.format(args.run_name))

    evalCurve([(PRC_raw, 'Full Alg.', 'green'),
               (PRC_baseline, 'Freq. Only', 'red')],
               os.path.join(RESULTS_DIR,'qc/prc.png'),
               title='PRC ({})'.format(args.run_name))

    evalCurve([(WROC_raw, 'Full Alg.', 'green'),
               (WROC_baseline, 'Freq. Only', 'red')],
               os.path.join(RESULTS_DIR,'qc/wroc.png'),
               title='Weighted ROC ({})'.format(args.run_name))
