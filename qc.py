import numpy as np
from matplotlib import pyplot as plt
from collections import defaultdict
import os
import copy
import sklearn.metrics
import tqdm
import logging
import sys

def compute_quality(PUTATIVE_FILE, IMGT_CLASSIFIED, cutoff=None):
    all_genes = set()

    logging.info("Opening file")
    with open(IMGT_CLASSIFIED, 'r') as inputfile:
        lines = inputfile.readlines()

    logging.info("Parsing lines")
    for line in tqdm.tqdm(lines):
        if line[0] == '>': 
            all_genes.add(line.split(';')[1].split(' ')[0])

    found_genes = defaultdict(int)
    ROC_data = []
    PRC_data = []
    FP = 0.0
    TP = 0.0
    
    lines = open(PUTATIVE_FILE, 'r').readlines()
    if cutoff is not None:
        lines = lines[:cutoff]

    for line in tqdm.tqdm(lines):
        if line[0] == '>': 
            IMGT_class, matchlen, dist = line.split(';')[1].split(' ')
            if int(dist) < 4:
                found_genes[IMGT_class] += 1
                if found_genes[IMGT_class] <= 1:
                    TP += 1
            else:
                FP += 1
            
            FN = (len(all_genes) - len(found_genes)) * 1.0
            ROC_data.append([TP, FP])
            PRC_data.append([TP / (TP + FP), TP / (TP + FN)])

    ROC_raw = np.array(ROC_data, dtype='float32')
    PRC_raw = np.array(PRC_data, dtype='float32')
    
    ROC_raw[:, 0] /= len(all_genes)
    ROC_raw[:, 1] /= np.max(ROC_raw[:, 1])

    return ROC_raw, PRC_raw

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

    logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s %(levelname)s %(process)d: ' +
                    '%(message)s',
                    level=logging.INFO)

    RESULTS_DIR = '/Users/rishi/dbio/mnt/final-paper/human/091/'
    # RESULTS_DIR = '/data/projects/rishi/final-paper/human/091/'
    FULL_ALGO_GENES = os.path.join(RESULTS_DIR, 'qc/full-alg-v-segs.VDJ.H3.L3.CH1.dnaH3.fa')
    FREQ_ONLY_GENES = os.path.join(RESULTS_DIR, 'qc/most-freq-v-segs.VDJ.H3.L3.CH1.dnaH3.fa')
    IMGT_CLASSIFIED = os.path.join(RESULTS_DIR, 'input_seqs.fa')

    ROC_raw, PRC_raw = compute_quality(FULL_ALGO_GENES, IMGT_CLASSIFIED, cutoff=750)
    ROC_baseline, PRC_baseline = compute_quality(FREQ_ONLY_GENES, IMGT_CLASSIFIED, cutoff=750)

    evalCurve([(ROC_raw, 'Full Alg.', 'green'),
               (ROC_baseline, 'Freq. Only', 'red')], 
               os.path.join(RESULTS_DIR,'qc/roc.png'),
              title='ROC (Human 091)')

    evalCurve([(PRC_raw, 'Full Alg.', 'green'),
               (PRC_baseline, 'Freq. Only', 'red')],
               os.path.join(RESULTS_DIR,'prc.png'),
               title='PRC (Human 091)')
