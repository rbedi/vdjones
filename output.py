from collections import defaultdict

def dump_connectivity(adj_list, out_file):
    with open(out_file, 'w') as of:
        for v in adj_list:
            for nb in adj_list[v]:
                of.write('{}\t{}\t{}\n'.format(v, nb, 1))


def dump_fasta(seqs, out_file, data=None):
    with open(out_file, 'w') as of:
        for idx, v in enumerate(seqs):
            if data is None:
                of.write('>{}\n{}\n'.format(idx, v))
            else:
                of.write('>{};{}\n{}\n'.format(idx, data[idx], v))

def dump_list(seqs, out_file):
    with open(out_file, 'w') as of:
        for idx, v in enumerate(seqs):
            of.write('{}\n'.format(v))

def recover_adj_list(in_file):
    adj_list = defaultdict(list)
    with open(in_file, 'r') as inf:
        lines = inf.readlines()
        for line in lines:
            parts = line.strip().split()
            adj_list[int(parts[0])].append(int(parts[1]))
    return adj_list

