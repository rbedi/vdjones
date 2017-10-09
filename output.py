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
