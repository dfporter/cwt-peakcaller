#import callpeaks
from callpeaks import get_sequences
from rc import rc
import config
import glob
import os
import pandas
import re
import argparse


def score_metrics(dir_name, config, include_ncrna=False):
    """Score metrics for every peaks file in a directory.
dir_name: directory name
config: dict of lib info
    """
    if os.path.isdir(dir_name):
        fname_list = glob.glob(dir_name + '/*.txt')
    elif os.path.isfile(dir_name):
        fname_list = [dir_name]
    else:
        fname_list = []
    li = ""
    for filename in sorted(fname_list, key=lambda x: os.path.basename(x)):
        print "Scoring metrics for %s" % filename
        li += score_metric(filename, config=config, include_ncrna=include_ncrna)
    with open('score_metrics_%s.txt' % config['experiment_name'], 'w') as f:
        f.write(li)
    return li


def score_metric(filename, label="", given_peaks=False, peaks=False,
                 config=None, include_ncrna=False):
    if not label:
        label = os.path.basename(filename)
    if not given_peaks:
        if len(open(filename).readlines()) < 2:
            return "No peaks."
        else:
            peaks = pandas.read_csv(filename, sep='\t')
            if not include_ncrna:
                if 'biotype' not in peaks.columns:
                    print('Asked to not count ncRNA, but not biotype column found')
                else:
                    len_all = len(peaks.index)
                    peaks = peaks[peaks['biotype']=='protein_coding']
                    print "Removed ncRNA: %i peaks input > %i after ncRNA removal" % (
                        len_all, len(peaks.index))
    if len(peaks.index) == 0:
        return "No peaks."
    get_sequences(peaks, fasta_filename=config['fasta'])
    score_binding_site(peaks, config=config)
    #run_dreme(peaks, label)
    positives = score_positives(peaks, config=config)
    return write_metrics(peaks, positives, label)


def write_metrics(peaks, positives, label):
    li = """
Dataset: {label}
Number of peaks: {df_size}
Number of genes: {n_genes}
With FBE: {with_fbe}, {fbe_perc}%
Without FBE: {without_fbe}
Positive controls: {observed}/{expected}
Missing positives: {missing}
""".format(label=label,
    df_size=len(peaks), n_genes=len(list(set(peaks['gene_name']))),
           with_fbe=len(peaks[peaks['has_motif']==1]),
           fbe_perc= float(100 * len(peaks[peaks['has_motif']==1])/len(peaks)),
           without_fbe=len(peaks[peaks['has_motif']==0]),
           observed=positives['number of observed positives'],
           expected=positives['expected'],
           missing=positives['missing positives'])
    return li


def score_binding_site(peaks, config=None):
    if config is None or 'motif' not in config:
        motif = 'tgt\w\w\wat'
    else:
        motif = config['motif'].lower()
    pat = re.compile(motif, re.IGNORECASE)
    for index, peak_row in peaks.iterrows():
        if pat.search(peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'has_motif'] = 1
        else:
            peaks.loc[index, 'has_motif'] = 0


def score_positives(peaks, config=None):
    if config is None or 'positive_control_genes' not in config:
        return {'observed positives': set([]), 'number of observed positives': 0,
            'missing positives': set([]), 'number of missing positives': 0,
            'expected': 0}
    positives = config['positive_control_genes']
    print 'score_positives'
    print positives
    known_pos = set(positives)
    obs_genes = set(peaks['gene_name'])
    obs_pos = known_pos & obs_genes
    missing_pos = known_pos - obs_genes
    obs_pos_n = len(list(obs_pos))
    missing_pos_n = len(list(missing_pos))
    return {'observed positives': obs_pos, 'number of observed positives': obs_pos_n,
            'missing positives': missing_pos, 'number of missing positives': missing_pos_n,
            'expected': len(list(known_pos))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--config', default='clip_config',
    )
    parser.add_argument(
        '-p', '--peak_dir', default='peaks',
    )
    parser.add_argument(
        '-n', '--include_ncrna',
        action='store_true', default=False
    )
    args = parser.parse_args()
    lib = config.config(filepath=args.config)
    score_metrics(args.peak_dir, lib, include_ncrna=args.include_ncrna)
