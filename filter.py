"""
The foolish original definitions of the configfile are felt here:
the replicates were given as paths to .wig files (as though nonstranded
bedgraphs were used, but they aren't), and the control datasets to
bed files (that do exist), and every path is sort of filled out backwards
from that. 

A second version of the config file settings would be the config.ini format:

[library]

top: ./
lib: /groups/Kimble/Common/fog_iCLIP/calls/lib/
gtf_raw: %(lib)s/Caenorhabditis_elegans.WBcel235.78.noheader.gtf
fasta: %(lib)s/c_elegans.WS235.genomic.fa
fai: %(lib)s/c_elegans.WS235.genomic.fa.fai
chr_sizes: %(lib)s/chr_sizes.txt
gtf_one_txpt_per_gene: %(lib)s/gtf_one_txpt_per_gene.txt
gtf: %(lib)s/gtf_with_names_column.txt

bedgraphs_folder: %(top)s/../bedgraphs_unnorm/combined_unnorm/

# Used by subpeaks.py
bedgraph_exp_plus:  %(bedgraphs_folder)s/exp_+.wig
bedgraph_exp_minus: %(bedgraphs_folder)s/exp_-.wig
read_beds:  %(top)s/../bed_collapsed/combined/
#
figs: %(top)s/figs/
#
control_bed1: control_n2.bed
#control_bed2: n2_oo_lane1_rt16.bed
#control_bed3: n2_oo_lane1_rt3.bed
exp_bed1: fbf_rep_1.bed
exp_bed2: fbf1_oo_lane2_rt6.bed
exp_bed3: fbf1_oo_lane2_rt9.bed
#
clusters_dir: %(top)s/clusters/
permutation_peaks_dir: %(top)s/permutation_peaks/
"""
import re
import os
import sys
import pandas
import argparse
import HTSeq
import collections
import numpy as np

def get_val(ga, iv):
    return np.max(np.fromiter(ga[iv], dtype=np.float))



def add_heights_to_peak_file(peak_fname, bed_ga_dict):
    peaks = pandas.read_csv(peak_fname, sep='\t')
    ivs = zip(peaks['chrm'].tolist(), peaks['left'].tolist(),
        peaks['right'].tolist(), peaks['strand'].tolist())
    print ivs
    for bedgraph_name in bed_ga_dict:
        peaks[bedgraph_name] =  [
            get_val(bed_ga_dict[bedgraph_name], HTSeq.GenomicInterval(*iv)
                    ) for iv in ivs]
    return peaks


def load_bedgraph_file(fname, add_strand=True):
    fname = fname.rstrip('.bed')
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
    if add_strand:
        plus_file = fname.partition('.wig')[0] + '_+.wig'
        add_strand_to_ga_from_bedgraph_file(plus_file, ga, '+')
        minus_file = fname.partition('.wig')[0] + '_-.wig'
        add_strand_to_ga_from_bedgraph_file(minus_file, ga, '-')
    else:
        if re.match('.*\+.*', os.path.basename(fname)) is not None:
            add_strand_to_ga_from_bedgraph_file(fname, ga, '+')
        elif re.match('.*-.*',  os.path.basename(fname)) is not None:
            add_strand_to_ga_from_bedgraph_file(fname, ga, '-')
    return ga


def add_strand_to_ga_from_bedgraph_file(fname, ga, strand):
    with open(fname, 'r') as f:
        next(f)
        for line in f:
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), strand)] = float(s[3])
    return ga


def read_config(filename):
    """Expect:
experiment_name\tname  # Optional.
clip_replicate\tfilename1.wig
clip_replicate\tfilename2.wig
clip_replicate\tfilename3.wig
gtf_filename\tgtf_filename
rna_seq_filename\tfilename
neg_ip_filename\tfilename
positive_control_genes\tname1;name2;name3;
motif\ttgt\w\w\wat
    """
    config = collections.defaultdict(list)
    with open(filename, 'r') as f:
        for li in f:
            li = li.partition('#')[0]  # Skip comments.
            if li == '': continue
            s = li.rstrip('\n').split('\t')
            try: config[s[0]].append(s[1])
            except: print "Error parsing line %s in config file. Wrong length?" % li
    for key in config:
        if len(config[key]) == 1:
            config[key] = config[key][0]
        if key == 'positive_control_genes':
            config[key] = config[key].split(';')
    if 'experiment_name' not in config:
        config['experiment_name'] = os.path.basename(filename)
    return config

def get_bed_size(fname):
    return float(len(open(fname).readlines()))


def add_read_columns(args, config):
    non_control = {}
    bedfiles = {'control': "{a}/{b}.wig".format(
        a=os.path.dirname(config['clip_replicate'][0]),
        b=os.path.basename(config['neg_ip_filename']).rstrip('.bed'))}
    for x in config['clip_replicate']:
        name = os.path.basename(x).rstrip('.bed').rstrip('.wig')
        bedfiles[name] = x
        non_control[name] = x
    ga_d = {}
    print "Loading bedgraphs..."
    for name in bedfiles:
        ga_d[name] = load_bedgraph_file(bedfiles[name])
    peaks = add_heights_to_peak_file(args.peaks_fname, ga_d)
    size = {}
    for bedgraph in bedfiles:
        if re.search('control', bedgraph):
            size[bedgraph] = get_bed_size(config['neg_ip_filename'])
            col_name = 'depth_control_' + bedgraph
            peaks[col_name] = [
                1e6 * x/size[bedgraph] for x in peaks[bedgraph].tolist()]
            continue
        print config['bed_dir']
        print bedgraph
        print os.path.basename(bedgraph).rstrip('.wig')
        bed_file = config['bed_dir'] +'/' + \
                   os.path.basename(bedgraph).rstrip('.wig') + '.bed'
        size[bedgraph] = get_bed_size(bed_file)
        col_name = 'depth_exp_' + bedgraph
        peaks[col_name] = [
            1e6 * x/size[bedgraph] for x in peaks[bedgraph].tolist()]
    return peaks


def add_sum_and_ratio_columns(peaks):
    exp_cols = [col for col in peaks.columns if re.search('depth_exp_', col)]
    control_cols = [col for col in peaks.columns\
                    if re.search('depth_control_', col)]
    for i, row in peaks.iterrows():
        peaks.loc[i, 'exp'] = sum([
            peaks.loc[i, col] for col in exp_cols])/float(len(exp_cols))
        peaks.loc[i, 'control'] = sum([
            peaks.loc[i, col] for col in control_cols])/float(len(control_cols))
        peaks.loc[i, 'ratio'] = float(peaks.loc[i, 'exp'])\
                                /float(max([1., row['control']]))
    return peaks

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--peaks_fname',
                        help='Filename of peaks file.')
    parser.add_argument('-c', '--config')
    parser.add_argument('-r', '--ratio_cutoff',
                        help='Enrichment ratio cutoff.')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    args.ratio_cutoff = float(args.ratio_cutoff)
    config = read_config(args.config)
    header = open(args.peaks_fname).readline()
    if (re.search('depth_exp', header) is not None) and (
        re.search('depth_control', header) is not None):
        print "Peaks file appears to already have read count columns.\
 Overwrite them?"
        answer = raw_input(prompt)
        print answer
        if answer[0].upper() == 'Y':
            peaks = add_read_columns(args, config)
            peaks = add_sum_and_ratio_columns(peaks)
            peaks.to_csv(args.peaks_fname, sep='\t', index=False)
        else:
            print "Using the existing columns then..."
    else:
        peaks = add_read_columns(args, config)
        peaks = add_sum_and_ratio_columns(peaks)
        peaks.to_csv(args.peaks_fname, sep='\t', index=False)
    peaks = peaks[peaks['ratio']>=args.ratio_cutoff]
    peaks.to_csv(args.output, sep='\t', index=False)

