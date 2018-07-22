import HTSeq
import numpy as np
import pickle
import pandas
import re
import os
import glob
import sys
import scipy as sp
import peak_calling_tools
from statsmodels.sandbox.stats.multicomp import multipletests
import logging
import time

logger = logging.getLogger(__name__)

_debug = True
_verbose = True


def do_statistics(
        peak_objs, config, clip_replicate_filename,
        nb_fits=None, skip_nb=False):
    bedfiles = {
            'clip': clip_replicate_filename,
            'rna_seq': config['rna_seq_filename'],
            'neg_ip': config['neg_ip_filename']}
    logger.info("do_statistics() with %i peaks..." % len(peak_objs))
    if len(peak_objs) == 0:
        logger.warn('No peaks called!')
        return {}
    logger.info("Peak 1: %s" % str(peak_objs[0].__dict__))
    for bedfile in list(bedfiles.values()):
        if not os.path.exists(bedfile): logger.error('Expected file %s does not exist.' % bedfile)
    bedfile_sizes = read_bedfile_sizes(bedfiles)
    norm_factors = get_norm_files(bedfile_sizes)
    replicate_start_time = time.time()
    block_start_time = time.time()
    if nb_fits is None: nb_fits = {}
    nb_pvals_gene = {'neg_ip': [], 'clip': [], 'rna_seq': []}
    nb_pvals_local = {'neg_ip': [], 'clip': [], 'rna_seq': []}
    for n, peak in enumerate(peak_objs):
        nb_fits.setdefault(peak.gene_name, {'Exons': {}, 'Locals': {}})
        #if peak.gene_name != 'gld-1': continue
        peak.find_max_bin()
        peak.local_poisson(bedfiles, norm_factors)
        peak.local_norm(bedfiles, norm_factors)
        res = peak.local_nb(bedfiles, norm_factors, known=nb_fits, skip=skip_nb)
        for key in res:
            nb_pvals_local[key].append(res[key])
        peak.gene_poisson(bedfiles, norm_factors)
        peak.gene_norm(bedfiles, norm_factors)
        res = peak.gene_nb(bedfiles, norm_factors, known=nb_fits, skip=skip_nb)
        for key in res:
            nb_pvals_gene[key].append(res[key])
        if not n % 1e2:
            if n == 0: continue
            exp_time_min = (time.time() - replicate_start_time)/60.
            block_time_min = (time.time() - block_start_time)/60.
            time_per_peak_min = block_time_min/100.
            time_to_go_hours = (len(peak_objs) - n) * time_per_peak_min/60.
            logger.info("do_statistics(): Example of peak.pvalue: " + str(peak.pvalues))
            logger.info('''Time expired (replicate {r}, peak {pk_num}/{total}): {t} m.
            Per peak {pp} m. To go {g} h.'''.format(
                r=clip_replicate_filename, pk_num=n, total=len(peak_objs),
                t=exp_time_min,
                pp=time_per_peak_min, g=time_to_go_hours
            ))
            logger.info('''Local NB p values: neg_ip {n} rna seq {r} clip {c}'''.format(
                n=percentage_split(nb_pvals_local['neg_ip'], [0.25, 0.25, 0.25, 0.25]),
                r=percentage_split(nb_pvals_local['rna_seq'], [0.25, 0.25, 0.25, 0.25]),
                c=percentage_split(nb_pvals_local['clip'], [0.25, 0.25, 0.25, 0.25])
            ))
            logger.info('''By gene NB p values: neg_ip {n} rna seq {r} clip {c}'''.format(
                n=percentage_split(nb_pvals_gene['neg_ip'], [0.25, 0.25, 0.25, 0.25]),
                r=percentage_split(nb_pvals_gene['rna_seq'], [0.25, 0.25, 0.25, 0.25]),
                c=percentage_split(nb_pvals_gene['clip'], [0.25, 0.25, 0.25, 0.25])
            ))
            logger.info('Genes in nb_fits {s}.'.format(
                s=len(nb_fits)
            ))
            block_start_time = time.time()
    return nb_fits

def percentage_split(seq, percentages):
    # http://stackoverflow.com/questions/14280856/separate-a-list-into-four-parts-based-on-percentage-even-if-the-list-is-not-divi
    cdf = np.cumsum(percentages)
    seq = sorted(seq)
    assert cdf[-1] == 1.0  # Test if percentages sum to 1.0
    stops = list(map(int, cdf * len(seq)))
    sep_l = [seq[a:b] for a,b in zip([0]+stops, stops)]
    return [np.median(x) for x in sep_l]

def get_norm_files(bedfile_sizes):
    norm_factors = {}
    norm_factors['clip'] = 1.0
    norm_factors['neg_ip'] = float(bedfile_sizes['clip'])/float(
        bedfile_sizes['neg_ip'])
    norm_factors['rna_seq'] = float(bedfile_sizes['clip'])/float(
        bedfile_sizes['rna_seq'])
    logger.info('get_norm_files(): bedfile_sizes: {bs}\nnorm_factors: {nm}'.format(
        bs=str(bedfile_sizes), nm=str(norm_factors)
    ))
    return norm_factors


def read_bedfile_sizes(bedfiles):
    bedfile_sizes = {}
    for bedfile in bedfiles:
        print('reading %s' % bedfile)
        fname = bedfiles[bedfile]
        plus_file = fname.partition('.wig')[0] + '_+.wig'
        minus_file = fname.partition('.wig')[0] + '_-.wig'
        bedfile_sizes[bedfile] = read_bedfile_size_on_strand(fname)
        #bedfile_sizes[bedfile] += read_bedfile_size_on_strand(minus_file)
    return bedfile_sizes


def read_bedfile_size_on_strand(fname):
    with open(fname, 'r') as f: wc_out = len(f.readlines())
    return wc_out


def fdr_correction(config, peak_table, clip_replicate_filename, alpha=0.01):
    # Need FDRs.
    fdrs = {}
    bamfiles = {
            'clip': config['bed_dir'] + '/' + os.path.basename(clip_replicate_filename).partition('wig')[0] + 'bed',
            'rna_seq': config['rna_seq_filename'],
            'neg_ip': config['neg_ip_filename']}
    for bamfile in bamfiles:
        # Local Poisson (only use CLIP).
        fdrs["%s_local_poisson" % bamfile] = multipletests(peak_table["%s_local_poisson" % bamfile].astype(float),
                                      alpha=alpha, method='fdr_bh')
        #print sorted(fdrs["%s_local_poisson" % bamfile][0])
        #print sorted(fdrs["%s_local_poisson" % bamfile][1])

        peak_table["%s_local_poisson_rej" % bamfile] = pandas.Series(
            fdrs["%s_local_poisson" % bamfile][0], index=peak_table.index)
        peak_table["%s_local_poisson_cor" % bamfile] = pandas.Series(
            fdrs["%s_local_poisson" % bamfile][1], index=peak_table.index)
        # Gene Poisson (only use CLIP).
        fdrs["%s_gene_poisson" % bamfile] = multipletests(peak_table["%s_gene_poisson" % bamfile].astype(float),
                                      alpha=alpha, method='fdr_bh')
        peak_table["%s_gene_norm_rej" % bamfile] = pandas.Series(
            fdrs["%s_gene_poisson" % bamfile][0], index=peak_table.index)
        peak_table["%s_gene_poisson_cor" % bamfile] = pandas.Series(
            fdrs["%s_gene_poisson" % bamfile][1], index=peak_table.index)
        # Local normal.
        fdrs["%s_local_norm" % bamfile] = multipletests(peak_table["%s_local_norm" % bamfile].astype(float),
                                      alpha=alpha, method='fdr_bh')
        peak_table["%s_local_norm_rej" % bamfile] = pandas.Series(
            fdrs["%s_local_norm" % bamfile][0], index=peak_table.index)
        peak_table["%s_local_norm_cor" % bamfile] = pandas.Series(
            fdrs["%s_local_norm" % bamfile][1], index=peak_table.index)
        # Gene normal.
        fdrs["%s_gene_norm" % bamfile] = multipletests(peak_table["%s_gene_norm" % bamfile].astype(float),
                                      alpha=alpha, method='fdr_bh')
        peak_table["%s_gene_norm_rej" % bamfile] = pandas.Series(
            fdrs["%s_gene_norm" % bamfile][0], index=peak_table.index)
        peak_table["%s_gene_norm_cor" % bamfile] = pandas.Series(
            fdrs["%s_gene_norm" % bamfile][1], index=peak_table.index)
        try:
            # Local negative binomial.
            fdrs["%s_local_nb" % bamfile] = multipletests(peak_table["%s_local_nb" % bamfile].astype(float),
                                          alpha=alpha, method='fdr_bh')
            peak_table["%s_local_nb_rej" % bamfile] = pandas.Series(
                fdrs["%s_local_nb" % bamfile][0], index=peak_table.index)
            peak_table["%s_local_nb_cor" % bamfile] = pandas.Series(
                fdrs["%s_local_nb" % bamfile][1], index=peak_table.index)
            # Gene negative binomial.
            fdrs["%s_gene_nb" % bamfile] = multipletests(peak_table["%s_gene_nb" % bamfile].astype(float),
                                          alpha=alpha, method='fdr_bh')
            peak_table["%s_gene_nb_rej" % bamfile] = pandas.Series(
                fdrs["%s_gene_nb" % bamfile][0], index=peak_table.index)
            peak_table["%s_gene_nb_cor" % bamfile] = pandas.Series(
                fdrs["%s_gene_nb" % bamfile][1], index=peak_table.index)
        except:
            print("Skipping NB FDR calculations...")

def evaluate_hypothesis(peak_table, clip_bed_filename, config, alpha=0.01):
    '''Apply an FDR cutoff and write output files.
    Input:
    peak_table: a pandas.DataFrame of peaks.
    clip_bed_filename: the name of the folder to write to in exp_name/peaks/
    config: config dict.
    Output: files named exp_name/peaks/replicate_name/null_hyp_*.txt
    '''
    out_dir = config['experiment_name'] + '/peaks/'
    if not os.path.exists(config['experiment_name']):
        os.system('mkdir ' + config['experiment_name'])
    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)
    replicate_dir = out_dir + os.path.basename(clip_bed_filename).partition('.')[0]
    if not os.path.exists(replicate_dir):
        os.system('mkdir ' + replicate_dir)
#    with open(replicate_dir + '/%s' % (config['experiment_name'], label), 'w') as f:
#        combined.to_csv(f, sep='\t')
    # Null hypothesis 2: Signal is not a local peak relative to itself.
    sub = peak_table[peak_table['clip_local_poisson_cor']<alpha]
    sub.to_csv('%s/null_hyp_2.txt' % replicate_dir, sep='\t', index=False)
    # Null hypothesis 3: Signal is not a peak relative to CLIP signal in the gene.
    sub = peak_table[peak_table['clip_gene_poisson_cor']<alpha]
    sub.to_csv('%s/null_hyp_3.txt' % replicate_dir, sep='\t', index=False)

    # Null hypothesis 4: Signal is not enriched relative to the negative in the local region.
    sub = peak_table[peak_table['neg_ip_local_norm_cor']<alpha]
    sub.to_csv('%s/null_hyp_4.txt' % replicate_dir, sep='\t', index=False)
    before = sorted(peak_table['neg_ip_local_norm_cor'].tolist())
    after = sorted(sub['neg_ip_local_norm_cor'].tolist())
    # Null hypothesis 5: Signal is not enriched relative to the negative for the gene.
    sub = peak_table[peak_table['neg_ip_gene_norm_cor']<alpha]
    sub.to_csv('%s/null_hyp_5.txt' % replicate_dir, sep='\t', index=False)
    try:
        # Null hypothesis 6: Signal is not enriched relative to Neg IP locally.
        sub = peak_table[peak_table['neg_ip_local_nb_cor']<alpha]
        sub.to_csv('%s/null_hyp_6.txt' % replicate_dir, sep='\t', index=False)
        # Null hypothesis 7: Signal is not enriched relative to Neg IP for the gene.
        sub = peak_table[peak_table['neg_ip_gene_nb_cor']<alpha]
        sub.to_csv('%s/null_hyp_7.txt' % replicate_dir, sep='\t', index=False)
    except:
        print("Skipping NB file output...")
        tmp = pandas.DataFrame(columns=peak_table.columns)
        tmp.to_csv('%s/null_hyp_6.txt' % replicate_dir, sep='\t', index=False)
        tmp.to_csv('%s/null_hyp_7.txt' % replicate_dir, sep='\t', index=False)
    # Null hypothesis 8: Signal is not enriched relative to RNA-seq locally.
    sub = peak_table[peak_table['rna_seq_local_norm_cor']<alpha]
    sub.to_csv('%s/null_hyp_8.txt' % replicate_dir, sep='\t', index=False)
    # Null hypothesis 9: Signal is not enriched relative to RNA-seq for the gene.
    sub = peak_table[peak_table['rna_seq_gene_norm_cor']<alpha]
    sub.to_csv('%s/null_hyp_9.txt' % replicate_dir, sep='\t', index=False)


