import HTSeq
from scipy import signal
import numpy as np
import pickle
import pandas
import re
import os
import glob
import sys
from itertools import groupby
from operator import itemgetter

import scipy as sp
from peak import peak
import time
import logging

_debug = True
_verbose = True

logger = logging.getLogger(__name__)


def find_peaks(coverage, chrm='I', strand='+', config=None,
               use_raw_height_not_cwt=0):
    if config is not None and 'cwt_params' in config:
        params = config['cwt_params']
        (low_w, high_w, step_w, noise_perc) = params.split(':')
        (low_w, high_w, step_w) = (int(low_w), int(high_w), int(step_w))
        noise_perc = float(noise_perc)
    else:
        low_w = 25
        high_w = 325
        step_w = 100
        noise_perc = 0.05
    if config is not None and 'min_read_cutoff' in config:
        min_read_cutoff = int(config['min_read_cutoff'])
    if config is not None and 'use_raw_height_not_cwt' in config:
        use_raw_height_not_cwt = int(config['use_raw_height_not_cwt'])
    else:
        min_read_cutoff = 10  # FBF paper used a 10 here.
    logger.info("%s, %s..." % (chrm, strand))
    #if chrm != 'I':
    #    return []
    last_value = None
    for iv in coverage[chrm][strand].steps():
        last_value = iv[0].start
    if last_value is None or (last_value < 2): return []
    peak_pos = []
    width = 1e3
    start_time = time.time()
    n_windows = 0
    total_windows = last_value/int(width)
    logger.info('Chrm %s Strand %s Len %s' % (chrm, strand, last_value))
    f_max = np.max(np.fromiter(
        coverage[HTSeq.GenomicInterval(chrm, 0, min([1e6, last_value -1]), strand)],
        dtype=np.float
    ))
    logger.info('Max coverage in first 1e6 bp %f' % float(f_max))
    if f_max < 5:
        logger.warning('Low coverage!')
    for start in range(0, last_value, int(width)):
        n_windows += 1
        coverage_array = np.array([0])
        window = HTSeq.GenomicInterval(chrm, start, start + width, strand)
        coverage_array = np.fromiter(coverage[window], dtype=np.float)
        if use_raw_height_not_cwt:
            v = np.where(coverage_array >= use_raw_height_not_cwt,
                         range(len(coverage_arrcay)), 0 * len(coverage_array))
            for k, g in groupby(enumerate(v), lambda (i, x): i-x):
                _range = map(itemgetter(1), g)
                if len(_range) > 1 and _range[0] > 0:
                    _range = sorted(_range, key=lambda xpos: coverage_array[xpos])
                    peak_pos.append(start + _range[-1])
            continue
        if not n_windows % 1e3:
            cur_time = float(time.time() - start_time)/60.
            speed = cur_time/float(max([1., n_windows]))
            finish_genome = 10 * total_windows*speed
            logger.info(
                "time elapsed: %.3f m, estimate for chrom %.3f m, genome %.3f m, all replicates %.3f m, window %i/%i" % (
                    cur_time, total_windows*speed, finish_genome, finish_genome * 3, n_windows, total_windows)
            )
        if max(coverage_array) < min_read_cutoff:
            continue
        for a_pos in signal.find_peaks_cwt(
                coverage_array, np.arange(low_w, high_w, step_w), noise_perc=noise_perc):
            if coverage_array[a_pos] >= min_read_cutoff:
                peak_pos.append(start + a_pos)
    li = "peak_calling_tools.find_peaks(): Number of peaks called %i" % len(
        peak_pos)
    logger.info(li)
    logger.info('Windows used: %i' % n_windows)
    if n_windows < 10: logger.warn("Low window number!")
    return peak_pos


# From: http://stackoverflow.com/questions/17118350/how-to-find-nearest-value-that-is-greater-in-numpy-array
def argfind(array, predicate):
    for i in xrange(array.shape[0]):
        if predicate(array[i]):
            return i
    return False


def find_nearest_below(array, value):
    return argfind(array, lambda x: x <= value)


def find_borders(peak_pos, coverage, chrm, strand):
    # What are the borders of each peak?
    peak_objs = []
    for a_peak_pos in peak_pos:
        width = 1e3
        left = max(0, a_peak_pos - width)
        right = a_peak_pos + width
        window = HTSeq.GenomicInterval(chrm, left, right, strand)
        peak_center = HTSeq.GenomicPosition(chrm, a_peak_pos, strand)
        if coverage[peak_center] < 5:
            continue
        wincvg = np.fromiter(coverage[window], dtype='i')
        cutoff = coverage[peak_center]*0.2#max(wincvg)*.1
        rv = wincvg[::-1]
        left_border = find_nearest_below(rv[1000:], cutoff)
        right_border = find_nearest_below(wincvg[1000:], cutoff)
        left_border = a_peak_pos - left_border
        right_border = a_peak_pos + right_border
        if left_border == right_border:
            left_border -= 20
            right_border += 20
        peak_obj = peak(chrm, left_border, right_border, strand)
        window = HTSeq.GenomicInterval(chrm, left_border, right_border, strand)
        wincvg = np.fromiter(coverage[window], dtype='i')
        peak_obj.height = int(max(wincvg))
        if peak_obj.height == 'na':
            print "Set a peak height to 'na': %s..." % str(peak_obj.__dict__)
        peak_objs.append(peak_obj)
    return peak_objs


# Merge overlapping peaks.
# From: http://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
def merge_overlapping_on_chrm_and_strand(intervals, coverage):
    """Merge in a given chrom and strand.
    """
    sorted_by_lower_bound = sorted(intervals, key=lambda x: x.left)
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if higher.left <= lower.right:
                upper_bound = int(max(lower.right, higher.right))
                new_peak = peak(lower.chrm, lower.left, upper_bound, lower.strand)
                new_peak.height = 0
                window = HTSeq.GenomicInterval(lower.chrm, lower.left, upper_bound, lower.strand)
                wincvg = np.fromiter(coverage[window], dtype='i')
                new_peak.height = int(max(wincvg))
                merged[-1] = new_peak  # replace by merged interval
            else:
                merged.append(higher)
    return merged

def merge_overlapping(intervals):
    by_chrm = {}
    merged = {}
    for apeak in intervals:
        by_chrm.setdefault(apeak.chrm, {})
        merged.setdefault(apeak.chrm, {})
        by_chrm[apeak.chrm].setdefault(apeak.strand, [])
        merged[apeak.chrm].setdefault(apeak.strand, [])
        by_chrm[apeak.chrm][apeak.strand].append(apeak)
    for chrm in by_chrm:
        for strand in by_chrm[chrm]:
            merged[chrm][strand] = merge_overlapping_on_chrm_and_strand(by_chrm[chrm][strand])
            # Check.
            check_overlap(merged[chrm][strand])
    return merged


def call_peaks_from_wig(coverage, fname, config, load_data=False):
    gtf_filename = config['gtf_filename']
    datafile = 'data/%s/peak_objs_by_chrm_%s.p' % (
                config['experiment_name'], os.path.basename(fname))
    if load_data:
        logger.info(
            "call_peaks_from_wig(): Loading existing peak_objs_by_chrm.p object.")
        with open(datafile, 'rb') as f:
            peak_objs_by_chrm = pickle.load(f)
        return peak_objs_by_chrm
    gtf_df = pandas.read_csv(gtf_filename, sep='\t')
    peaks_by_chrm = {}
    peak_objs_by_chrm = {}
    for chrm in dict(gtf_df['0'].value_counts()).keys():
        peaks_by_chrm[chrm] = {}
        peak_objs_by_chrm[chrm] = {}
        try:
            a = coverage[chrm]
        except:
            logger.warn("%s chrom not in coverage" % chrm)
            continue
        for strand in ['+', '-']:
            peaks_by_chrm[chrm][strand] = find_peaks(coverage, chrm=chrm, strand=strand, config=config)
            peak_objs_by_chrm[chrm][strand] = find_borders(
                peaks_by_chrm[chrm][strand], coverage, chrm, strand)
            peak_objs_by_chrm[chrm][strand] = merge_overlapping_on_chrm_and_strand(
                peak_objs_by_chrm[chrm][strand], coverage)
            assign_to_gene(peak_objs_by_chrm, chrm, strand, gtf_df)
            logger.info('Called %i peaks on chrm %s strand %s' % (
                len(peak_objs_by_chrm[chrm][strand]), chrm, strand
            ))
    if not os.path.exists('data/'): os.system('mkdir data')
#    if not os.path.exists('data/pickled_peak_obj'): os.system('mkdir data/pickled_peak_obj')
    with open(datafile, 'wb') as f:
        pickle.dump(peak_objs_by_chrm, f)
    return peak_objs_by_chrm


def convert_peak_objs_to_table(peak_objs_by_chrm):
    # Convert to a table and write.
    peak_list = []
    output_cols = ['chrm', 'left', 'right', 'strand',
                   'height', 'gene_name']  # Simple columns.
    for chrm in peak_objs_by_chrm:
        for strand in peak_objs_by_chrm[chrm]:
            for p in peak_objs_by_chrm[chrm][strand]:
                if not hasattr(p, 'gene_name'): p.gene_name = 'NA'
                row = [p.chrm, p.left, p.right, p.strand, p.height, p.gene_name]
                peak_list.append(row)
    peak_table = pandas.DataFrame(peak_list, columns=output_cols)
    peak_table.to_csv('peak_table.txt', sep='\t')
    return peak_table


def load_bed_folder(folder_name, config=None):
    logger.info("Loading .bed files in %s/" % folder_name)
    ga = {}
    if config is not None:
        exp_name = config['experiment_name']
        exp_name = re.sub('\w+', '', exp_name)  # Delete numbers for no reason.
    for fname in glob.glob(folder_name + '/*.bed'):
        exp = os.path.basename(fname).partition('.bed')[0]
        if re.match('fog', exp) is not None: comb_exp = 'fog'
        elif re.match('control', exp) is not None: comb_exp = 'control'
        elif re.match('n2', exp) is not None: comb_exp = 'n2'
	elif re.match('.*' + exp_name + '.*', exp) is not None:
            comb_exp = exp_name
        else: comb_exp = 'Unknown'
        ga.setdefault(exp, HTSeq.GenomicArrayOfSets(chroms='auto', stranded=True))
        ga[exp] = load_bed_file(fname)
    return ga


def load_bed_file(fname):
    if not os.path.exists(fname):
        fname = re.sub('fog3', 'fog', fname)
    logger.info("peak_calling_tools.load_bed_file({n}) running...".format(n=fname))
    #ga = HTSeq.GenomicArrayOfSets(chroms='auto', stranded=True)
    ga = HTSeq.GenomicArray(chroms='auto', stranded=True)
#    os.system('wc -l {c} > tmp'.format(c=fname))
    start_time = time.time()
    with open(fname, 'r') as f:
        for n, line in enumerate(f):
            s = line.rstrip('\n').split('\t')
            ga[HTSeq.GenomicInterval(
                s[0], int(s[1]), int(s[2]), s[5])] += 1
            if not n % 1e3: print "Loading {i}: line {n}.".format(i=fname, n=n)
    li = "\tTook %.3f m to read bed file." % float((time.time() - start_time)/60.)
    logger.info(li)
    return ga


def load_bedgraph_file(fname, add_strand=True):
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


def load_bedgraph_folder(folder_name):
    ga = {}
    for fname in glob.glob(folder_name + '/*.wig'):
        exp = os.path.basename(fname).partition('.bed')[0]
        ga[exp] = load_bedgraph_file(fname)
    return ga


def find_reproducible_peaks():
    pass


def fill_in_peaks(coverage, peak_table, output_filename='all_peaks.txt'):
    _table = []
    logger.info("Filling in peak information...")
    # Create object of peak ranges with which to count.
    n = 0
    ht_peaks = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    for index, row in peak_table.iterrows():
        n += 1
        iv = HTSeq.GenomicInterval(row['chrm'], row['left'], row['right'], row['strand'])
        ht_peaks[iv] += str(n)
    # Count reads in peak ranges.
    n = 0
    for index, row in peak_table.iterrows():
        n += 1
        _table.append({
            'name': n,
            'gene_name': row['gene_name'],
            'chrm': row['chrm'], 'left': row['left'],
            'right': row['right'], 'strand': row['strand']})
        iv = HTSeq.GenomicInterval(
            row['chrm'], int(row['left']), int(row['right']), row['strand'])
        for exp in coverage:
            if exp == 'fog': continue
            if exp == 'control': continue
            read_nums = set()
            for iv, val in coverage[exp][iv].steps():
                read_nums |= val
            counts = len(list(read_nums))
            _table[-1][exp] = counts
            #_table[-1][exp] = np.max(np.fromiter(coverage[exp][iv], dtype='i'))
    df = pandas.DataFrame(_table)
    if not os.path.exists('filled_counts/'): os.system('mkdir filled_counts')
    df.to_csv(output_filename, sep='\t', index=False)
    for exp in coverage:
        if exp == 'fog': continue
        if exp == 'control': continue
        with open('filled_counts/%s.txt' % exp, 'w') as f:
            for row in _table:
                f.write('%s\t%s\n' % (row['name'], row[exp]))


def check_overlap(intervals):
    sorted_by_lower_bound = sorted(intervals, key=lambda x: x.left)
    merged = []
    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher.left <= lower.right:
                print "Overlapping peak ranges..\n%s\n%s\n*" % (str(lower), str(higher))


def get_distance_pandas(_gene, apeak):
    """Used to assign peaks to a gene.
    """
    if _gene['4'] < apeak.left:
        return apeak.left - _gene['4']
    if apeak.right < _gene['3']:
        return apeak.right - _gene['3']
    return 0


def find_nearby_genes(sub, apeak):
    """Assigns a gene for a given peak.
    """
    cutoff = 1e3
    genes_for_peak = []
    for index, _gene in sub.iterrows():
        dist = get_distance_pandas(_gene, apeak)
        if abs(dist) < cutoff:
            genes_for_peak.append(_gene)
    return genes_for_peak


def resolve_tie(ties):
    biotypes = []
    for tup in ties:
        _gene = tup[1]
        m = re.search('gene_biotype "([^"]+)"', _gene['8'])
        if m is not None:
            biotypes.append((tup[0], tup[1], m.group(1)))
        else:
            biotypes.append((tup[0], tup[1], 'Unknown'))
    # If there is a non-coding RNA, assign to that.
    non_coding = []
    for tup in biotypes:
        if tup[2] != 'protein_coding':
            non_coding.append(tup)
    if len(list(non_coding)) == 1:
        return (non_coding[0][0], non_coding[0][1])
    if len(list(non_coding)) > 1:
        # Multiple non-coding RNAs. Pick randomly.
        return (non_coding[0][0], non_coding[0][1])
    if len(list(non_coding)) == 0:
        # No non-coding RNA. Pick randomly.
        return (ties[0][0], ties[0][1])


# Find the closet gene for each peak, and resolve ties.
def assign_to_gene(merged_peaks, chrm, strand, gtf_df, verbose=False):
    if chrm not in merged_peaks:
        return False
    if strand not in merged_peaks[chrm]:
        return False
    sub = gtf_df[(gtf_df['0']==chrm) & (gtf_df['6']==strand)]
    sub = sub[sub['2']=='transcript']
    for index, apeak in enumerate(merged_peaks[chrm][strand]):
        if not index % 1e3:
            logger.info("Assigning %i genes for chrm/strand %s/%s" % (
                len(merged_peaks[chrm][strand]), chrm, strand))
        asub = sub[(abs(sub['3'] - apeak.left) < 1e5
                    ) | (abs(sub['4'] - apeak.left) < 1e5)]
        apeak.genes_for_peak = find_nearby_genes(asub, apeak)
        closest = (1e4, None)
        ties = []
        for _gene in apeak.genes_for_peak:
            dist = get_distance_pandas(_gene, apeak)
            if abs(dist) < abs(closest[0]):
                ties = []
                closest = (dist, _gene)
            elif abs(dist) == abs(closest[0]):
                ties.append((dist, dict(_gene)))
                ties.append((dist, dict(closest[1])))
        if len(ties) > 0:
            for tup in ties:
                dist = tup[0]
                _gene = tup[1]
                #dist = get_distance_pandas(_gene, apeak)
                if dist == closest[0]:
                    closest = resolve_tie(ties)
                    break
        apeak.gene_dist = closest[0]
        apeak.gene = closest[1]
        if closest[1] is not None:
            apeak.gene_name = closest[1]['gene_name']
        else:
            apeak.gene_name = 'NA'
        del apeak.genes_for_peak


def load_bedgraphs_and_call_peaks(in_dir, output_filename='all_peaks.txt',
                                  reads_dir='bed_collapsed/'):
    config = {'gtf_filename': 'lib/gtf_with_names_column.txt',
              'gtf_filename_noheader': 'lib/gtf_with_names_column_no_header.txt'
    }
    ga = load_bedgraph_folder(in_dir)
    if not os.path.exists('peaks/'): os.system('mkdir peaks')
    peak_objs_by_chrm = call_peaks_from_wig(ga['fog'], 'fog', config)
    peak_table = convert_peak_objs_to_table(peak_objs_by_chrm)
    #with open('peak_objs_by_chrm_%s.p' % os.path.basename(clip_bam_filename), 'wb') as f:
    #    pickle.dump(peak_objs_by_chrm, f)
    peak_table.to_csv('peaks/%s.txt' % 'fog', sep='\t')
    ga_raw = load_bed_folder(reads_dir)
    fill_in_peaks(ga_raw, peak_table, output_filename=output_filename)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='''Peak calling tools.''')
    parser.add_argument(
        '-b', '--bedgraph', default='wig/',
        help='''Folder of bedgraph files of raw coverage data.'''
    )
    parser.add_argument(
        '-o', '--output_peaks_file', default='all_peaks.txt',
        help='''Output filename for peaks called from bedgraphs.'''
    )
    parser.add_argument(
        '-r', '--reads_dir', default='../../pre_calls/bed_collapsed/',
        help='Directory of bed files of all reads.'
    )
    args = parser.parse_args()
    #load_bed_folder(args.bedgraph)
    load_bedgraphs_and_call_peaks(
        args.bedgraph + '/', output_filename=args.output_peaks_file,
        reads_dir=args.reads_dir)
