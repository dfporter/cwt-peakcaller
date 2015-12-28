import HTSeq
import numpy as np
import pickle
import pandas
import re
import os
import glob
import sys
import scipy as sp
from statsmodels.sandbox.stats.multicomp import multipletests
import argparse
import logging
import time
import datetime

import peak_calling_tools
import combine_replicates
import stats_to_call_peaks
import add_signal


def start_logger(exp_name):
    global logger
    logger = logging.getLogger(__name__)
    _debug = True
    if not os.path.exists('logs'): os.system('mkdir logs')
    logging.basicConfig(
        filename='logs/%s_%s_callpeaks.log' % (
            datetime.datetime.now().strftime('%dd%Hh%Mm%Ss'), exp_name),
        level=logging.DEBUG)
    logging.info('Module called %s' % str(time.localtime()))
    return logger


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
        cutoff = max(wincvg)*.1
        rv = wincvg[::-1]
        left_border = find_nearest_below(rv[1000:], cutoff)
        right_border = find_nearest_below(wincvg[1000:], cutoff)
        left_border = a_peak_pos - left_border
        right_border = a_peak_pos + right_border
        if left_border == right_border:
            left_border -= 10
            right_border += 10
        peak_obj = peak(chrm, left_border, right_border, strand)
        window = HTSeq.GenomicInterval(chrm, left_border, right_border, strand)
        wincvg = np.fromiter(coverage[window], dtype='i')
        peak_obj.height = int(max(wincvg))
        if peak_obj.height == 'na':
            logger.warn("Set a peak height to 'na': %s..." % str(peak_obj.__dict__))
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


def any_have_na(peak_objs):
    has_na = False
    for apeak in peak_objs:
        if apeak.height == 'na':
            logger.warn("Has na height: %s" % str(apeak.__dict__))
    if not has_na:
        logger.info("No peaks have na height.")


def call_peaks_by_cwt_on_replicate(clip_wig_filename, config, load_data=False):
    """Calls or loads CWT peaks for a wig filename.
    In: filename of .wig, config data.
    Writes: data/expname/*.p datafiles of CWT peaks.
    Writes: data/expname/cwt_calls/ tables of CWT peaks.
    Returns: pandas.DataFrame of the CWT table written to cwt_calls/
    """
    loaded_data = False
    if not os.path.exists('data/'): os.system('mkdir data')
    if not os.path.exists('data/%s' % config['experiment_name']):
        os.system('mkdir data/%s' % config['experiment_name'])
    datafile = 'data/%s/peak_objs_by_chrm_%s.p' % (
                config['experiment_name'], os.path.basename(clip_wig_filename))
    if load_data:
        try:
            with open(datafile, 'rb') as f:
                peak_objs_by_chrm = pickle.load(f)
            li = "Loaded peaks called previously and saved as file %s." % datafile
            logger.info(li)
            peak_table = convert_peak_objs_to_table(peak_objs_by_chrm)
            loaded_data = True
            print li
        except:
            li = "Failed to load data for %s (Expected %s). Regenerating..." % (
                clip_wig_filename,
                datafile)
            logger.warn(li)
            print li
            loaded_data = False
    if not loaded_data:
        coverage = peak_calling_tools.load_bedgraph_file(clip_wig_filename)
        peak_objs_by_chrm = peak_calling_tools.call_peaks_from_wig(
            coverage, clip_wig_filename, config)
        peak_table = convert_peak_objs_to_table(peak_objs_by_chrm)
        with open(datafile, 'wb') as f:
            pickle.dump(peak_objs_by_chrm, f)
    if not os.path.exists('data/%s/cwt_calls/' % config['experiment_name']):
        os.system('mkdir data/%s/cwt_calls' % config['experiment_name'])
    peak_table.to_csv('data/%s/cwt_calls/%s' % (
        config['experiment_name'], os.path.basename(clip_wig_filename)),
                      sep='\t', index=False)
    return peak_table


def convert_peak_objs_to_table(peak_objs_by_chrm):
    total_peaks = 0
    for chrm in peak_objs_by_chrm:
        for strand in peak_objs_by_chrm[chrm]:
            total_peaks += len(peak_objs_by_chrm[chrm][strand])
    logger.info(
        'Converting %i peak objects in %i chromosomes to a pandas.DataFrame.' % (
            total_peaks, len(peak_objs_by_chrm)))
    # Convert to a table and write.
    peak_list = []
    output_cols = ['chrm', 'left', 'right', 'strand',
                   'height']#, 'max_bin']  # Simple colums.
    for key in peak_objs_by_chrm.keys():
        if peak_objs_by_chrm[key]['+'] is not None:
            try:
                output_cols += dict(peak_objs_by_chrm[key]['+'][0].gene).keys()
                output_cols += dict(peak_objs_by_chrm[key]['+'][0].pvalues).keys()
                break
            except: return pandas.DataFrame()
    for chrm in peak_objs_by_chrm:
        for strand in peak_objs_by_chrm[chrm]:
            for p in peak_objs_by_chrm[chrm][strand]:
                row = [p.chrm, p.left, p.right, p.strand, p.height]#, p.max_bin]
                if not hasattr(p, 'gene') or p.gene is None:
                    continue
                row += [str(p.gene[key]) for key in dict(p.gene).keys()]
                row += [str(p.pvalues[key]) for key in dict(p.pvalues).keys()]
                #else:
                #    row += ['na' for x in range(1,14)]
                #    row += ['na' for x in range(1,10)]
                peak_list.append(row)
    peak_table = pandas.DataFrame(peak_list, columns=output_cols)
    peak_table = peak_table[peak_table['transcript_id']!='na']
    peak_table.to_csv('peak_table.txt', sep='\t', index=False)
    logger.info(
        '''Resulting pandas.DataFrame of peaks contains %i peaks and %i columns.
        ''' % (
            len(peak_table.index), len(peak_table.columns)
        )
    )
    return peak_table


def read_config(filename):
    """Expect:
experiment_name\tname  # Optional.
clip_replicate\tfilename1.wig
clip_replicate\tfilename2.wig
clip_replicate\tfilename3.wig
gtf_filename\tgtf_filename
rna_seq_filename\tfilename
neg_ip_filename\tfilename
    """
    config = {}
    with open(filename, 'r') as f:
        for li in f:
            li = li.partition('#')[0]  # Skip comments.
            if li == '': continue
            s = li.rstrip('\n').split('\t')
            config.setdefault(s[0], [])
            try: config[s[0]].append(s[1])
            except: print "Error parsing line %s in config file. Wrong length?" % li
    for key in config:
        if len(config[key]) == 1:
            config[key] = config[key][0]
    if 'experiment_name' not in config:
        config['experiment_name'] = os.path.basename(filename)
    return config


def complement(s):
    #basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    basecomplement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def rc(s):
    s = s[::-1]
    s = complement(s)
    return s


def score_metrics(dir_name, config):
    print "score_metrics(): given %s directory." % dir_name
    fname_list = []
    for filename in glob.glob(dir_name + '/*'):
        fname_list.append(filename)
    li = ""
    for filename in sorted(fname_list, key=lambda x: int(os.path.basename(x))):
        print "Scoring metrics for %s" % filename
        li += score_metric(filename)
    with open('score_metrics_%s.txt' % config['experiment_name'], 'w') as f:
        f.write(li)
    return li


def score_metric(filename, label="", given_peaks=False, peaks=False):
    if not label:
        label = os.path.basename(filename)
    if not given_peaks:
        peaks = pandas.DataFrame.from_csv(filename, sep='\t')
    if len(peaks.index) == 0:
        logger.warn('No peaks in file %s.' % filename)
        return "No peaks."
    get_sequences(peaks)
    score_binding_site(peaks)
    #run_dreme(peaks, label)
    positives = score_positives(peaks)
    return write_metrics(peaks, positives, label)


def run_dreme(combined, label):
    if not os.path.exists('combined_fasta'):
        os.system('mkdir combined_fasta')
    fasta_filename = 'combined_fasta/%s.fa' % label
    as_fasta = ""
    for index, peak_row in combined.iterrows():
        start = combined.loc[index, 'left']
        end = combined.loc[index, 'right']
        chrm = combined.loc[index, 'chrm']
        name = combined.loc[index, 'gene_name']
        name += ":%s:%i:%i" % (chrm, start, end)
        as_fasta += ">%s\n%s\n" % (name, combined.loc[index, 'seq'])
    with open(fasta_filename, 'w') as f:
        f.write(as_fasta)
    out_dir = 'dreme_results/%s' % label
    os.system('/home/dp/meme/bin/dreme -e 0.000000001 -p %s -norc -maxk 10 -oc %s' % (fasta_filename, out_dir))


def get_sequences(combined):
    fasta_filename = 'lib/c_elegans.WS235.genomic.fa'
    sequences = dict((re.sub('CHROMOSOME_', '', p.name), p.seq) for p in HTSeq.FastaReader(fasta_filename))
    for index, peak_row in combined.iterrows():
        start = combined.loc[index, 'left']
        end = combined.loc[index, 'right']
        chrm = combined.loc[index, 'chrm']
        seq = sequences[chrm][start:end]
        if combined.loc[index, 'strand'] == '-':
            seq = rc(seq)
        combined.loc[index, 'seq'] = seq


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
           with_fbe=len(peaks[peaks['has_fbe']==1]),
           fbe_perc= float(100 * len(peaks[peaks['has_fbe']==1])/len(peaks)),
           without_fbe=len(peaks[peaks['has_fbe']==0]),
           observed=positives['number of observed positives'],
           expected=positives['expected'],
           missing=positives['missing positives'])
    return li


def score_binding_site(peaks):
    for index, peak_row in peaks.iterrows():
        if re.search('tgt\w\w\wat', peaks.loc[index, 'seq']) is not None:
            peaks.loc[index, 'has_fbe'] = 1
        else:
            peaks.loc[index, 'has_fbe'] = 0


def score_positives(peaks):
    known_pos = set(['gld-1', 'htp-1', 'htp-2', 'mpk-1', 'him-3',
                         'fbf-1', 'lip-1', 'syp-2', 'fbf-2', 'fog-1',
                         'fem-3', 'syp-3', 'gld-3', 'fog-3', 'egl-4'])
    obs_genes = set(peaks['gene_name'])
    obs_pos = known_pos & obs_genes
    missing_pos = known_pos - obs_genes
    obs_pos_n = len(list(obs_pos))
    missing_pos_n = len(list(missing_pos))
    return {'observed positives': obs_pos, 'number of observed positives': obs_pos_n,
            'missing positives': missing_pos, 'number of missing positives': missing_pos_n,
            'expected': len(list(known_pos))}


def load_and_combine_replicates(config):
    combined = {}
    hyp = {}
    rep_dir_names = [
        "%s/peaks/%s/" % (config['experiment_name'], os.path.splitext(os.path.basename(rep_bam_path))[0]) for rep_bam_path in config['clip_replicate']]
    logging.info('After calling peaks and doing statistics, \
    replicates are loaded to be combined, from directory names: %s' % str(rep_dir_names))
    all_hyp_nums = set()
    if 'min_rep_number' in config:
        min_rep_number = int(config['min_rep_number'])
    else:
        min_rep_number = 3
    for rep_name in rep_dir_names:
        for hyp_num in range(1, len(glob.glob(rep_name + '/*')) + 2):
            filename = '%s/null_hyp_%i.txt' % (rep_name, hyp_num)
            if os.path.exists(filename):
                all_hyp_nums |= set([hyp_num])
    for hyp_num in all_hyp_nums:
        hyp[hyp_num] = {}
        for rep_name in rep_dir_names:
            filename = '%s/null_hyp_%i.txt' % (rep_name, hyp_num)
            if not os.path.exists(filename):
                logger.error('Missing peaks filename %s' % filename)
                continue
            hyp[hyp_num][rep_name] = get_peaks(filename)
            li = "Peaks list %s for hyp number %i comprises %i peaks." % (
                filename, hyp_num, len(hyp[hyp_num][rep_name]))
            logging.info(li)
            print li
        # A dict with key=gene is returned.
        combined[hyp_num] = combine_replicates.combine_peaks_not_pandas(
            hyp[hyp_num], min_num_reps=min_rep_number,
            output_dir='data/{v}/'.format(v=config['experiment_name']),
            one_experiment=True)
        logging.info("Combined peaks list for hyp number %i comprises %i peaks." % (
            hyp_num, len(combined[hyp_num])
        ))
        all_rows = []
        for gene in combined[hyp_num]:
            all_rows.extend(combined[hyp_num][gene])
        combined[hyp_num] = pandas.DataFrame(all_rows)
        write_combined(combined[hyp_num], str(hyp_num), config)
    return combined


def get_peaks(filename):
    peaks = pandas.DataFrame.from_csv(filename, sep='\t')
    return peaks


def write_combined(combined, label, config):
    out_dir = config['experiment_name'] + '/peaks/'
    if not os.path.exists(config['experiment_name']):
        os.system('mkdir ' + config['experiment_name'])
    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)
    if not os.path.exists('combined'):
        os.system('mkdir ' + out_dir + 'combined_%s/' % config['experiment_name'])
    with open(out_dir + 'combined_%s/%s' % (config['experiment_name'], label), 'w') as f:
        combined.to_csv(f, sep='\t')


def to_list(pob):
    pl = []
    for p in pob:
        row = p.__dict__
        for key in row['pvalues']:
            row[key] = row['pvalues'][key]
        pl.append(row)
    return pl


def load_tables_of_cwt_peak_calls(config, args):
    logger.info(
        '''CWT peak calls will be made using the configuration data %s.
        If successful, tables will be written to data/cwt_calls/.''' % str(config))
    peak_tables = {}
    for clip_replicate in config['clip_replicate']:
        if not args.just_stats:
            peak_tables[clip_replicate] = call_peaks_by_cwt_on_replicate(
                clip_replicate, config, load_data=not args.overwrite)
        else:
            peak_tables[clip_replicate] = pandas.read_csv(
                'data/cwt_calls/' + os.path.basename(clip_replicate), sep='\t', index=False)
    return peak_tables


def get_gtf(args, config):
    if args.load_gtf:
        with open('lib/gtf_as_dict.p', 'rb') as f:
            as_d = pickle.load(f)
            return as_d
    else:
        gtf_df = pandas.read_csv(config['gtf_filename'], sep='\t')
        gtf_l = add_signal.df_to_dict(gtf_df)
        with open('lib/gtf_as_dict.p', 'wb') as f:
            pickle.dump(gtf_l, file=f)
        return gtf_l


def add_reads_to_peaks_and_return_peak_objects(peak_table, config, gtf_l, clip_replicate):
    if len(peak_table.index) == 0:
        logger.warn('No peaks called in peak table!')
    ga_raw = peak_calling_tools.load_bed_file(config['bed_dir'] + '/' + os.path.basename(clip_replicate).partition('wig')[0] + 'bed')
    pob = add_signal.add_signal_to_replicate(ga_raw, peak_table, gtf_l, clip_replicate, config)
    return pob


def do_stats_apply_fdr_and_output_filtered_files(
        pob, config, clip_replicate_filename, nb_fits=None):
    '''Apply statistical cutoffs and evaluate the result.
    Input:
    pob: Peak objects
    config: config dict
    clip_replicate_filename:
    nb_fits: Known nb fits.
    Output:

    '''
    if (nb_fits is None):
        nb_fits = stats_to_call_peaks.do_statistics(
            pob, config, clip_replicate_filename)
    else:
        stats_to_call_peaks.do_statistics(
            pob, config, clip_replicate_filename,
            nb_fits=nb_fits)
    peak_list = to_list(pob)
    peak_df = pandas.DataFrame(peak_list)
    if not os.path.exists('data/%s/peaks_with_stats' % config['experiment_name']):
        os.system('mkdir data/%s/peaks_with_stats' % config['experiment_name'])
    peak_df.to_csv('data/%s/peaks_with_stats/%s' % (
        config['experiment_name'], os.path.basename(clip_replicate_filename)),
        sep='\t', index=False)
    if len(peak_list) == 0: return None
    print "FDR correction..."
    stats_to_call_peaks.fdr_correction(config, peak_df, clip_replicate_filename)
    # Write peaks and evaluate for each of the seven different null hypothesis.
    stats_to_call_peaks.evaluate_hypothesis(peak_df, clip_replicate_filename, config)
    return nb_fits


def load_peaks_with_stats_and_apply_fdr_and_write(
        clip_replicate_filename, config, alpha=0.01):
    peak_df = pandas.read_csv('data/{exp}/peaks_with_stats/{rep}'.format(
        exp=config['experiment_name'], rep=clip_replicate_filename
    ))
    stats_to_call_peaks.fdr_correction(
        config, peak_df, clip_replicate_filename, alpha=alpha)
    stats_to_call_peaks.evaluate_hypothesis(
        peak_df, clip_replicate_filename, config, alpha=alpha)


def call(args, config, gtf_l, ga_raw, do_start_logger=True):
    print 'calling peaks...'
    if do_start_logger: start_logger(config['experiment_name'])
    logger.info('Experiment bedfiles {d}'.format(
        d=ga_raw.keys()[0]))
    # Outputs pickled CWT call objects to data/raw_* and peak tables to data/cwt_calls.
    # If args.overwrite is False and the .p objects exist, they will be loaded.
    peak_tables = load_tables_of_cwt_peak_calls(config, args)
    logger.info('call(): Finished loading/creating %i cwt tables: %s' % (
        len(peak_tables), str([str(peak_tables[x]) for x in peak_tables])
    ))
    nb_pvals_local, nb_pvals_gene, nb_fits = (None, None, None)
    for clip_replicate_filename in peak_tables:
        peak_table = peak_tables[clip_replicate_filename]
        logger.info('call(): Calling add_signal.add_signal_to_replicate for %s' % clip_replicate_filename)
        pob = add_signal.add_signal_to_replicate(
            ga_raw, peak_table, gtf_l, clip_replicate_filename, config)
        nb_fits = do_stats_apply_fdr_and_output_filtered_files(
            pob, config, clip_replicate_filename, nb_fits=nb_fits)
    load_and_combine_replicates(config)
    score_metrics('%s/peaks/combined_%s/' % (config['experiment_name'], config['experiment_name']), config)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    Call peaks in CLIP-seq data.
    Outputs first pickled objects to data/raw_peak_objs_by_chrm_from_callpeaks_*,
    and equivalent dataframe-based tables to data/cwt_calls/.
    Requires one gtf library file, two bed files of negative ip and rna-seq controls,
    a folder of CLIP bed files, and a folder of CLIP bedgraph files. Paths to these
    files, and an experiment name, must be given in a config file referenced with
    -c <config_file>.
    ''')
    parser.add_argument(
        '-c', '--config', default='clip_config',
    )
    parser.add_argument(
        '-w', '--just_cwt', action='store_true', default=False,
        help='''Only do a CWT call, without statistics (Default: False)'''
    )
    parser.add_argument(
        '-s', '--just_stats', action='store_true', default=False,
        help='Only do a stats call, without CWT (Default: False)'
    )
    parser.add_argument(
        '-x', '--overwrite', action='store_true', default=False,
        help='Overwrite existing peak object pickle objects (Default: False)'
    )
    parser.add_argument(
        '-l', '--load_gtf', action='store_true', default=False,
        help='Load an existing gtf object to speed up (Default: False)'
    )
    args = parser.parse_args()
    config = read_config(args.config)
    # Load high RAM data. Slow.
    print "Loading GTF..."
    gtf_data = get_gtf(args, config)
    ga_raw = {}
    print "Loading bed files..."
    for clip_replicate in config['clip_replicate']:
        ga_raw[clip_replicate] = peak_calling_tools.load_bed_file(
            config['bed_dir'] + '/' + os.path.basename(clip_replicate).partition('wig')[0] + 'bed')
    call(args, config, gtf_data, ga_raw)
