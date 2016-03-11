import HTSeq
import numpy as np
import pickle
import pandas
import re
import os
import glob
import argparse
import logging
import time
import datetime
import collections
import sys
from rc import rc
import config
from get_sequences import get_sequences
from score_metrics import score_metrics

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


def any_have_na(peak_objs):
    has_na = False
    for apeak in peak_objs:
        if apeak.height == 'na':
            logger.warn("Has na height: %s" % str(apeak.__dict__))
    if not has_na:
        logger.info("No peaks have na height.")

import inspect
def call_peaks_by_cwt_on_replicate(clip_wig_filename, config, load_data=False):
    """Calls or loads CWT peaks for a wig filename.
    In: filename of .wig, config data.
    Writes: data/expname/*.p datafiles of CWT peaks.
    Writes: data/expname/cwt_calls/ tables of CWT peaks.
    Returns: pandas.DataFrame of the CWT table written to cwt_calls/
    """
    print inspect.currentframe().f_code.co_name
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
            if peak_table is not None:
                loaded_data = True
                print li
            else:
                print "Failed to load data for %s (got None). \
Regenerating..." % (
                clip_wig_filename)
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
        if peak_table is not None:
            with open(datafile, 'wb') as f:
                print "Writing %i peaks to %s" % (len(peak_table), datafile)
                pickle.dump(peak_objs_by_chrm, f)
        else:
            print "Tried to call peaks again, but got a None value."
    if not os.path.exists('data/%s/cwt_calls/' % config['experiment_name']):
        os.system('mkdir data/%s/cwt_calls' % config['experiment_name'])
    if peak_table  is None:
        print "Got None value for peaks bedgraph file %s / datafile %s" % (
            clip_wig_filename, datafile)
        print "Bedgraph file + has %i lines minus has %i lines, datafile has %i" % (
            sum([1 for x in open(re.sub('\.wig', '_+.wig', clip_wig_filename)).readlines()]),
            sum([1 for x in open(re.sub('\.wig', '_-.wig', clip_wig_filename)).readlines()]),
            sum([1 for x in open(datafile).readlines()]))
        return None
    peak_table.to_csv('data/%s/cwt_calls/%s' % (
        config['experiment_name'], os.path.basename(clip_wig_filename)),
                      sep='\t', index=False)
    return peak_table


def define_output_cols(peak_objs_by_chrm):
    output_cols = ['chrm', 'left', 'right', 'strand',
                   'height']#, 'max_bin']  # Simple colums.
    print inspect.currentframe().f_code.co_name
    for key in peak_objs_by_chrm.keys():
        print key
        print peak_objs_by_chrm[key]
        if ('+' in peak_objs_by_chrm[key]) and (
            len(peak_objs_by_chrm[key]['+']) > 0):
            try:
                output_cols += dict(peak_objs_by_chrm[key]['+'][0].gene).keys()
                output_cols += dict(peak_objs_by_chrm[key]['+'][0].pvalues).keys()
                print "+ used"
                return output_cols
            except: return output_cols
        elif ('-' in peak_objs_by_chrm[key]) and (
            len(peak_objs_by_chrm[key]['-']) > 0):
            #try:
            output_cols += dict(peak_objs_by_chrm[key]['-'][0].gene).keys()
            output_cols += dict(peak_objs_by_chrm[key]['-'][0].pvalues).keys()
            print "- used"
            return output_cols
        print "---"
            ##except: return output_cols
    print "No peaks..."
    return output_cols


def convert_peak_objs_to_table(peak_objs_by_chrm):
    print inspect.currentframe().f_code.co_name
    total_peaks = 0
    for chrm in peak_objs_by_chrm:
        for strand in peak_objs_by_chrm[chrm]:
            total_peaks += len(peak_objs_by_chrm[chrm][strand])
    if total_peaks == 0: return None
    logger.info(
        'Converting %i peak objects in %i chromosomes to a pandas.DataFrame.' % (
            total_peaks, len(peak_objs_by_chrm)))
    # Convert to a table and write.
    peak_list = []
    output_cols = define_output_cols(peak_objs_by_chrm)
    for chrm in peak_objs_by_chrm:
        for strand in peak_objs_by_chrm[chrm]:
            for i, p in enumerate(peak_objs_by_chrm[chrm][strand]):
#                if i < 10: print p
                row = [p.chrm, p.left, p.right, p.strand, p.height]#, p.max_bin]
                if not hasattr(p, 'gene') or p.gene is None:
                    continue
                row += [str(p.gene[key]) for key in dict(p.gene).keys()]
                row += [str(p.pvalues[key]) for key in dict(p.pvalues).keys()]
                if len(row) != len(output_cols):
                    print "row len %i col len %i" % (len(row), len(output_cols))
                    continue
                peak_list.append(row)
    print "Loaded %i peaks" % len(peak_list)
    if len(peak_list) < 1:
        print "Empty peak objects object!"
        return None
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


def load_and_combine_replicates(config):
    print inspect.currentframe().f_code.co_name
    combined = {}
    peaks_by_hypothesis = collections.defaultdict(dict)
    rep_dir_names = [
        "%s/peaks/%s/" % (config['experiment_name'], os.path.splitext(os.path.basename(rep_bam_path))[0]) for rep_bam_path in config['clip_replicate']]
    logging.info('After calling peaks and doing statistics, \
    replicates are loaded to be combined, from directory names:\
    %s' % str(rep_dir_names))
    all_hypothesis = set()
    if 'min_rep_number' in config:
        min_rep_number = int(config['min_rep_number'])
    else:
        min_rep_number = 3
    for rep_name in rep_dir_names:
        for filename in glob.glob(rep_name + '/*.txt'):
            all_hypothesis |= set([os.path.basename(filename)])
    for hypothesis in all_hypothesis:
        for rep_name in rep_dir_names:
            filename = '%s/%s' % (rep_name, hypothesis)
            if not os.path.exists(filename):
                logger.error('Missing peaks filename %s' % filename)
                continue
            peaks_by_hypothesis[hypothesis][
                rep_name] = pandas.read_csv(filename, sep='\t')
            li = "Peaks list %s comprises %i peaks." % (
                filename, len(peaks_by_hypothesis[hypothesis][rep_name]))
            logging.info(li)
            print li
        # A dict with key=gene is returned.
        combined[hypothesis] = combine_replicates.combine_peaks_not_pandas(
            peaks_by_hypothesis[hypothesis], min_num_reps=min_rep_number,
            output_dir='data/{v}/'.format(v=config['experiment_name']),
            one_experiment=True)
        logging.info("Combined peaks list for hypothesis %s comprises %i peaks." % (
            hypothesis, len(combined[hypothesis])
        ))
        all_rows = []
        for gene in combined[hypothesis]:
            all_rows.extend(combined[hypothesis][gene])
        combined[hypothesis] = pandas.DataFrame(all_rows)
        write_combined(
            combined[hypothesis], str(hypothesis), config)
    return combined


#def get_peaks(filename):
#    peaks = pandas.DataFrame.from_csv(filename, sep='\t')
#    return peaks


def write_combined(combined, label, config):
    out_dir = config['experiment_name'] + '/peaks/'
    if not os.path.exists(config['experiment_name']):
        os.system('mkdir ' + config['experiment_name'])
    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)
    out_dir += 'combined_%s/' % config['experiment_name']
    if not os.path.exists(out_dir):
        os.system('mkdir ' + out_dir)
#    if len(combined.index) == 0:
#        return False
    combined.index = range(len(combined.index))
    combined.to_csv(out_dir + label, sep='\t', index=False)


def to_list(pob):
    pl = []
    for p in pob:
        row = p.__dict__
        for key in row['pvalues']:
            row[key] = row['pvalues'][key]
        pl.append(row)
    return pl


def load_tables_of_cwt_peak_calls(config, args):
    print inspect.currentframe().f_code.co_name
    logger.info(
        '''CWT peak calls will be made using the configuration data %s.
        If successful, tables will be written to data/cwt_calls/.''' % str(config))
    peak_tables = {}
    for clip_replicate in config['clip_replicate']:
        clip_replicate = re.sub('_[-+]\.wig', '.wig', clip_replicate)
        if args.just_stats:
            peak_tables[clip_replicate] = pandas.read_csv(
                'data/%s/cwt_calls/' % config['experiment_name'] \
                + os.path.basename(clip_replicate),
                sep='\t')
        else:
            a_table = call_peaks_by_cwt_on_replicate(
                clip_replicate, config, load_data=not args.overwrite)
            if a_table is not None:
                peak_tables[clip_replicate] = a_table
            else:
                print "Got a None value for %s" % clip_replicate
    return peak_tables


def get_gtf(args, config):
    if args.load_gtf:
        with open('lib/gtf_as_dict.p', 'rb') as f:
            as_d = pickle.load(f)
            return as_d
    else:
        gtf_d = add_signal.df_to_dict(
            pandas.read_csv(config['gtf_filename'], sep='\t'))
        with open('lib/gtf_as_dict.p', 'wb') as f:
            pickle.dump(gtf_d, file=f)
        return gtf_d


def add_reads_to_peaks_and_return_peak_objects(peak_table, config, gtf_l, clip_replicate):
    if len(peak_table.index) == 0:
        logger.warn('No peaks called in peak table!')
    ga_raw = peak_calling_tools.load_bed_file(config['bed_dir'] + '/' + os.path.basename(clip_replicate).partition('wig')[0] + 'bed')
    pob = add_signal.add_signal_to_replicate(ga_raw, peak_table, gtf_l, clip_replicate, config)
    return pob


def do_stats_apply_fdr_and_output_filtered_files(
        pob, config, clip_replicate_filename, nb_fits=None,
        skip_nb=False):
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
            pob, config, clip_replicate_filename, skip_nb=skip_nb)
    else:
        stats_to_call_peaks.do_statistics(
            pob, config, clip_replicate_filename,
            nb_fits=nb_fits, skip_nb=skip_nb)
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


def call(args, config, gtf_l, ga_raw, do_start_logger=True,
         skip_nb=False):
    print inspect.currentframe().f_code.co_name
    #args.overwrite = False
    print 'Calling peaks...'
    if do_start_logger: start_logger(config['experiment_name'])
    logger.info('Experiment bedfiles {d}'.format(
        d=ga_raw.keys()[0]))
    #skip = '''
    if args.overwrite: os.system(
        'rm -r data/{a}/cwt_calls/'.format(a=config['experiment_name']))
    if args.overwrite: os.system(
        'rm -r data/{a}/peaks_with_stats/'.format(a=config['experiment_name']))

    # Outputs pickled CWT call objects to data/raw_* and peak tables to data/cwt_calls.
    # If args.overwrite is False and the .p objects exist, they will be loaded.
    peak_tables = load_tables_of_cwt_peak_calls(config, args)
    print "Made peak tables:"
    if peak_tables is None:
        print "Peak tables is None. Ending."
        return
    for k in peak_tables:
        print "\nTable for %s:\n" % k
        if peak_tables[k] is None:
            print "Peak table for %s is None." % k
            continue
        if len(peak_tables[k].index) > 0:
            print peak_tables[k].head(1)
        else:
            print "Peak table for %s is empty." % k
    logger.info('call(): Finished loading/creating %i cwt tables: %s' % (
        len(peak_tables), str([str(peak_tables[x]) for x in peak_tables])
    ))
    nb_pvals_local, nb_pvals_gene, nb_fits = (None, None, None)
    for clip_replicate_filename in peak_tables:  # bedgraph
        print clip_replicate_filename
        if clip_replicate_filename in config['clip_replicate']:
            i = config['clip_replicate'].index(
                clip_replicate_filename)
            print 'in at index %i' % i
            bedname = config['clip_replicate_bed'][i]
            print bedname
        #sys.exit()
        peak_table = peak_tables[clip_replicate_filename]
        logger.info('call(): Calling add_signal.add_signal_to_replicate for %s' % clip_replicate_filename)
        li = ["A\n" + str(x) for x in (
            ga_raw.keys(), peak_table, \
            bedname, config)]
        pob = add_signal.add_signal_to_replicate(
            ga_raw, peak_table, gtf_l, \
            bedname, config)
        nb_fits = do_stats_apply_fdr_and_output_filtered_files(
            pob, config, bedname, nb_fits=nb_fits,
            skip_nb=skip_nb)#'''
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
    parser.add_argument(
        '-m', '--skip_nb', action='store_true', default=False,
        help='Skip the NB calculation for speed (Default: False)'
    )
    args = parser.parse_args()
    lib = config.config(filepath=args.config)
    # Load high RAM data. Slow.
    print "Loading GTF..."
    gtf_data = get_gtf(args, lib)
    ga_raw = {}
    print "Loading bed files..."
    for clip_replicate in lib['clip_replicate']:
        ga_raw[clip_replicate] = peak_calling_tools.load_bed_file(
            lib['bed_dir'] + '/' + os.path.basename(clip_replicate).partition('wig')[0] + 'bed')
    call(args, lib, gtf_data, ga_raw)
