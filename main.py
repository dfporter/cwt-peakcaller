import sys
import callpeaks
import argparse
import peak_calling_tools
import pickle
import os
import traceback
import peak_calling_tools
import combine_replicates
import stats_to_call_peaks
import add_signal
import peak
import NBin


def load_args():
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
        '-n', '--no_ui', action='store_true', default=False,
        help='No user input (Default: False)'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = load_args()
    config = callpeaks.read_config(args.config)
    callpeaks.start_logger(config['experiment_name'])
    print "Loading gtf..."
    gtf_data = callpeaks.get_gtf(args, config)
    ga_raw = {}
    print "Loading bed files..."
    ga_raw['neg_ip'] = peak_calling_tools.load_bed_file(config['neg_ip_filename'])
    ga_raw['rna_seq'] = peak_calling_tools.load_bed_file(config['rna_seq_filename'])
    for clip_replicate in config['clip_replicate']:
        ga_raw[clip_replicate] = peak_calling_tools.load_bed_file(
            config['bed_dir'] + '/' + os.path.basename(clip_replicate).partition('wig')[0] + 'bed')
    if args.no_ui:
        callpeaks.call(args, config, gtf_data, ga_raw)
        print "Finished calling peaks. Exiting."
        sys.exit()
    while True:
        try:
            callpeaks.call(args, config, gtf_data, ga_raw)
        except:
            print "Error in callpeaks()."
            print traceback.format_exc()
        print "Hit enter to reload callpeaks(), CTRL-C to close."
        sys.stdin.readline()
        reloaded = False
        while not reloaded:
            try:
                reload(callpeaks)
                reload(peak_calling_tools)
                reload(combine_replicates)
                reload(stats_to_call_peaks)
                reload(add_signal)
                reload(peak)
                reload(NBin)
                print "Successfully recompiled callpeaks et al."
                reloaded = True
            except:
                print "Crash when trying to compile callpeaks."
                print traceback.format_exc()
            print "Hit enter to re-run the script, CTRL-C to close."
            sys.stdin.readline()

