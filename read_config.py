import collections
import os
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
