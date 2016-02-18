import HTSeq
import pandas
import os
import sys
import pickle
import logging
import peak
import logging
import time
import collections
_verbose = True

logger = logging.getLogger(__name__)

# For each exon, get the 5' read ends.
def read_region(ga, iv):
    return total_reads_in_bin(ga, iv)


def total_reads_in_bin(ga, iv):
    total_reads = 0
    for iv, value in ga[iv].steps():
        total_reads += value
    return total_reads


def split_into_bins(ga, iv):
    bin_size = 50
    bins_list = []
    for bin_left in range(iv.start, iv.end, bin_size):
        bin_right = bin_left + bin_size
        iv = HTSeq.GenomicInterval(iv.chrom, bin_left, bin_right, iv.strand)
        bins_list.append(total_reads_in_bin(ga, iv))
    return bins_list


def add_local_signal(p, bedfile_cat_to_ga_key, ga):
    width = 250
    for bedfile in bedfile_cat_to_ga_key:  # Values are filenames.
        bedfile = os.path.basename(bedfile)
        # Determine local signal.
        iv = HTSeq.GenomicInterval(
            p.chrm,
            max(0, p.left-width),
            p.right+width,
            p.strand)  # Chr, start, end, strand.
        bins_list = split_into_bins(ga[bedfile_cat_to_ga_key[bedfile]], iv)
        p.local[bedfile] = bins_list


def add_gene_signal(peak, gtf, bedfile_cat_to_ga_key, ga):
    genes = {}
    if (peak.gene_name is None): return False
    if (peak.transcript_id not in gtf):
        li = "\tNo ID in gtf for: %s.." % str(peak.transcript_id)
        logger.warn(li)
        return False
    txpt_id = peak.transcript_id
    _exons = [x for x in gtf[txpt_id] if x['2']=='exon']
    _exons.sort(key=lambda x: int(x['3']))
#    exons_df = gtf_df[(gtf_df['transcript_id']==txpt_id) & (gtf_df['2']=='exon')]
    exons = [HTSeq.GenomicInterval(exon['0'], exon['3'], exon['4'], exon['6']) for exon in _exons]
    for bedfile in bedfile_cat_to_ga_key:  # Values are filenames.
        if txpt_id in genes and bedfile in genes[txpt_id]:
            peak.exons[bedfile] = genes[txpt_id][bedfile]
        else:
            # Determine gene signal.
            peak.exons[bedfile] = []
            genes.setdefault(txpt_id, {})
            genes[txpt_id][bedfile] = []
            for exon in exons:  # Each element is an interval tuple.
                bins_list = split_into_bins(ga[bedfile_cat_to_ga_key[bedfile]], exon)
#                reads = read_region(ga[bedfile_cat_to_ga_key[bedfile]], exon)
                genes[txpt_id][bedfile].extend(bins_list)
                peak.exons[bedfile].extend(bins_list)
    return True


def df_to_dict(gtf):
    as_d = collections.defaultdict(list)
    gtf_d = gtf.to_dict('records')
    for row in gtf_d:
        as_d[row['transcript_id']].append(row)
    return as_d


def add_signal(ga_raw, peak_table, gtf, config):
    if type(gtf) == pandas.DataFrame:
        gtf = df_to_dict(gtf)
    for exp in peak_table:
        add_signal_to_replicate(ga_raw, peak_table[exp], gtf, exp, config)


def add_signal_to_replicate(ga, peak_table, gtf, clip_replicate, config):
    li = "\tadd_signal.add_signal_to_replicate() called on %s..." % clip_replicate
    logger.info(li)
    pob = []
    bedfile_cat_to_ga_key = {
            'clip': clip_replicate,
            'rna_seq': 'rna_seq',
            'neg_ip': 'neg_ip'}
    for index, _peak in peak_table.iterrows():
        #if _peak.gene_name != 'gld-1': continue
        iv = (_peak['chrm'], _peak['left'], _peak['right'], _peak['strand'])
        pob.append(peak.peak(*iv))
        pob[-1].gene_name = _peak['gene_name']
        pob[-1].transcript_id = _peak['transcript_id']
        pob[-1].transcript_name = _peak['transcript_name']
        add_gene_signal(pob[-1], gtf, bedfile_cat_to_ga_key, ga)
        add_local_signal(pob[-1], bedfile_cat_to_ga_key, ga)
    logger.info("Adding gene signal ...")
    return pob

