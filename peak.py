import HTSeq
from scipy import signal
import numpy as np
import re
import os
import scipy as sp
import logging
import NBin
import peak

logger = logging.getLogger(__name__)

class peak:

    def __init__(self, chrm, left, right, strand):
        self.chrm = str(chrm)
        self.left = int(left)
        self.right = int(right)
        self.strand = str(strand)
        self.iv = HTSeq.GenomicInterval(chrm, left, right, strand)
        self.local = {}
        self.exons = {}
        self.pvalues = {}
        self.height = 0

    def __repr__(self):
        #return "__dict__:%s" % str(self.__dict__)
        if hasattr(self, 'height'):
            return "%s:%i-%i (height=%i)" % (self.chrm, self.left, self.right, self.height)
        return "%s:%i-%i" % (self.chrm, self.left, self.right)

    def find_max_bin(self):
        """The peak bin is the maximum bin in the center 40% of bins.
        40% * 1kbp = 400 bp = 8 (50 bp) bins.
        """
        x = self.local['clip']
        midrange = x[int(len(x)*0.2):int(len(x)*.2 + int(len(x)*.5))]
        self.max_bin = max(midrange)

    def local_poisson(self, bamfiles, norm_factors):
        mean_local = {}
        for bamfile in bamfiles:
            if len(self.local[bamfile]) == 0:
                self.pvalues['%s_local_poisson' % bamfile] = 1.0
                continue
            mean_local[bamfile] = float(sum(self.local[bamfile]))/float(
                len(self.local[bamfile]))
            mean_local[bamfile] *= norm_factors[bamfile]
            self.pvalues['%s_local_poisson' % bamfile] = sp.stats.poisson.sf(
                self.max_bin, mean_local[bamfile])
            if np.isnan(self.pvalues['%s_local_poisson' % bamfile]):
                if mean_local[bamfile] == 0.0 and self.max_bin > 0:
                    self.pvalues['%s_local_poisson' % bamfile] = 0
                if mean_local[bamfile] == 0.0 and self.max_bin == 0:
                    self.pvalues['%s_local_poisson' % bamfile] = 1.0
            # Bonferroni to make p value the significance for local region.
            self.pvalues['%s_local_poisson' % bamfile] *= len(self.local[bamfile])
            self.pvalues['%s_local_poisson' % bamfile] = min(1., self.pvalues['%s_local_poisson' % bamfile])
            if np.isnan(self.pvalues['%s_local_poisson' % bamfile]):
                if self.max_bin == 0:
                    self.pvalues['%s_local_poisson' % bamfile] = 1.
                else:
                    self.pvalues['%s_local_poisson' % bamfile] = 0.
                logger.error("Error in local_poisson(): %s" % str(self.__dict__))

    def local_norm(self, bamfiles, norm_factors):
        mean_local = {}
        for bamfile in bamfiles:
            norm_list = [norm_factors[bamfile] * x for x in self.local[bamfile]]
            self.local_norm_mu, self.local_norm_std = sp.stats.norm.fit(norm_list)
            self.pvalues['%s_local_norm' % bamfile] = 1 - sp.stats.norm.cdf(
                self.max_bin, self.local_norm_mu, self.local_norm_std)
            if np.isnan(self.pvalues['%s_local_norm' % bamfile]):
                if self.local_norm_mu == 0 and self.max_bin > 0:
                    self.pvalues['%s_local_norm' % bamfile] = 0
                elif self.local_norm_mu == 0 and self.max_bin == 0:
                    self.pvalues['%s_local_norm' % bamfile] = 1.0
            # Bonferroni to make p value the significance for the gene.
            self.pvalues['%s_local_norm' % bamfile] *= len(self.local[bamfile])
            self.pvalues['%s_local_norm' % bamfile] = min(1., self.pvalues['%s_local_norm' % bamfile])

    def close_enough(self, bedfile, known):
        try:
            borders_list = known[self.gene_name]['Locals'][bedfile]
        except:
            return False
        for borders in borders_list:
            if (abs(borders['left'] - self.left) < 100) and (
                    abs(borders['right'] - self.right) < 100):
                logger.info('Could recover local nb signal for bedfile {b} with parameters {p}'.format(
                    b=bedfile, p=borders['params']
                ))
                return borders['params']
        return False

    def local_nb(self, bedfiles, norm_factors, known=None):
        nb_results = {}
        for bedfile in bedfiles:
            if bedfile == 'clip' or len(self.local[bedfile]) == 0:
                self.pvalues['%s_local_nb'% bedfile] = 1.
                nb_results[bedfile] = 1.
                continue
            norm_list = [norm_factors[bedfile] * x for x in self.local[bedfile]]
            params = self.close_enough(bedfile, known)
            if type(params) == bool and (params==False):
                params = NBin.fit_params(norm_list, just_positives=True)
                known[self.gene_name]['Locals'].setdefault(bedfile, [])
                known[self.gene_name]['Locals'][bedfile].append(
                    {'left': self.left, 'right': self.right, 'params': params}
                )
            if params is None:
                pval = np.nan
            else:
                frac_with_signal = \
                    float(len([x for x in norm_list if x>0]))/float(len(norm_list))
                pval = sp.stats.nbinom.cdf(self.max_bin*100, params[0], params[1])
                pval = frac_with_signal * (1. - pval)
            self.pvalues['%s_local_nb' % bedfile] = pval
            if np.isnan(pval):
                if (np.max(self.local[bedfile]) <= 0) and (self.max_bin > 0):
                    self.pvalues['%s_local_nb' % bedfile] = 0.
                elif (np.max(self.local[bedfile]) <= 0) and (self.max_bin == 0):
                    self.pvalues['%s_local_nb'% bedfile] = 1.
                else:
                    self.pvalues['%s_local_nb'% bedfile] = 1.
#            self.pvalues['%s_local_nb' % bedfile] *= len(self.local[bedfile])
            self.pvalues['%s_local_nb' % bedfile] = min(1., self.pvalues['%s_local_nb' % bedfile])
            nb_results[bedfile] = min(1., self.pvalues['%s_local_nb' % bedfile])
        return nb_results

    def gene_poisson(self, bamfiles, norm_factors):
        mean_gene = {}
        for bamfile in bamfiles:
            if not hasattr(self, 'exons') or (bamfile not in self.exons) or (
                len(self.exons[bamfile]) == 0
            ):
                if not hasattr(self, 'exons'):
                    logger.warn("No exons for %s." % str(self.__dict__))
                elif bamfile not in self.exons:
                    logger.warn(
                        "No bamfile %s in self.exons %s" % (bamfile, str(self.exons)))
                self.pvalues['%s_gene_poisson' % bamfile] = 1.
                continue
            mean_gene[bamfile] = float(sum(self.exons[bamfile]))/float(
                len(self.exons[bamfile]))
            mean_gene[bamfile] *= norm_factors[bamfile]
            self.pvalues['%s_gene_poisson' % bamfile] = sp.stats.poisson.sf(
                self.max_bin, mean_gene[bamfile])
            if np.isnan(self.pvalues['%s_gene_poisson' % bamfile]):
                if mean_gene[bamfile] == 0.0 and self.max_bin > 0:
                    self.pvalues['%s_gene_poisson' % bamfile] = 0
                elif mean_gene[bamfile] == 0.0 and self.max_bin == 0:
                    self.pvalues['%s_gene_poisson' % bamfile] = 1.0
            # Bonferroni to make p value the significance for the gene.
            self.pvalues['%s_gene_poisson' % bamfile] *= len(self.exons[bamfile])
            self.pvalues['%s_gene_poisson' % bamfile] = min(1., self.pvalues['%s_gene_poisson' % bamfile])

    def gene_norm(self, bamfiles, norm_factors):
        for bamfile in bamfiles:
            if not hasattr(self, 'exons') or (bamfile not in self.exons):
                self.gene_norm_mu, self.gene_norm_std = (np.nan, np.nan)
                self.pvalues['%s_gene_norm' % bamfile] = np.nan
                continue
            norm_list = [norm_factors[bamfile] * x for x in self.exons[bamfile]]
            self.gene_norm_mu, self.gene_norm_std = sp.stats.norm.fit(norm_list)
            self.pvalues['%s_gene_norm' % bamfile] = 1 - sp.stats.norm.cdf(
                self.max_bin, self.gene_norm_mu, self.gene_norm_std)
            if np.isnan(self.pvalues['%s_gene_norm' % bamfile]):
                if self.gene_norm_mu == 0 and self.max_bin > 0:
                    self.pvalues['%s_gene_norm' % bamfile] = 0
                elif self.gene_norm_mu == 0 and self.max_bin == 0:
                    self.pvalues['%s_gene_norm' % bamfile] = 1.0
            # Bonferroni to make p value the significance for the gene.
            self.pvalues['%s_gene_norm' % bamfile] *= len(self.exons[bamfile])
            self.pvalues['%s_gene_norm' % bamfile] = min(1., self.pvalues['%s_gene_norm' % bamfile])

    def gene_nb(self, bedfiles, norm_factors, known=None):
        nb_results = {}
        for bedfile in bedfiles:
            if bedfile == 'clip':
                self.pvalues['%s_gene_nb' % bedfile] = 1.
                nb_results[bedfile] = 1.
                continue
            if not hasattr(self, 'exons') or bedfile not in self.exons:
                self.gene_norm_mu, self.gene_norm_std = (np.nan, np.nan)
                self.pvalues['%s_gene_norm' % bedfile] = np.nan
                continue
            norm_list = [norm_factors[bedfile] * x for x in self.exons[bedfile]]
            try:
                params = known[self.gene_name]['Exons'][bedfile]
            except:
                params = NBin.fit_params(norm_list, just_positives=True)
                if self.gene_name not in known: known.setdefault(self.gene_name, {})
                if 'Exons' not in known[self.gene_name]: known[self.gene_name].setdefault('Exons', {})
                if bedfile not in known[self.gene_name]['Exons']:
                    known[self.gene_name]['Exons'][bedfile] = params
            if params is None:
                pval = np.nan
            else:
                frac_with_signal = \
                    float(len([x for x in norm_list if x>0]))/float(len(norm_list))
                pval = sp.stats.nbinom.cdf(self.max_bin*100, params[0], params[1])
                pval = frac_with_signal * (1. - pval)
            self.pvalues['%s_gene_nb' % bedfile] = pval
            if np.isnan(pval):
                if (np.max(self.local[bedfile]) <= 0) and (self.max_bin > 0):
                    self.pvalues['%s_gene_nb' % bedfile] = 0.
                elif (np.max(self.local[bedfile]) <= 0) and (self.max_bin == 0):
                    self.pvalues['%s_gene_nb'% bedfile] = 1.
                else:
                    self.pvalues['%s_gene_nb' % bedfile] = 1.
            self.pvalues['%s_gene_nb' % bedfile] = min(1., self.pvalues['%s_gene_nb' % bedfile])
            nb_results[bedfile] = min([1., self.pvalues['%s_gene_nb' % bedfile]])
        return nb_results
