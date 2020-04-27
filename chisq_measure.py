#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jonathantaklee@gmail.com

Test for linkage disequilibrium between all pairs of segregating loci in a cross
"""

import numpy as np
import pandas as pd
from scipy.stats.distributions import chi2


class ChiSquareHMM:
    
    def __init__(self, infile):
        self.infile
        self.df = None
        self.results = None
        self.run()
    
    def unique_sites(self, df):
        """condenses dataframe to independent segregating regions of the genome
        
        input table columns - chr, pos, gp, [HMM calls], freq
        allele calls from HMM should be 0 and 1
        """
        
        df = df[(df['freq'] > .1) & (df['freq'] < .9)].iloc[:, :-1]
        df.index = df['gp']
        df = df.iloc[:, 3:]
        
        cols = df.columns.values
        df_unq = df[cols].loc[(df[cols].shift() != df[cols]).any(axis=1)]
        
        site_calls  = df_unq.T.to_dict(orient='list')
        print('{} unique sites'.format(len(site_calls)))
        return site_calls
    
    def get_pairs(self, sites):
        """constructs a list of pairs of sites"""
        gps = list(sites.keys())
        pairs = []
        for i in range(len(gps)):
            j = i
            while j < len(gps):
                pairs.append((gps[i], gps[j]))
                j += 1
        return pairs
    
    def test_chisq(self, loci, genos1, genos2):
        """compare counts of observed to expected haplotypes at a locus"""
        if loci[0] == loci[1]: return (loci[0], loci[1], 0.0, 1.0)
        
        # count alleles and haplotypes
        ac1 = {0:0, 1:0}
        ac2 = {0:0, 1:0}
        obs = {'00':0, '01':0, '10':0, '11':0}
        
        n_loci = len(genos1)
        for i in range(n_loci):
            allele1 = genos1[i]
            allele2 = genos2[i]
            
            haplo = '{0}{1}'.format(allele1, allele2)
            obs[haplo] += 1
            
            ac1[allele1] += 1
            ac2[allele2] += 1
        
        # observed haplotype counts
        observed = np.array( list(obs.values()) )
        
        # expected haplotype counts
        allele_counts = [ac1[0], ac1[1], ac2[0], ac2[1]]
        af = np.array(allele_counts)/float(n_loci)
        expected = np.array([af[0]*af[2], 
                             af[0]*af[3], 
                             af[1]*af[2], 
                             af[1]*af[3]]) * n_loci
        
        # perform chisquare test
        hap_chisq = ((observed - expected)**2)/expected
        chisq_tot = hap_chisq.sum()
        p = chi2.sf(chisq_tot, 1)
        return (loci[0], loci[1], chisq_tot, p)
    
    def test_loci(self, locus_pairs, genotypes):
        chisq_results = []
        for pair in locus_pairs:
            site1 = genotypes[pair[0]]
            site2 = genotypes[pair[1]]
            test_output = self.test_chisq(pair, site1, site2)
            chisq_results.append(test_output)
        return chisq_results
    
    def run(self):
        df = pd.read_csv(self.infile, sep = '\t')
        sites = self.unique_sites(df)
        site_pairs = self.get_pairs(sites)
        print('running tests...')
        test_results = self.test_loci(site_pairs, sites)
        
        df_pvals = pd.DataFrame(test_results, 
                                columns = ['locus1', 'locus2', 'chisq', 'pval'])
        
        self.df = df
        self.results = df_pvals
        print('done')
    
    def export(self, outfile):
        self.results.to_csv(outfile, sep = '\t', index = False)

