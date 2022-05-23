
import pandas as pd
from sklearn.linear_model import LinearRegression

import numpy as np

import scipy.stats

from . import io_utils, pyro_utils

import torch
import matplotlib.pyplot as plt

torch_matmul = lambda x,y : (torch.tensor(x) @ torch.tensor(y)).numpy() # do we need this? apparently yes!?

def deconvolve(geno, dat, sample_inds = range(5,16), total_thres = 100, plot = True):
    
    # join genotype data and input allele counts
    merged = geno.merge(dat, on = ["variantID", "refAllele", "altAllele"])

    # consider different defitions of ref vs alt
    geno_flip = geno.rename(columns={"altAllele" : "refAllele", "refAllele":"altAllele"})
    geno_flip.iloc[:,sample_inds] = 1. - geno_flip.iloc[:,sample_inds]
    merged_flip = geno_flip.merge(dat, on = ["variantID", "refAllele", "altAllele"])
    combined = pd.concat((merged,merged_flip), axis=0) # this handles the misordering of alt/ref correctly
    
    # remove any rows with missigness genotypes
    to_keep = np.isnan(combined.iloc[:,sample_inds]).mean(1) == 0. # keep 96%
    combined = combined[to_keep].copy()

    combined["allelic_ratio"] = combined.altCount / combined.totalCount
    
    # only perform deconv using SNPs with >total_thres total counts
    comb_sub = combined[combined.totalCount >= total_thres].copy()

    X = comb_sub.iloc[:,sample_inds].to_numpy() # dosage matrix
    y = comb_sub.allelic_ratio.to_numpy() # observed allelic proportions

    reg_nnls = LinearRegression(positive=True, fit_intercept=False)
    reg_nnls.fit(X, y)
    w = reg_nnls.coef_

    if plot: 
        print("sum(w)=%f ideally would be 1" % w.sum())
        combined["pred"] = torch_matmul( combined.iloc[:,5:16].to_numpy(), w )
        #combined["pred"] = combined.iloc[:,sample_inds].to_numpy() @  w 

        plt.bar(x = range(len(w)), height=w)
        plt.show()

        combined_30 = combined[combined.totalCount >= 30]
        corr,_ = scipy.stats.pearsonr(combined_30.pred, combined_30.allelic_ratio)
        R2 = corr*corr

        plt.scatter(combined_30.pred, combined_30.allelic_ratio, alpha = 0.05)
        plt.title("R2=%.3f" % R2)
        plt.xlabel("Predicted from genotype")
        plt.ylabel("Observed in input")
        plt.show() 
        
    return w

def merge_geno_and_counts(sanger, 
                          dat, 
                          dat_IP, 
                          w, 
                          sample_inds = range(5,16),
                          num_haploids = 18,
                          input_total_min = 10, 
                          allele_count_min = 4, 
                          ip_total_min = 30,
                          plot = True):
    # have to match on rsID because sanger.vcf is hg19 and allelic counts are on hg38
    imp_merged = sanger.rename(columns = {"SNP" : "variantID"}).merge(dat, on = ["variantID", "refAllele", "altAllele"]) # sanger is hg19
    # there are only 0.08% flipped alleles so not worth doing.
    # np.isnan(imp_merged.iloc[:,5:16]).any() all False
    imp_merged["allelic_ratio"] = imp_merged.altCount / imp_merged.totalCount
    X = 0.5 * imp_merged.iloc[:,sample_inds].to_numpy().copy()
    # p = X @ w # WTF doesn't this work!? 
    # p = np.dot(X,w) # doesn't work either
    #p_ = np.array([ X[i] @ w for i in range(X.shape[0]) ])
    imp_merged["pred"] = torch_matmul(X, w)
    
    if plot:
        imp_merged_30 = imp_merged[imp_merged.totalCount >= 30]
        corr,_ = scipy.stats.pearsonr(imp_merged_30.pred, imp_merged_30.allelic_ratio)
        R2 = corr*corr
        plt.scatter(imp_merged_30.pred, imp_merged_30.allelic_ratio, alpha = 0.005) 
        plt.title("R2=%.3f" % R2)
        plt.xlabel("Predicted from genotype")
        plt.ylabel("Observed in input")
        plt.show()

    # merge (imp_geno+input) with IP
    merged = imp_merged.drop(labels=imp_merged.columns[range(5,16)], axis=1).rename(columns={"position_y":"position","contig_x":"contig"}).merge(dat_IP, on = ("contig", "position", "variantID", "refAllele", "altAllele"), suffixes = ("_input", "_IP"))
    #merged = merged.drop(labels=["contig_y", "position_x" ], axis=1)

    dat_sub = merged[merged.totalCount_input >= input_total_min].rename(columns = {"pred" : "pred_ratio"})
    dat_sub = dat_sub[dat_sub.refCount_input >= allele_count_min]
    dat_sub = dat_sub[dat_sub.altCount_input >= allele_count_min]
    dat_sub = dat_sub[dat_sub.totalCount_IP >= ip_total_min]
    dat_sub = dat_sub[dat_sub.pred_ratio >= 0.5/num_haploids]
    dat_sub = dat_sub[dat_sub.pred_ratio <= (1.-0.5/num_haploids)]

    return merged,dat_sub
