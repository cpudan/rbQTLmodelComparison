# %%
import torch
import pyro

#import seaborn as sns
#import matplotlib.pyplot as plt
#import scipy.stats
import pandas as pd

import numpy as np

from pooledQTL import beta_model
from pooledQTL import interval_utils

from pathlib import Path, PosixPath
#import plotnine as p9

fdr_threshold = 0.05
device = "cuda:0" if torch.cuda.is_available() else "cpu"
use_structured_guide = True # performs much better

# %% [markdown]
# # Load input and IP RNA allelic counts

# %%
rbps = ["24a-hNIL-control-tdp", "24a-hNIL-c9-tdp", "24a-hNIP-control-tdp", "24a-hNIP-c9-tdp"]
donors = [
    "CTRL-NEUHE723FGT-02545-G", "CASE-NEUFV237VCZ-01369-G",   
    "CTRL-NEUHE723FGT-02545-G", "CASE-NEUFV237VCZ-01369-G"]   
vcf = "/gpfs/commons/groups/knowles_lab/bindingQTL_share/genotype/chrAll_QCFinished_full_2sample.anno.vcf"
#vcf = "/gpfs/commons/groups/knowles_lab/bindingQTL_share/genotype/chrAll_QCFinished_full.anno.vcf"
basedir = Path("../data")

# %%
input_files = [basedir / f"{rbp}_input_allelic.out" for rbp in rbps]
#input_files_filtered = [basedir / rbp / f"allelic/{rbp}_input_allelic.filtered_epsilon0.3.txt" for rbp in rbps]
input_files_filtered = [basedir / f"{rbp}_input_allelic.10readsInput_filtered_epsilon0.3_sharedhet.txt" for rbp in rbps]
IP_files = [basedir / f"{rbp}_ip_allelic.out" for rbp in rbps]
input_counts = [ pd.read_csv(f, sep = "\t", usecols = range(8), index_col = False) for f in input_files]
input_counts_filtered = [ pd.read_csv(f, sep = "\t", usecols = range(8), index_col = False) for f in input_files_filtered]
IP_counts = [ pd.read_csv(f, sep = "\t", usecols = range(8), index_col = False) for f in IP_files]

# %% [markdown]
# # Run Beta Model

# %%
i=0
dat = input_counts_filtered[i]
dat_IP = IP_counts[i]
sample_ind = donors[i]

dat = dat.rename(columns={ "chrom":"contig" })
dat = dat[dat.variantID != "."]
dat = dat[~dat.variantID.duplicated()]
dat_IP = dat_IP[dat_IP.variantID != "."]
dat_IP = dat_IP[~dat_IP.variantID.duplicated()]

# %%
dat

# %%
for i in range(len(rbps)):
    dat = input_counts_filtered[i]
    dat_IP = IP_counts[i]
    sample_ind = donors[i]
    
    dat = dat.rename(columns={ "chrom":"contig" })
    dat = dat[dat.variantID != "."]
    dat = dat[~dat.variantID.duplicated()]
    dat_IP = dat_IP[dat_IP.variantID != "."]
    dat_IP = dat_IP[~dat_IP.variantID.duplicated()]
    
    #imp_merged = sanger.merge(
    #    dat, 
    #    on = ["contig", "position", "variantID", "refAllele", "altAllele"]#, suffixes = ["_hg19",""]
    #) # sanger is GRCh38
    # there are only 0.08% flipped alleles so not worth doing.
    # np.isnan(imp_merged.iloc[:,5:16]).any() all False
    dat["input_ratio"] = dat.altCount / dat.totalCount
    #dat = dat[~dat[sample_ind].isna()]
    dat["pred_ratio"] = 0.5 #dat.loc[:,sample_ind].to_numpy()
    merged = dat.merge(
        dat_IP, 
        on = ("contig", "position", "variantID", "refAllele", "altAllele"), 
        suffixes = ("_input", "_IP"))
    #merged = merged.drop(labels=["contig_y", "position_x" ], axis=1)
    
    merged["IP_ratio"] = merged.altCount_IP / merged.totalCount_IP
    
    print(merged.shape)
    
    # filter for reasonable counts -- previous method
    #input_total_min = 10
    #allele_count_min = 4 
    #ip_total_min = 30
    #dat_sub = merged[merged.totalCount_input >= input_total_min]
    #dat_sub = dat_sub[dat_sub.refCount_input >= allele_count_min]
    #dat_sub = dat_sub[dat_sub.altCount_input >= allele_count_min]
    #dat_sub = dat_sub[dat_sub.totalCount_IP >= ip_total_min]
    #dat_sub = dat_sub[dat_sub.pred_ratio >= 0.45]
    #dat_sub = dat_sub[dat_sub.pred_ratio <= 0.55]
    
    ip_total_min = 30
    dat_sub = merged[merged.totalCount_IP >= ip_total_min]
    print(dat_sub.shape)
    
    epsilon = 0.3
    results_dir = PosixPath(f'../results')
    results_file = results_dir / (f"{rbps[i]}_beta_10readsInput_filtered_epsilon0.3_ipreads30_allhet_struct.tsv.gz")
    results = beta_model.fit_and_save(dat_sub, 
                                     results_file,
                                     use_structured_guide = True,
                                     iterations = 1000,
                                     device = device)

# %% [markdown]
# ## Add in annotations

# %%
exons = pd.read_csv("../data/gencode.v38.exons.txt.gz", 
                   sep = "\t",  
                   index_col = False, usecols = range(3)).rename(columns = {"chr" : "chrom"})
exons = exons[(exons.end - exons.start) >= 9] # remove super short exons
exons_tree = interval_utils.to_interval_trees(exons, exons.chrom.unique())
genc = pd.read_csv("../data/gencode.v38.annotation.gtf.gz",
                   sep = "\t",  index_col = False, header=None, comment="#", usecols = range(5))
genc = genc.rename(columns={k:x for (k,x) in enumerate(["chrom", "gene", "annotation", "start", "end"])})

genc = genc[(genc.end - genc.start) >= 9] # remove super short annotations
def f(x):
    return interval_utils.to_interval_trees(x, genc.chrom.unique())
gene_tree = f(genc[genc.annotation == "gene"][genc.columns[[0,3,4]]])
transcript_tree = f(genc[genc.annotation == "transcript"][genc.columns[[0,3,4]]])
utr_tree = f(genc[genc.annotation == "UTR"][genc.columns[[0,3,4]]])
def get_overlap_window(peaks, snps, window=100000): 
    return [min(len(peaks[row.chrom][row.position-window] | peaks[row.chrom][row.position] | peaks[row.chrom][row.position+window]), 1) for row in snps.itertuples() ]
for i in range(len(rbps)):
    rbp=rbps[i]
    results_dir = PosixPath(f'../results')
    results_file = results_dir / (f"{rbp}_beta_10readsInput_filtered_epsilon0.3_ipreads30_allhet_struct.tsv.gz")
    results = pd.read_csv(results_file, sep="\t")
    
    peaks_file_pos = f"../data/{rbp}-pos_peaks.narrowPeak"
    peaks_file_neg = f"../data/{rbp}-neg_peaks.narrowPeak"
    peaks_pos = pd.read_csv(peaks_file_pos, 
                        header = None,
                        names= ["chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"],
                        sep = "\t",  
                        index_col = False)
    peaks_neg = pd.read_csv(peaks_file_neg, 
                        header = None,
                        names= ["chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"],
                        sep = "\t",  
                        index_col = False)
    peak_tree_pos = interval_utils.to_interval_trees(peaks_pos, peaks_pos.chrom.unique())
    peak_tree_neg = interval_utils.to_interval_trees(peaks_neg, peaks_neg.chrom.unique())
    res2 = results.rename(columns={"contig": "chrom"})

    results["in_peak_pos"] = interval_utils.get_overlap(peak_tree_pos, res2)
    results["in_peak_neg"] = interval_utils.get_overlap(peak_tree_neg, res2)
    results["in_peak"] = results["in_peak_neg"] | results["in_peak_neg"]
    results["near_peak_100k"] = np.array(get_overlap_window(peak_tree_pos, res2)) |  np.array(get_overlap_window(peak_tree_neg, res2))
    results["in_exon"] = interval_utils.get_overlap(exons_tree, results.rename(columns={"contig": "chrom"}))
    results["in_transcript"] = interval_utils.get_overlap(transcript_tree, results.rename(columns={"contig": "chrom"}))
    results["in_gene"] = interval_utils.get_overlap(gene_tree, results.rename(columns={"contig": "chrom"}))
    results["in_utr"] = interval_utils.get_overlap(utr_tree, results.rename(columns={"contig": "chrom"}))
    
    results_file_annotated = results_dir / (f"{rbps[i]}_beta_10readsInput_filtered_epsilon0.3_ipreads30_allhet_struct_with_peaks.tsv.gz")
    results.rename(columns={"contig": "chrom"}).to_csv(results_file_annotated, sep="\t", index=False)


