#library(ashr)
library(doMC)
#library(feather)
#library(tidyverse)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(readr)


#vcf_file <- "/gpfs/commons/home/dmeyer/bindingQTL_share/genotype/chrAll_QCFinished_full_2sample.anno.vcf.gz"
#input_files <- c('/gpfs/commons/home/dmeyer/bindingQTL_share/24a-hNIL-control-tdp/allelic/24a-hNIL-control-tdp_input_allelic.out',
#    '/gpfs/commons/home/dmeyer/bindingQTL_share/24a-hNIL-c9-tdp/allelic/24a-hNIL-c9-tdp_input_allelic.out',
#    '/gpfs/commons/home/dmeyer/bindingQTL_share/24a-hNIP-control-tdp/allelic/24a-hNIP-control-tdp_input_allelic.out',
#    '/gpfs/commons/home/dmeyer/bindingQTL_share/24a-hNIP-c9-tdp/allelic/24a-hNIP-c9-tdp_input_allelic.out')
#ip_files <- c('/gpfs/commons/home/dmeyer/bindingQTL_share/24a-hNIL-control-tdp/allelic/24a-hNIL-control-tdp_ip_allelic.out',
#               '/gpfs/commons/home/dmeyer/bindingQTL_share/24a-hNIL-c9-tdp/allelic/24a-hNIL-c9-tdp_ip_allelic.out',
#               '/gpfs/commons/home/dmeyer/bindingQTL_share/24a-hNIP-control-tdp/allelic/24a-hNIP-control-tdp_ip_allelic.out',
#               '/gpfs/commons/home/dmeyer/bindingQTL_share/24a-hNIP-c9-tdp/allelic/24a-hNIP-c9-tdp_ip_allelic.out')

donors <- c('CTRL-NEUHE723FGT-02545-G', 'CASE-NEUFV237VCZ-01369-G', 'CTRL-NEUHE723FGT-02545-G', 'CASE-NEUFV237VCZ-01369-G')
rbps = c("24a-hNIL-control-tdp", "24a-hNIL-c9-tdp", "24a-hNIP-control-tdp", "24a-hNIP-c9-tdp")

input_files <- paste0("data/",rbps, "_input_allelic.out")
ip_files <- paste0("data/",rbps, "_ip_allelic.out")

#vcf_case <- read_feather("/gpfs/commons/home/dmeyer/bindingQTL_share/genotype/CTRL-NEUHE723FGT-02545-G.genotype.feather")
#vcf_ctrl <- read_feather("/gpfs/commons/home/dmeyer/bindingQTL_share/genotype/CASE-NEUFV237VCZ-01369-G.genotype.feather")
#
#vcf_combined <-
#  vcf_case%>%filter(gt %in% c("0/1", "1/0"))%>%
#  inner_join(vcf_ctrl%>%filter(gt %in% c("0/1", "1/0")),
#             by=c('chrom','position', 'variantID', 'refAllele', 'altAllele')
#  )
#
#het_snps <- vcf_combined$variantID
het_snps <- readLines("het_snps.txt")

registerDoMC(4)

# SNP filtering done here and written out to new file at the end
res <- list()
for (i in 1:4) {
    # Used for getting per-donor heterozygous SNPs
    #if (i == 1 || i == 3) {
    #  vcf = vcf_ctrl
    #} else {
    #  vcf = vcf_case 
    #}
    input_file = input_files[i]
    ip_file = ip_files[i]
    sample_index <- donors[i]
    #sum_stats_file <- str_replace(input_file, ".out$", ".sum_stats.tsv")
    input = read_tsv(input_file)%>%rename(chrom=1)
    res[[i]] <- nrow(input)
    input = input[(input$variantID != "."),]
    
    # This is to get heterozygous variants per sample
    #input = input[(input$variantID %in% vcf[vcf$altCount == 1,]$variantID),]
    input = input[(input$variantID %in% het_snps),]
    res[[i]] <- c(res[[i]],  nrow(input))
    
    input = input[input$totalCount >= 10,]
    res[[i]] <- c(res[[i]],  nrow(input))
    
    ip_data = read_tsv(ip_file)%>%rename(chrom=1)%>%
        filter(totalCount >= 30)
    joined=inner_join(input, ip_data, by = join_by(chrom, position, variantID, refAllele, altAllele))
    res[[i]] <- c(res[[i]],  nrow(joined))
    
    input <- input[input$variantID %in% joined$variantID,]
    sum_stats = foreach(i = 1:nrow(input), .combine = bind_rows) %dopar% {
        here = input[i,]
        
        bt = binom.test(c(here$altCount, here$refCount), p = 0.5)
        lor_ci = log(bt$conf.int)
        or_ci = bt$conf.int
        tibble( 
                ci1 = or_ci[1],
                ci2 = or_ci[2],
                or_mean = mean(bt$conf.int),
                or_se = (or_ci[2] - or_ci[1])/4,
                lor_mean = mean(lor_ci),     
                lor_se = (lor_ci[2] - lor_ci[1])/4 )
    }
    
    p = 0.5
    epsilon = 0.3
    #ci = sum_stats$or_mean
    idx = apply(sum_stats, 1, function(x) {
            (x['ci1'] >= p - epsilon) &
            (x['ci1'] <= x['ci2']) & 
            (x['ci2'] <= p + epsilon)
    })
    
    res[[i]] <- c(res[[i]],  sum(idx))
    
    filtered_input_file <- str_replace(input_file, "out$", paste0("10readsInput_filtered_epsilon",epsilon,"_sharedhet.txt"))
    #filtered_input_file <- str_replace(input_file, "out$", paste0("filtered_epsilon",epsilon,".txt"))
    write_tsv(input[idx,], filtered_input_file)
}

