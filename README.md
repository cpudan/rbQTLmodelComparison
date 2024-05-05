# rbQTL Model Comparison

Code to detect allele specific binding (ASB) of RNA binding proteins (RBP) to RNA from pooled (from multiple individuals) RNA
immunoprecipation & sequencing (RIP-seq) data.

## Setup and file organization

Analysis pipeline is in /src. Underlying code is in the package /src/poolQTL.

For code underlying beta model see /src/pooledQTL/beta_model.py

Counts data can be found at https://drive.google.com/drive/folders/1YgU54SKnJL-zgb5c5ehvrFBPJDPy4-Pt?usp=drive_link.
Code works on the assumption that data is downloaded files are all extracted into a subdirectory of the github repo called /data

You will also need to download two gencode v38 files for generating annotations. These files can also be found in the google drive link and should be put in the /data directory.

Results generated by each model can be found in the /results directory

Due to restrictions on human subjects data, WGS information cannot be shared for the cell lines used to generate this data, however, a list of heterozygous SNPs is included in /het_snps.txt.

## Dependencies

Python package requirements listed in setup.cfg and can be installed by running:
```
pip install plotnine intervaltree scikit-learn torch seaborn pyro-ppl
```

R packages can be installed by running:
```
install.packages(c("doMC", "dplyr", "stringr", "tidyr", "ggplot2", "readr"))
```
Further requirements may need to be installed in on your machine depending on your development environment.
Please check logs of R for errors in setup to guide you through installing the correct packages on your system.
I have tested on CentOS 7.9.2009 and Fedora 39.

## Usage

1. Install any dependencies
2. `git clone https://github.com/cpudan/rbQTLmodelComparison.git`
3. Change directory to git repo `cd rbQTLmodelComparison`
4. Download data from [Google Drive](https://drive.google.com/drive/folders/1YgU54SKnJL-zgb5c5ehvrFBPJDPy4-Pt?usp=drive_link)
5. Move zip file downloaded from google drive to current directory: `mv ~/Downloads/rbQTL_project_data*.zip .`
6. Extract zip and tidy up `unzip rbQTL_project_data*.zip; mv rbQTL_project_data data; rm rbQTL_project_data*.zip`
7. Run the R script `Rscript src/snp_filtering.R` to do preprocessing / SNP filtering
8. Run the python script `python3 src/Beta_model_ASB_TDP43.py` to generate results from beta model
9. Run Chi-squared model and generate plots in `model_comparison.Rmd`
