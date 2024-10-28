# TIRTL
Code for paired TCR data analysis from TIRTL-seq experiments.

## Installation
### Requirements

For GPU functionality (Nvidia only at the moment) you will need to have CUDA installed. You can install it from [here](https://developer.nvidia.com/cuda-downloads).

R dependencies are _Matrix_ and _data.table_ packages, available by default with most distributions. If not, you can install them by running from R console:
```R
install.packages("Matrix")
install.packages("data.table")
```

Python3 dependencies are _numpy_, _pandas_, _cupy_ (for GPU runs, optional). You can install them by running:
```bash
pip install numpy pandas
pip install cupy
```
## Quick Start

For the basic usage of pairing pipeline, you can run the following commands in R:
```R
source("TIRTL_functions.R")
# to run on cpu
test<-run_single_point_analysis_sub_gpu("data/")
# to run on gpu
test_gpu<-run_single_point_analysis_sub_gpu("data/",backend="cupy")
# to run on 192 wells, first 12 columns of 384-well plate
test_half_plate<-run_single_point_analysis_sub_gpu("data/",wellset1=get_well_subset(1:16,1:12))

print(test) # we got alpha-beta pairs!
table(test$method) # note that same alpha-beta pairs can be called by different methods, so there are duplicates
```

### Input
The input data should be a folder with TCRalpha and TCRbeta repertoires from one plate in separate tab separated files from standard _mixcr_ output.

### Output
This will take TCRalpha and TCRbeta repertoires from the folder "data/" and pair them using the default parameters. The output will be saved in the working directory under _tmp_TIRTLoutput.tsv_ and also returned to R as _data.table_.

_wi_ - number of wells where alpha is found, but not beta

_wj_ - number of wells where beta is found, but not alpha

_wij_ - number of wells where both alpha and beta are found

_wa_ - number of wells where alpha is found

_wb_ -number of wells where beta is found

_alpha_nuc_seq_, _alpha_nuc_ - alpha CDR3nuc sequence

_beta_nuc_seq_, _beta_nuc_ - beta CDR3nuc sequence

_alpha_beta_ - alpha-beta pair concatenated nucleotide sequence

_method_ - method used for pairing t-shell, or MAD-HYPE

_r_ - for t-shell, the Pearson correlation coefficient between alpha and beta frequencies

_ts_ - for t-shell, the t-statistic of the correlation

_pval_ - for t-shell, the p-value of the correlation

_pval_adj_ - for t-shell, the adjusted p-value of the correlation

_loss_a_frac_ - the fraction of wells with lost alpha chain

_loss_b_frac_ - the fraction of wells with lost beta chain

_score_ - for MAD-HYPE, the score of the pairing (higher is better)

_cdr3a_ - alpha CDR3aa sequence

_va_ - alpha V gene

_ja_ - alpha J gene

_cdr3b_ - beta CDR3aa sequence

_vb_ - beta V gene

_jb_ - beta J gene



