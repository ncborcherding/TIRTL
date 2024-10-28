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

For basic usage, you can run the following commands in R:
```R
source("TIRTL_functions.R")
test<-run_single_point_analysis_sub_gpu("data/")
test_gpu<-run_single_point_analysis_sub_gpu("data/",backend="cupy")
test_half_plate<-run_single_point_analysis_sub_gpu("data/",wellset1=get_well_subset(1:16,1:12)) # to run on 192 wells, first 12 columns
```

