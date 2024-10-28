# TIRTL
Code for paired TCR data analysis from TIRTL-seq experiments.

## Installation
### Requirements

For GPU functionality (NVidia only) you will need to have CUDA installed. You can install it from [here](https://developer.nvidia.com/cuda-downloads).

R dependencies are _Matrix_ and _data.table_ packages, available by default with most distributions. If not, you can install them by running:
```R
install.packages("Matrix")
install.packages("data.table")
```

Python3 dependencies are _numpy_, _pandas_, _cupy_ (for GPU runs, optional). You can install them by running:
```bash
pip install numpy pandas
pip install cupy
```
## Usage


