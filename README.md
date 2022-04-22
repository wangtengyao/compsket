# compsket
Implementation of the complementary sketching algorithm for two sample testing of high-dimensional regression coefficients. 

# Description of files
## Python code
In the `./python/` folder
* `compsket.py`: file for the main algorithms
* `realdata.py`: file for implementing the real data example in Gao and Wang (2020)
* `example.ipynb`: IPython Notebook for the real data example
## R package
In `./R/` and `./man/` folders. Can be installed via `devtools::install_github('wangtengyao/compsket')` in `R`.
## matlab code 

## data
* `CD4_TREG_in_thymus.csv`: preprocessed data for the real data example in Section 5 of Gao and Wang (2020). 

# Reference
Gao, F. and Wang, T. (2020) Two-sample testing of high-dimensional linear regression coefficients via complementary sketching. Preprint, arxiv:2011.13624.
