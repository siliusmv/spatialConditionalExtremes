# Setup

Start by calling 
```r
renv::load()
```
and follow the instructions that appear on the screen to download/install the correct libraries in
order to ensure reproducible results. You might have to install `renv` first, if you have not
already done so. Then you need to call
```r
devtools::load_all()
```
to compile the Rcpp functions correctly. Finally, you need to call
```r
make_cgeneric("all")
```
in order to compile and link the necessary `cgeneric` models. Before calling this function, you might
need to change the C compiler `GCC` (by default equal to `gcc-12`) in the makefile
`cgeneric/Makefile` to a compiler that is available on your system.

# Running the code

All scripts for reproducing the results found in the paper are available in the `exec/`
folder. This folder contains three subfolders: `simulation-studies/`, `case-study/` and
`tests/`. 

## The `case-study/` folder

The `case-study/` folder contains all the scripts necessary for reproducing the results from the
case-study in the paper. The scripts are enumerated to show the suggested order they should be
executed in. The case-study scripts are:

- `1-download-data.R`  
  This script downloads all the data needed for running the case study. It is not strictly
  necessary to run this script, as the processed data are already available in the file
  `inst/extdata/downloads/radar.rds`.
- `2-process-data.R`  
  This script processes the data from `1-download-data.R` to create the file
  `inst/extdata/downloads/radar.rds`. It also creates a map plot containing the data used for
  the case study, found in the file `inst/extdata/images/height-map.jpg`.
- `3-examine-marginal-distributions.R`  
  In this script, the marginal empirical cumulative distributions functions of the radar data are
  computed using sliding aggregation windows of different sizes. Different quantiles of the
  empirical marginals are then plotted for all locations, to examine properties of the different 
  aggregation sizes. The plots are saved in the file
  `inst/extdata/images/marginal-distributions.pdf`.
- `4-model-selection.R`  
  In this script, the model selection from Section 5.3 of the paper is performed. Results of the
  model selection are displayed in the files `inst/extdata/images/model-selection1.pdf` and
  `inst/extdata/images/model-selection2.pdf`. 
- `5-final-modelling.R`  
  In this script, the inference and model evaluation described in Sections 5.4 and 5.5 of the paper
  are performed. All model fits are saved in the file `inst/extdata/results/final-modelling.rds`,
  and all the log-scores are saved in the file
  `inst/extdata/results/final-modelling-log-scores.rds`. In addition, several different plots are
  created and saved in the `inst/extdata/images/` folder. These are the files
  `design-of-experiment.pdf`, `case-study-densities.pdf`, `case-study-densities2.pdf` and
  `log-score-rankings.pdf`.
  
## The `simulation-studies/` folder

The `simulation-studies/` folder contains all scripts for running the simulation studies described
in the paper. The scripts are enumerated to show the suggested order they should be executed
in. The simulation-study scripts are:

- `1-univariate-gaussian.R`  
  This script contains the code for running the simulation study described in Section 4.1 of the
  paper. It also contains the code for creating the table of coverage percentages displayed in that
  Section. The results of the simulation study are stored in the file
  `inst/extdata/results/univariate-gaussian.rds`.
- `2-low-rank.R`  
  This script contains the code for running the simulation study described in Section 4.2 of the
  paper. It also contains the code for creating the table of coverage percentages displayed in that
  Section. The results of the simulation study are stored in the file
  `inst/extdata/results/low-rank.rds`.
- `3-block-likelihood.R`  
  This script contains the code for running the simulation study described in Section 4.3 of the
  paper. It also contains the code for creating the table of coverage percentages displayed in that
  Section. The results of the simulation study are stored in the file
  `inst/extdata/results/block-likelihood.rds`.
- `4-conditional-extremes-parameter-recovery.R`  
  This script contains the code for evaluating the ability to recover parameters from simulated data
  with the implemented conditional extremes model in `R-INLA` (see Section 4.4 of the paper). The
  script depends on the results of the case study in `case-study/5-final-modelling.R`. Results of
  the simulation study are saved in the file `inst/extdata/results/parameter-recovery.rds`. 
- `5-conditional-extremes-adjustment.R`  
  This script contains the code for completing the simulation study described in Section 4.4 of the
  paper, and for creating the table of coverage percentages displayed in that Section. The
  script depends on the results of the case study in `case-study/5-final-modelling.R`. Results of
  the simulation study are saved in the two files `inst/extdata/results/conditional-theta-star.rds`
  and `inst/extdata/results/conditional-adjustment.rds`. 
- `6-gaussian-conditional-extremes.R`  
  This script contains the code for fitting the conditional extremes model to observations from a
  spatial Gaussian random field, and then evaluating frequency properties of unadjusted and
  adjusted posteriors. Results of the simulation study are saved in the file
  `inst/extdata/results/gaussian-conditional-extremes.rds`.
- `7-conditional-extremes-fixed-rho_b.R`  
  This script is almost a replicate of `5-conditional-extremes-adjustment.R`, but we fix `rho_b`
  instead of computing its posterior. The results of the simulation study are saved in the file
  `inst/extdata/results/conditional-adjustment-fixed-rho_b.rds`.
- `8-self-inconsistency.R`  
  This script demonstrates the problems that are caused by the lack of self-consistency of the
  conditional extremes model, as described in Appendix B of the paper.
- `9-constraining-Z.R`  
  This script estimates the correlation structure of a random field that has been constrained by
  subtraction, using Monte Carlo estimation. The estimated correlation structure is saved in the
  file `inst/extdata/images/constrained-correlation.pdf`.

## The `tests/` folder

The `tests/` folder contains two scripts for demonstrating that the implemented `Rcpp`-functions and
`inla.cgeneric`-functions work as they should. The test scripts are:

- `1-rcpp.R`  
  This script demonstrates that our implemented `Rcpp` functions for computing the log-likelihood of
  the spatial conditional extremes model are both correct and faster than plain `R` functions for
  computing the same thing.
- `2-cgeneric.R`  
  This script demonstrates the correctness of the implemented `cgeneric` functions used for computing
  the precision matrix of the SPDE approximation, by comparing the precision matrix created by the
  `cgeneric` functions with that created using the `INLA::inla.spde2.precision()` function in `R`.
