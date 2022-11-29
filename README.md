This repository contains the code used for creating all the results in the paper CONTINUE!!!

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
folder. This folder contains two subfolders: `case_study/` and
`tests/`, and the script `simulation_study.R`

## The `case-study/` folder

The `case-study/` folder contains all the scripts necessary for reproducing the results from the
case-study in the paper. The scripts are enumerated to show the suggested order they should be
executed in. The case-study scripts are:

- `download_data.R`  
  This script downloads all the data needed for running the case study. It is not strictly
  necessary to run this script, as the processed data are already available in the file
  `inst/extdata/downloads/radar.rds`.
- `process_data.R`  
  This script processes the data from `download_data.R` to create the file
  `inst/extdata/downloads/radar.rds`. It also creates a map plot containing the data used for
  the case study, found in the file `inst/extdata/images/height-map.jpg`.
- `case_study.R`  
  In this script, the inference and model evaluation described in Section 6 of the paper
  are performed. All results of the script are saved in the file `inst/extdata/results/case-study.rds`.
  In addition, several different plots are
  created and saved in the `inst/extdata/images/` folder. These are the files
  `case-study_bad_properties.pdf`, `case-study_posteriors.pdf` and `case-study_properties.pdf`.
  
## The `simulation_study.R` script

This script contains all the code for the simulation study of Section 5 of the paper.
All results of the script are saved in the file `inst/extdata/results/simulation.rds`.
In addition, several different plots are
created and saved in the `inst/extdata/images/` folder. These are the files
`design-of-experiment.pdf`, `simulation-posteriors.pdf` and `simulation-properties.pdf`.

## The `tests/` folder

The `tests/` folder contains two scripts for demonstrating that the implemented `Rcpp`-functions and
`inla.cgeneric`-functions work as they should. The test scripts are:

- `rcpp.R`  
  This script demonstrates that our implemented `Rcpp` functions for computing the log-likelihood of
  the spatial conditional extremes model are both correct and faster than plain `R` functions for
  computing the same thing.
- `cgeneric.R`  
  This script demonstrates the correctness of the implemented `cgeneric` functions used for computing
  the precision matrix of the SPDE approximation, by comparing the precision matrix created by the
  `cgeneric` functions with that created using the `INLA::inla.spde2.precision()` function in `R`.
