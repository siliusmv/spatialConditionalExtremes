This repository contains the code used for creating all the results in the paper "An Efficient
Workflow for Modelling High-Dimensional Spatial Extremes", available as a preprint on
[https://arxiv.org/abs/2210.00760](https://arxiv.org/abs/2210.00760).

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
to load all `R` functions in the package and to compile and load the `Rcpp` functions
correctly. Finally, you need to call
```r
make_cgeneric("all")
```
in order to compile and link the necessary `cgeneric` models. Before calling this function, you might
need to change the C compiler `GCC` (by default equal to `gcc-12`) in the makefile
`cgeneric/Makefile` to a compiler that is available on your system.
In order to run the script for downloading all data from the online repositories (which is not
strictly necessary), you must first
install the Climate Data Operators (CDO) program. This is freely available at
[https://code.mpimet.mpg.de/projects/cdo](https://code.mpimet.mpg.de/projects/cdo).

# Running the code

All scripts for reproducing the results found in the paper are available in the `exec/`
folder. This folder contains two subfolders: `case_study/` and `simulation_studies/`.

## The `case-study/` folder

The `case-study/` folder contains all the scripts necessary for reproducing the results from the
case-study in the paper. The case-study scripts are:

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
  
## The `simulation_studies/` folder

The `simulation_studies/` folder contains the scripts necessary for running the two simulation
studies in the paper. These scripts are:

- `spde_toy_example.R`  
  This script contains the code for running the simulation study described in Section S1 of the
  supplementary material in the paper.
  It also contains the code for creating the table of coverage percentages displayed in that
  Section. The results of the simulation study are stored in the file
  `inst/extdata/results/low-rank.rds`.
- `simulation_study.R`  
  This script contains all the code for the simulation study of Section 5 of the paper.
  All results of the script are saved in the file `inst/extdata/results/simulation.rds`.
  In addition, several different plots are
  created and saved in the `inst/extdata/images/` folder. These are the files
  `design-of-experiment.pdf`, `simulation-posteriors.pdf` and `simulation-properties.pdf`.
