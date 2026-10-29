This repository contains the data and scripts necessary to reproduce all analyses and figures for:

Schakowski, A., Deffner, D., Kortet, R., Niemelä, P.T., Kavelaars, M.M., Pykälä, M., & Kurvers, R.H.J.M (preprint) Social foragers rely more on dynamic cues than prior beliefs in unpredictable natural environments.

Paper here: link to preprint

To reproduce analyses, an additional folder "data" has to be downloaded from ZENODO (10.5281/zenodo.23769205), and it has to be placed into the same working directory as the remaining folders.

"wrapper_script_ecological_knowledge.r" sources all helper scripts and runs analyses. Scripts source raw data stored in "utils/data" and store processed data in the folder "utils/data/processed_data". 

To store figures, a folder named "6_figures" has to be created and placed in the main directory.

Software requirements

The analysis code was written in R (v4.0.3.) Statistical models are fit using the Stan MCMC engine (v2.32.6) and cmdstanr (v0.8.1) packages, which require a C++ compiler. Installation instructions are available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started and https://mc-stan.org/cmdstanr/articles/cmdstanr.html. See also the Stan user guide at https://mc-stan.org/users/documentation. All other required R packages can be found in "functions/library.R"

Detailed package versions can be found in "session_info.txt".
