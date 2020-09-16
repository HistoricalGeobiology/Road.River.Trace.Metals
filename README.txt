This repository contains the required R code and data files to produce the all of the main text figures in ‘An exceptional record of early to mid-Paleozoic redox change from the Road River Group, Yukon, Canada’ 
Erik A. Sperling, Michael J. Melchin, Tiffani Fraser, Richard G. Stockey, Una C. Farrell, Liam Bhajan, Tessa N. Browne, Devon B. Cole, Benjamin C. Gill, Alfred Lenz, David K. Loydell, Joseph Malinowski, Austin J. Miller, Stephanie Plaza-Torres, Beatrice Rodewald, Alan D. Rooney, Sabrina A. Tecklenburg, Jacqueline M. Vogel, Noah J. Planavsky, Justin V. Strauss.

Before running the R scripts included, download this folder and set it as your R working directory. Required data files should then load when called in each script, and pdf files will save within the same folder.

To replicate the analyses and plots presented here, the following R packages are required:
binom
deeptime
dplyr
egg
ggplot2
simpleboot

The deeptime package is currently only available on GitHub, and will need to be installed from there to reproduce Figures 3 and 4. This can be achieved by running the following commands in your R console (ignore the first line if you already have devtools installed).
install.packages("devtools")
devtools::install_github("willgearty/deeptime")

All other packages can be installed from CRAN. These scripts have been tested using R version 3.6.2 - Copyright (C) 2019 The R Foundation for Statistical Computing.
