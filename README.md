# Graph-structured variable selection with Gaussian Markov random field horseshoe prior.

M. Denis and M.G. Tadesse

07/28/2021

![](logo.png)

This repository contains the R code for fitting the HS-GMRF model along with a detailed simulation example to simulate data as presented in the Annals of Applied Statistics paper "Graph-structured variable selection with Gaussian Markov random field horseshoe prior" by M. Denis and M.G. Tadesse (submitted).

`MainFunction.R` file contains the R code for the HS-GMRF model.

`UtilFunctions.R` file contains the R code for constructing the contrast and adjacency matrices and for computing the Matthews correlation coefficient and mean squared error.

`call.R` ffile is a R script to run HS-GMRF on dataset simulated as in the paper.
