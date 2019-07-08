# model-sensitivity-analysis
Latin hypercube sampling and partial rank correlation coefficients for analyzing model parameter sensitivity.

The procedure for sensitivity analysis is reviewed in the pdf slides... This repository contains code to conduct this analysis in either matlab or python, depending on preference. 

The Matlab file LHSPRCC.m is the main code file which calls the functions DrawSamples.m, VariedPRCC.m and UnvariedPRCC.m, as well as any user specified model funtion files for running simulations with the sampled parameters. To conduct the Latin hypercube sampling, start by editing LHSPRCC.m to take your model function when running simulations (see line XX). Currently the code solves the linear function y=mx+b contained in testlinear.m which has the parameters m and b as a trivial example. 
VariedPRCC.m and/or UnvariedPRCC.m may need to be adjusted to ensure that the model outputs of interest are those examined for correlation with parameter values.

The Jupyter notebook LHS-PRCC.ipynb does similarly... blah edit after done with the notebook code.
