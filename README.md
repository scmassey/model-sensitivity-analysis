# model-sensitivity-analysis
## Latin hypercube sampling and partial rank correlation coefficients for analyzing model parameter sensitivity.

LHS + PRCC is a useful method for investigating the sensitivity of a mathematical model to it's parameters. This can be useful in developing the model to understand how it behaves in various parameter regimes, as well as to understand better how uncertainty in your parameter estimates may impact the results given by the model.  

An overview of the procedure is provided as a pdf slide deck. The LHS method for parameter sampling in Monte Carlo studies was first developed by [McKay, Beckman, and Conover, 1979](https://doi.org/10.1080/00401706.1979.10489755) and was applied in conjunction with partial rank correlation coefficients for use in biomathematical models in [Blower and Dowlatabadi 1994](http://doi.org/10.2307/1403510). *A brief illustration of utility of this method as applied to the [proliferation-invasion-recruitment model](https://github.com/scmassey/2D_Proliferation-Invasion-Recruitment) will be on BioRxiv (as part of the mathematical oncology channel) in the near future.*

This repository contains code to conduct LHS+PRCC analysis in either matlab or python, depending on user preference. 

### Matlab
The Matlab file LHSPRCC.m is the main code file which calls the function DrawSamples.m to perform the Latin hypercube sampling step, any user-specified model functions for completing the Monte-Carlo Simulations, and either UnariedPRCC.m or VariedPRCC.m to compute partial rank correlation coefficients (at a single time/location index or at all times/locations). LHSPRCC.m also calls the functions plotSampleHists.m, plotSimulationOutput.m and plotUnvariedPRCC.m or plotVariedPRCC.m to display results from these various steps. 

Specifics about the sampled parameters are requested as user inputs in the command line, but a few code adjustments will need to be made as well to specify the particular model to be investigated as well as the output of interest for examining correlation between parameter space and model results.

Presently the code solves the linear function y=mx+b as a trivial example for the Monte Carlo simulations step. This is defined by the function testlinear.m which has the sampled parameters m and b. Note that this has a simple single output for computing PRCCs, but for models that are comprised of systems of equations with multiple dependent variables, the user will need to specify the particular output that they would like to investigate (either a single variable, or a sum or ratio of variables perhaps).


### Python
The Jupyter notebook LHS-PRCC.ipynb does the same procedure but is contained in a single file. Some user inputs can be done through interactive modules, while specifying the model and output of interest will need to be specified in the code itself. Further, the LHS-PRCC.ipynb notebook can be accessed using Google Colab so that users who are new to python may use the code and try it out without need to install a local python distribution. 
