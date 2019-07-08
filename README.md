# model-sensitivity-analysis
##Latin hypercube sampling and partial rank correlation coefficients for analyzing model parameter sensitivity.

The procedure for sensitivity analysis is reviewed in the pdf slides... 

This repository contains code to conduct this analysis in either matlab or python, depending on user preference. 

### Matlab
The Matlab file LHSPRCC.m is the main code file which calls the functions DrawSamples.m, VariedPRCC.m and UnvariedPRCC.m, as well as any user specified model funtion files for running simulations with the sampled parameters.

To conduct the Latin hypercube sampling for your own particular model, start by editing LHSPRCC.m to take your model function when running simulations (see line XX). Currently the code solves the linear function y=mx+b contained in testlinear.m which has the parameters m and b as a trivial example. 
VariedPRCC.m and/or UnvariedPRCC.m may need to be adjusted to ensure that the model outputs of interest are those examined for correlation with parameter values. Further, these can be adjusted to accommodate ode or pde type models. 

When you run LHSPRCC.m, it will prompt (via the command line) for user inputs to specify the number of parameters to be sampled and the number of samples to draw. Because neither orthogonality nor space-filling are implemented, a larger number of samples is recommended to increase likelihood of better coverage of the parameter space.

Then LHSPRCC will call DrawSamples, which will prompt for user inputs to specify parameter names (for labeling in plots) and distribution details - i.e., the type of distribution and values to specify it. After the Latin hypercube sampling is complete, the histograms of the samples for each parameter are plotted to show the frequency of samples.  

Next, it will run Monte Carlo simulations using the randomly paired parameter samples using the model specified by the user; the output from these will be compared to the parameters using partial rank correlation coefficients. For multispecies models (WHAT is the phrase I actually want?), the user will need to decide whether a sum or ratio of outputs from the equations may be specified as the "output" to compare with PRCC. This will need to be specified within the PRCC function files themselves.

Further, there are two ways the code looks at PRCCs - one is for a single static point in time and/or space using a call to UnvariedPRCC, and another looks at variation over time or space using VariedPRCC. The user should specify indices for a static time and/or space as applicable in LHSPRCC.m prior to the call to VariedPRCC or UnvariedPRCC. 

Once the PRCC values have been completed, results will be plotted to compare the relative sensitivities of the parameters.


### Python
The Jupyter notebook LHS-PRCC.ipynb does similarly... blah edit after done with the notebook code. May also want to point to the google collab if use that.
