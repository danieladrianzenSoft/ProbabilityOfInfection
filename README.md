<h2>Description</h2>

ProbabilityOfInfection is a collection of MATLAB functions that aims to estimate the probability of HIV-1 infection in women after sexual exposure. The model solves a system of simultaneous partial and ordinary differential equations that characterize HIV viral transport, anti-HIV drug transport, cell migration, and infection kinetics. By defining parameters as distributions and using a Monte-Carlo simulation approach, a user-defined number of simulations are run, and the number of positive infections to ultimately estimate the probability of HIV infection (POI).

"complete_analysis.m" is the main file that sets up the Monte Carlo simulation parameters, and executes them in parallel using parfor. Model parameters are defined in "initializeParameters.m" and in "jmp_doe.xlsx" depending on the nature of the parameter. 

"POI_BinaryCR_integrated.m" defines each simulation itself, and calls the corresponding functions depending on the drug administration conditions.
