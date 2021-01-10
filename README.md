SDERGM

julia code used to obtain the results in Score-Driven Exponential RanUtilities Graphs: A New Class of Time-Varying Parameter Models for Dynamical Networks https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3394593


This repository is still a work in progress and can contain bugs.

For a simple application of the SD beta model discussed in the paper start from simulationsAndPlotsPaper/beta_model/SD_Beta_model_example.jl

The data used for the empirical application regarding the eMid financial network is not publicaly available, hence not included here.

We have not prepared a list of dependencies. All the packages imported should be installed by:
	- using Pkg
	- Pkg.add("NAME_OF_PACKAGE")

Some methods and scripts related with global stastistics (not the fitness models) require Julia to call multiple R Libraries. For this interface to work, a working installation of R and statnet is required. We aknowledge the importance of this library for our analysis. For further information, please refer to  
Statnet Development Team
(Pavel N. Krivitsky, Mark S. Handcock, David R. Hunter, Carter T. Butts, Chad Klumb, Steven M. Goodreau, and Martina Morris) (2003-2020).
statnet: Software tools for the Statistical Modeling of Network Data. 
URL http://statnet.org