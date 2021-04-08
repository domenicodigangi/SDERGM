# SDERGM Score Driven Exponential Random Graphs

This repository collects all code used to obtain the results in Score-Driven Exponential RanUtilities Graphs: A New Class of Time-Varying Parameter Models for Dynamical Networks https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3394593


Some methods and scripts require Julia to call multiple R Libraries. For this interface to work, a working installation of R and statnet is required. We aknowledge the importance of this library for our analysis. For further information, please refer to  
Statnet Development Team
(Pavel N. Krivitsky, Mark S. Handcock, David R. Hunter, Carter T. Butts, Chad Klumb, Steven M. Goodreau, and Martina Morris) (2003-2020).
statnet: Software tools for the Statistical Modeling of Network Data. 
URL http://statnet.org

## Caveats
- This repository is still a work in progress and can contain bugs.
- The paper has undergone one major revision that resulted in the rewriting of a large portion of the codebase. Not all the plots and analysis required a revision, and some applications were not updated to the last codebase version. Old non working scripts, that were used to generate the plots in the paper, are tagged by "OLD_REPO_VERSION" in their folders' names. They should run with the last commit of 2019. Atm there is no plan to update them.
- The data used for the empirical application regarding the eMid financial network is not publicaly available, hence not included here.