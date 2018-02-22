# epimutation_rate_analysis_R
This repository contains all the code for the analysis of the epimutation rates for the contexts CG, CHH and CHG. The biological problem tackled is described in the following paper: "Rate, spectrum, and evolutionary dynamics of spontaneous epimutations Adriaan van der Graaf a,1 , René Wardenaar a,1 , Drexel A. Neumann b , Aaron Taudt c , Ruth G. Shaw d , Ritsert C. Jansen a ,Robert J. Schmitz b,2 , Maria Colomé-Tatché c,2 , and Frank Johannes a,2". 

As an alternative approach this project estimated the rates using Markov chains with mean imputation using the EM - algorithm for missing values. For details see "report" document in repository concluding our work or see the following paper: "Estimation of the transition matrix of a discrete-time Markov chain Bruce A. Craig a, * and Peter P. Sendi b a b Department of Statistics, Purdue University, West Lafayette, USA Internal Medicine Outpatient Department, University of Basel, Switzerland", also in repository.

This work was done as part of the Topics in Computational Biology course at the Technical University of Munich (TUM), organized by Professor Fabian Theis. This particular project was supervised by Dr. Maria Colome Tatche, Junior Group Leader Computational Epigenomics Institute of Computational Biology, Helmholtz Zentrum Muenchen.

INSTRUCTIONS FOR REPRODUCABILITY:

The three R scripts in the general folder are the three main scripts for the contexts. Given that the data has been obtained, one needs to link these scripts to the helper functions (in helper functions folder) and the data and then simply execute. This will yield the result files. 

Once the result files are obtained the user can find the visualize results code in the helper function folder. With these script, yet again linked to the location of functions and data, the user can obtain the graphs saved in the "result" folder.

General comment on R-scripts "data_helper_functions" and "stat_functions". The stat_function script contains the EM-algorithm and bootstrap implementation. The "data_helper_function" contains general data preprocessing functions. These scripts are loaded into every main script in this project.

Data:

The data for this project needs to be abtained from the authoers of the paper mentioned in the first section of this README file. # epimutation_rate_analysis_R
