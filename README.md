# RM_ASCA
This is a repository for reproducing the results in the manuscript "Repeated measures ASCA+ for analysis of longitudinal intervention studies with multivariate outcome data". The scripts are in MATLAB-code (2019b), and require the Statistics and Machine Learning Toolbox to work. For data availability, see the preprint for more information.
The preprint is available here: doi: https://doi.org/10.1101/2020.12.03.20243097

UPDATE
The preprint has now been published in PLOS Computational Biology (https://doi.org/10.1371/journal.pcbi.1009585). The code is available in the file paper_code.m, which uses the function RM_ASCA.m

The function RM_ASCA.m expects a table which, in addition to the multivariate response data, also includes the following variables: 
- A column named "ID" (for subject ID information)
- A column named "timepoint", which is a categorical variable indicating which timepoint each row is measured at
- A column named "group", which indicates the (treatment) group membership for each row. 

For data availability, please see the the publication for more information.

For questions, please contact torfinn.s.madssen@ntnu.no.
