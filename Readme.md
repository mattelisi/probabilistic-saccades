This repository contains the data and code used in the statistical analyses of the upcoming article "Gain control of saccadic eye movements is probabilistic", by Matteo Lisi, Joshua A. Solomon and Michael J. Morgan.

Analysis are documented through Rmarkdown notebooks (`.Rmd`), rendered to Html output. The repository contains one file demonstrating the analysis of the perceptual task of experiment 1 (see `exp1_perception.html`) and one file demonstrating the main analyses of saccade data of experiment 1 (`exp1_saccade.html`); the analysis of saccade data for experiment 2 and 3 used the same steps and methods, and is not reported here. There is also one script demonstrating the analysis of the range-effect in Experiment 1 and 3, `rangeEffect.R`; the hierarchical Bayesian model used in this analysis was written in Stan and is defined in the file `range_model.stan`.
Finally, there is one last file demonstrating in detail the modelling of the cost function for all the three saccade experiments (`analysis_cost_asymmetry.R`), and another one that illustrate the analysis of secondary saccades  (`secondary_saccades_analysis.R`).

Some of the analyses use custom functions from a library (`mlisi`) that I made to keep organized my frequently used functions, which is available here: [https://github.com/mattelisi/mlisi](https://github.com/mattelisi/mlisi)

The `data` folders contains the raw datasets (and any intermediate output of the analysis scripts). In particular, the raw datasets are:
- `exp1_perception.txt`, data of perceptual task, exp. 1
- `exp1_saccade.txt`, data from saccade task, exp. 1
- `exp2_luminance-fixed.txt`, data for exp. 2, fixed-luminance condition
- `exp2_size-fixed.txt`, data for exp. 2, fixed-size condition
- `exp3.txt`, data for exp. 3
- `saccades_xp123_allfit.txt`, saccade dataset, all 3 experiments pooled together, used in the analysis of cost asymmetry
- `secondary_saccades_xp123_allfit.txt`, secondary saccade dataset (for all 3 saccade experiments)

All dataset are in tab-separated text format, and the first line of each dataset indicate the variable names, which should be (hopefully) self-explanatory. The raw gaze recording files (Eyelink `edf` format) are not included in this repository, but are available upon request. The matlab code used to process raw gaze recordings and identify saccades is available at this GitHub repository [github.com/mattelisi/gaussianblobnoise-saccade](https://github.com/mattelisi/gaussianblobnoise-saccade). 

For any questions, please feel free to drop a line at: matteo.lisi@city.ac.uk
