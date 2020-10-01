# FIA-Nfix-Estimate
Scripts used to estimate tree N fixation from FIA

Two methods were used: (1) Estimate from single time-point N-fixer abundances, and (2) estimate from N required to fulfil growth observed between two time points.

Data mergning and cleaning scripts:
<clean.R>
clean_sapling.R
functions.R

Method one (single time point abundances) scripts:
do_withBootstrapping.R

Method two (N required for growth) scripts:
doGR.R
doGR_PerakisConstants.R
ndfa_GRnearby.R
ndfa_bootstrapping.R

Further processing:
N_dep_proc2.R
sensitivity_new.R

Visualization scripts:
figs_new.R

Overview in fia-est-fullworkflow.Rmd, comparison to prior work in prev-estimates.Rmd
