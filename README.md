# FIA-Nfix-Estimate
Scripts used to estimate tree N fixation from FIA

Two methods were used: (1) Estimate from single time-point N-fixer abundances, and (2) estimate from N required to fulfil growth observed between two time points.

Data mergning and cleaning scripts:<br/>
`clean.R` <br/>
`clean_sapling.R` <br/>
`functions.R` <br/>

Method one (single time point abundances) scripts:<br/>
`do_withBootstrapping.R` <br/>

Method two (N required for growth) scripts:<br/>
`doGR.R` <br/>
`doGR_PerakisConstants.R` <br/>
`ndfa_GRnearby.R` <br/>
`ndfa_bootstrapping.R` <br/>

Further processing:<br/>
`N_dep_proc2.R` <br/>
`sensitivity_new.R` <br/>

Visualization scripts:<br/>
`figs_new.R` <br/>

Overview in `fia-est-fullworkflow.Rmd`, comparison to prior work in `prev-estimates.Rmd`
