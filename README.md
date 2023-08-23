# Introduction

This repo contains scripts exploring how biased exposure-response relationships can be created by the choice of exposure metric. This work resulted in a forthcoming publication, and we hope this repo will help an interested reader further exposure these issues. The primary script and analysis look at comparing average concentration until an event to cycle 1 average concentration in a hypothetical oncology study, and how a biased analysis and conclusions could occur. 

# Packages

The packages used are captured in a renv lock file. R version 4.1.3 was used.

# Key Scripts

## tte-demo.Rmd

This is the primary script for the results and figures presented in the manuscript. 

## tte-demo-replications.R

This script repeats the core parts of the methodology in tte-demo.Rmd 1000 times. This supports the supplementary material and conclusions that the primary analysis is not due to random variability. 

## model/nonmem/simmod/pk2cmt.cpp

The [mrgsolve](https://mrgsolve.org/) model file used to generate exposures from a 2 compartment model. 

# Supplemental and Work in Progress Scripts

- dose-modification.Rmd: Work in progress exploration of how dose modifications throughout a trial can also induce bias in the exposure-response relationships. This occurs because lower doses only happen later in the trial, for patients who remain in the study and don't have death, progression, or serious adverse events, etc. 
- cmax-demo.Rmd: Brief analysis that any time-dependent exposure metric can have bias in the exposure response relationship, for example where Cmax increases cycle by cycle because of accumulation. 
- dose-modification-dose-response.Rmd: Work in progress exploration of time-varying covariates and created bias in the exposure-response relationship. 

# Further References

- [Causal Inference: What If by Miguel Hernan and James Robins] (https://www.hsph.harvard.edu/miguel-hernan/wp-content/uploads/sites/1268/2023/07/hernanrobins_WhatIf_19jul23.pdf) This a thourough book explaining causal inference and controlling for confounding, among other issues, in analyses. 
