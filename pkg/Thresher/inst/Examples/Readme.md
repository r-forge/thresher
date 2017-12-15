---
author: Kevin R. Coombes and Min Wang
title: Code used for analyses in Thresher manuscript
date: 2017-12-08
---

# Overview
This folder contains the R code used to perform simulation and to
analyze breast cancer data sets in our paper titled _Thresher:
Determining the Number of Clusters while Removing Outliers_ and
published in ther jouirnal _BMC Bioinformatics_. 

## ExistingMethods

This subfolder contains scripts (and not packages) needed to implement
the existing methods we use for comparison. In particular, we needed
to make minor changes to the main `NbClust` driver function, since the
version in the package rather impolitely overwrites the seed for the
random number generator, thus breaking our simulation code.

## Simulations

This subfolder contains the code that performs the simulations
described in Section 3 of the manuscript.

## BRCA

This subfoilder contains the code to analyze GEO breast cancer data
sets, as described in Section 4 of the manuscript.
