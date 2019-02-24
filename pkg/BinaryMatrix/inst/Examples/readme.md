---
author: Kevin R. Coombes and Cailin E. Coombes
date: 24 February 2019
title: Preparing sample data for BinaryMatrix pacakge
---

## Inputs
The "raw" input data consists of two parts

1. **CML_Mitleman_2018.txt**  
  Contaisn the output of running the late-2018 early 2019 version (1.0)
  of the CVytoGPS code on all samples of chronic myelogenous leukemia,
  CML, from the public Mitelman database.
2. **ChrBands.txt**  
Contains the official names of all 916 cytogenetic bands (a.k.a.
"cytobands") that are partt of the latests ISCN standard.

## Processing
To convert the list of cytobands into an R object, run the script
`parseBands.R`. This creates a binary R data file in the `data` folder
called `lgfFeatures.rda`.

To convert the CML Mitelman data into an R  object, first run the
perl script `jsoinify.pl`. This script parses the (somewhat nasty)
output of the CytoGPS program into a well-behaved JSON file, called,
naturally enough. `CML_Mitelman_2018.json`. Next, run the R script
`Rify.R`. This parses the JSON file and converts it into a data
frame. It than creates a binary R data file, called `CMLData.rda`,
that includes thre object: `CMLData`, `CMLInfo`, and `lgfFeatures`.
