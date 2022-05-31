# Intro

This repo contains the Code and SampleData for our paper ***Mining hourly population dynamics by activity type based on decomposition of sequential snapshot data***.

Our paper propose a framework of population mixture analysis that uses sequential snapshot data to generate dynamic population distributions by activity type. The framework contains four steps：constructing the population mixture model, decomposing the population mixture, estimating the dynamic population distribution, and validating population estimation
  
![Framework](/img/framework.png)

# Table of contents

- [Intro](#intro)
- [Structure](#structure)
- [Usage](#usage)

# Structure
`Code` directory contains all the code to run the experiments:

+ `DecomposeMixture.m` is the code to run the framework with experimental parameters. It contains how to prepocess total population’s curve, decompose the population mixture and estimate dynamic population distributions by activity type.

+ `SSPP.m` and `GNMF_KL.m` are the code to define SSPP-NMF algorithme used for decomposing the population mixturerun.

+ `hyperNnls.m` and `hyperNormalize` are the code to support the estimate temporal population distributions for each activity.

`SampleData` directory contains all the sample data described in our paper.

+ `RTUD_Sample.tif` is the sequential snapshot data for sample region.

+ `RTUD_ISODATA.tif` is a cluster image for sample region.

# Usage
Step 1. Set your working director to `Code`

Step 2. Run `DecomposeMixture.m`. 
When finished running the code, `Result` directory cantains the population distributions for each type.




