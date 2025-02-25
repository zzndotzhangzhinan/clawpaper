# Introduction

This repository contains the R code and data to reproduce numerical experiments in our paper: ``A Conformalized Empirical Method for Multiple Testing with Side Information''. Note that the code in this repository may not exactly replicate some of the results in the paper because of the randomness in simulating the data. 

There are several folds in this repository:

1. Fold ``main'' contains code for the numerical results in the main text, i.e., Figures 1-2 the simulation study section. It also contains code for Figure D.1, Figures E.1-E.2 and Figure E.6 in the Supplement.
2. Fold ``adad-storey'' contains code for Figure D.2 in the Supplement.
3. Fold ``adaptivity_level'' contains code for Figure D.3 in the Supplement.
4. Fold ``dependence'' contains code for Figures E.3-E.5 in the Supplement.
5. Fold ``srandom'' contains code for Figure E.7 in the Supplement.
6. Fold ``twosample'' contains code for Figure E.8 in the Supplement.
7. Fold ``yeast'' contains code for Figure E.9 in the Supplement.

**HOW TO IMPLEMENT** 
For example, run the following code to reproduce Figure 1 of our paper:
```R
source("figure-1.R")
```
