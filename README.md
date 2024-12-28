# Introduction

This repository contains the R code and data to reproduce numerical experiments in our paper: ``A Conformalized Empirical Method for Multiple Testing with Side Information''. Note that the code in this repository may not exactly replicate some of the results in the paper because of the randomness in simulating the data. 

There are several folds in this repository:

1. Fold ``main'' contains code for the numerical results in the main text, i.e., Figures 1-2 and the application section. It also contains code for Figures 3-4, Figures 7-8 and Figure 11 in the Supplement.
2. Fold ``adad-storey'' contains code for Figure 5 in the Supplement.
3. Fold ``adaptivity_level'' contains code for Figure 6 in the Supplement.
4. Fold ``dependence'' contains code for Figures 9-10 in the Supplement.
5. Fold ``srandom'' contains code for Figure 12 in the Supplement.
6. Fold ``twosample'' contains code for Figure 13 in the Supplement.
7. Fold ``yeast'' contains code for Figure 14 in the Supplement.

**HOW TO IMPLEMENT** 
For example, run the following code to reproduce Figure 1 of our paper:
```R
source("figure-1.R")
```
