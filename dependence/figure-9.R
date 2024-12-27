source("exch.R")
grid_arrange_shared_legend(line1.1fdr,line1.2fdr,line1.3fdr,
                           line1.1tdp,line1.2tdp,line1.3tdp,nrow=2,ncol=3,position = "top")

source("exch-rho.R")
grid_arrange_shared_legend(line1.1fdr,line1.2fdr,line1.3fdr,
                           line1.1tdp,line1.2tdp,line1.3tdp,nrow=2,ncol=3,position = "top")
