library(ggplot2)
library(cowplot)

setwd(Sys.getenv("pca"))

plotList = list()

load("case_control_pca_plots_AFR.RData")
plotList[["AFR"]] = plots

load("case_control_pca_plots_AMR.RData")
plotList[["AMR"]] = plots

load("case_control_pca_plots_EUR.RData")
plotList[["EUR"]] = plots

load("case_control_pca_plots_FIN.RData")
plotList[["FIN"]] = plots


pdf("case_control_pca_plots_combined.pdf", width=10, height=10)
plot_grid(plotList[["AFR"]][[1]],plotList[["AFR"]][[2]],plotList[["AFR"]][[3]],
          plotList[["AMR"]][[1]],plotList[["AMR"]][[2]],plotList[["AMR"]][[3]],
          plotList[["EUR"]][[1]],plotList[["EUR"]][[2]],plotList[["EUR"]][[3]],
          plotList[["FIN"]][[1]],plotList[["FIN"]][[2]],plotList[["FIN"]][[3]],
          nrow=4, ncol=3, rel_heights=c(1,1,1,1))
dev.off()
