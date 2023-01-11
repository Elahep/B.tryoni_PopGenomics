We calculated directional relative migration rates and quantified asymmetry in patterns of gene-flow using the DivMigrate function of the R package diveRsity.

"""

library(diveRsity)
library(corrplot)

divMigrate(infile="./Qff_all.gen", outfile="Gst_Qff", stat="gst", para=TRUE, plot_network=FALSE)
data <- as.matrix(read.table("Gst_Qff"))
#Heatmap
corrplot(data, method="color", order = "alphabet", is.corr = FALSE)

"""




