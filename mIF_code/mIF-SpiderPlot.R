#Proportions for each marker were collected using ImageJ
#This script was used to create the spider plot for the mIF results 

#Load necessary libraries 
library(Seurat)
library(dplyr)
#install.packages("fmsb")
library(fmsb)
#devtools::install_github("ricardo-bion/ggradar", 
#                          dependencies = TRUE)
library('ggradar')

outdir <- '/results/mIF/'
dir.create(outdir)

#### scale data to highlight differences ####
##Octogonal plots with core and edge and the four markers and 2 groups 
df_octagon_2groups.combined <- data.frame(
  Macrophage_Core = c(0.01175, 0.04502), #macrophages values *2 
  Bcell_Core = c(0.00891, 0.08568), #B cell values *18
  Tcell_Core = c(0.05375, 0.0535), #T cell values *50
  Fibroblast_Core = c(0.1557, 0.17995),
  Fibroblast_Edge = c(0.15294, 0.17872),
  Tcell_Edge = c(0.04775, 0.175), #T cell values *50
  Bcell_Edge = c(0.04275, 0.17424), #B cell values *18
  Macrophage_Edge = c(0.03581, 0.18742) #macrophage values *2 
)

rownames(df_octagon_2groups.combined) <- c("Non-cancerous","Cancerous")

# Add max and min rows (needed for fmsb)
df_octagon_2groups.combined <- rbind(rep(0.25,8), rep(0,8), df_octagon_2groups.combined)

#rotate octagon so it is flat on the top 
#df_octagon_2groups.combined.rotated <- df_octagon_2groups.combined[, c(2:ncol(df_octagon_2groups.combined), 1)]

df_rot <- df_octagon_2groups.combined
df_rot$dummy <- 0.0
df_rot <- df_rot[,c(ncol(df_rot), 1:(ncol(df_rot)-1))]

pdf(paste0(outdir, "SpiderPlot_Octagon_2groups-combined_rotated_turquoiseANDdarkgrey.pdf"))
radarchart(df_rot, axistype=1,
           pcol=c("darkgrey","darkturquoise"),
           pfcol = c(scales::alpha("darkgrey", 0.4),
                     scales::alpha("darkturquoise",0.15)),
           plwd=2, plty=1,
           cglcol="grey", cglty=1,
           axislabcol="grey", caxislabels=seq(0,0.25,0.05),
           cglwd=0.4, vlcex=0.8)
legend(x = 0.5, y = 1.3, legend = c("Non-cancerous","Cancerous"),
       col=c("darkgrey","darkturquoise"), lty=1, lwd=2, bty = "n")
dev.off()
