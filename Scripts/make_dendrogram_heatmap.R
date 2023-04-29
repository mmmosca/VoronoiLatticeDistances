#
# Permission is granted to copy, distribute and/or modify the documents
# in this directory and its subdirectories unless otherwise stated under
# the terms of the GNU Free Documentation License, Version 1.1 or any later version 
# published by the Free Software Foundation; with no Invariant Sections, 
# no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
# is available at the website of the GNU Project.
# The programs and code snippets in this directory and its subdirectories
# are free software; you can redistribute them and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation; either version 2 of the License, or (at your option) any later
# version. This code is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
#
library(ggplot2)
library(reshape2)
library(ggdendro)
library(grid)

args <- commandArgs(trailingOnly=TRUE)

f2_t2l <- args[1]
outpath <- "."
dataframe_t2 <- read.csv(f2_t2l, header=TRUE, sep=",",dec=".")
n_crystals <- nrow(dataframe_t2)

rownames(dataframe_t2) <- dataframe_t2$ID
crystal_ids <- dataframe_t2$ID

vormetrics_distance_matrix <- as.matrix(dataframe_t2[, seq(2,ncol(dataframe_t2),1)])
maximum <- max(vormetrics_distance_matrix)

# Clustering 
hclusters <- hclust(d = as.dist(vormetrics_distance_matrix), method = "complete")

# Dendrogram
vormetrics_dendro <- as.dendrogram(hclusters)
dendrogram.plot <- ggdendrogram(data = vormetrics_dendro, rotate = F) +
                    theme(axis.text.y = element_text(size=16, hjust = 0.5, vjust = -15),
                        axis.text.x = element_text(size=14, vjust = 0.45, hjust = 1.0)) +
                    scale_y_continuous(breaks=seq(0, maximum, 1))
vormetrics_dendro_order <- order.dendrogram(vormetrics_dendro)
ggsave( filename=paste( outpath, "\\dendrogram.png", sep=""), plot = dendrogram.plot, width = 18, height = 11 )

# Heatmap with dendrogram order
#Scale to [0,255]

scale_fun <- function(v) ((v / maximum) * 255)
vormetrics_distance_matrix_colorcode <- vormetrics_distance_matrix
vormetrics_distance_matrix_colorcode[] <- vapply(vormetrics_distance_matrix_colorcode, scale_fun, numeric(1))

vormetrics_melted <- melt(vormetrics_distance_matrix_colorcode, factorsAsStrings = TRUE)
vormetrics_melted$Var1 <- factor(x = vormetrics_melted$Var1,
                            levels = crystal_ids[vormetrics_dendro_order],
                            ordered = T)
vormetrics_melted$Var2 <- factor(x = vormetrics_melted$Var2,
                            levels = crystal_ids[vormetrics_dendro_order],
                            ordered = T)

vormetrics.heatmap.plot <- ggplot(vormetrics_melted, aes(x=as.factor(Var1), y=as.factor(Var2))) +
    geom_tile(aes(fill = value)) +
    theme(axis.text.y = element_text(size=30), 
            axis.text.x = element_text(angle=90, vjust=0.2, hjust = 0.95, size=30),
            legend.position = "top",
            legend.title = element_text(size=20),
            legend.text = element_blank(),
            legend.key.size = unit(2,"line"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
    scale_fill_gradient( low="white", high="black" ) +
    labs(fill = "Distance Gradient", x = "Crystal IDs", y = "Crystal IDs")
ggsave( filename=paste( outpath, "\\heatmap.png", sep=""), plot = vormetrics.heatmap.plot, width = 20, height = 18 )
