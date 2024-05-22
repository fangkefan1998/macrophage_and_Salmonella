rm(list = ls())

library(pheatmap)
library(ggplot2)

raw_data <- read.table("D:/Documents/R/Data/sp1-1.5h-heatmap/dual_salm.CSV", header = TRUE, sep = ",")

data <- raw_data[2:4]
row.names(data) <- raw_data$Geneid

pheatmap(data, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         cellwidth = 12, cellheight = 4)

col_note <- data.frame(Treatment = factor(c("Control", "Raw_1.5h", "Raw_24h")))
row.names(col_note) <- colnames(data)

pheatmap(data, scale = "row", border_color = "white", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         cellwidth = 15, cellheight = 4.5, annotation_col = col_note)


row_note <- data.frame(Pathway = factor(raw_data$pathway))
row.names(row_note) <- rownames(data)

p <- pheatmap(data, scale = "row", border=FALSE, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
              cellwidth = 16, cellheight = 14, annotation_col = col_note, annotation_row = row_note, gaps_col = 3, gaps_row = c(13,20,24,26,33,39,41,48,54,69,80,88),
              legend = TRUE, legend_breaks = c(-1,0,1), annotation_names_row = FALSE, annotation_names_col = FALSE, fontsize = 15)

p

