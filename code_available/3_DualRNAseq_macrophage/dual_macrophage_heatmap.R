rm(list = ls())

library(pheatmap)
library(ggplot2)

raw_data <- read.table("D:/Documents/R/Data/D5-37-macrophage-gene/dualRNAseq.CSV", header = TRUE, sep = ",")

data <- raw_data[2:4]
row.names(data) <- raw_data$gene

pheatmap(data, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE,
         cellwidth = 12, cellheight = 5)

col_note <- data.frame(Treatment = factor(c("raw-control", "raw-infection-1.5h", "raw-infection-24h")))
row.names(col_note) <- colnames(data)

pheatmap(data, scale = "row", border_color = "white", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         cellwidth = 15, cellheight = 4.5, annotation_col = col_note)



row_note <- data.frame(Pathway = factor(raw_data$pathway))
row.names(row_note) <- rownames(data)

p <- pheatmap(data, scale = "row", border=FALSE, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
              cellwidth = 16, cellheight = 5, annotation_col = col_note, annotation_row = row_note, gaps_col = 0, gaps_row = c(13,26,29,46,74,86),
              legend = TRUE, legend_breaks = c(-1,0,1), annotation_names_row = FALSE, annotation_names_col = FALSE, fontsize = 15,
              color = colorRampPalette(c("darkgreen", "white", "darkorange1"))(100))

p





