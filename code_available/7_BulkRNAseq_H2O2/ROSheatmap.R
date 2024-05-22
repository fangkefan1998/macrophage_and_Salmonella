rm(list = ls())

library(pheatmap)
library(ggplot2)

raw_data <- read.table("D:/Documents/R/Data/ROSheatmap/ROSheatmap.CSV", header = TRUE, sep = ",")

data <- raw_data[2:7]
row.names(data) <- raw_data$Geneid

pheatmap(data, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         cellwidth = 20, cellheight = 2.2)


col_note <- data.frame(Treatment = factor(c("control", "control", "control", "ROS", "ROS", "ROS")))
row.names(col_note) <- colnames(data)


pheatmap(data, scale = "row", border_color = "white", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         cellwidth = 20, cellheight = 2.2, annotation_col = col_note)



row_note <- data.frame(Pathway = factor(raw_data$pathway))
row.names(row_note) <- rownames(data)


anno_pathway <- c("orange","green4","cadetblue2","deepskyblue3","darkorchid1","lightpink","yellow1","gray","antiquewhite","antiquewhite4", "aquamarine", "aquamarine4","coral")
names(anno_pathway) <- c("glycolysis", "TCA cycle", "oxidative phosphorylation", "amino acid metabolism", "pentose phosphate pathway", "nucleotide metabolism", "T3SS SP1", "T3SS SP2","Detoxifying enzymes", "Scavengers", "Repair systems", "Transcriptional regulators", "Thiol Chemistry")

ann_color <- list(Treatment = c(control = "slateblue3", ROS = "red2"), Pathway = anno_pathway)


p <- pheatmap(data, scale = "row", border="black", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
         cellwidth = 19, cellheight = 2.8, annotation_col = col_note, annotation_row = row_note, gaps_col = 3, gaps_row = c(15, 40, 58, 92, 95, 102, 130,163,176,178,186,197),
         legend = TRUE, legend_breaks = c(-1,0,1), annotation_names_row = FALSE, annotation_names_col = FALSE, fontsize = 15, annotation_colors = ann_color)

p


