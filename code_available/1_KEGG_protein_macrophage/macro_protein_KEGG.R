rm(list = ls())

library(clusterProfiler)
library(org.EcK12.eg.db)

macro <- read.table("D:/Documents/R/Data/MacroInCtl36/macrophage.CSV", header = TRUE, sep = ",")

LB <- read.table("D:/Documents/R/Data/MacroInCtl36/LB36h.CSV", header = TRUE, sep = ",")


MacroHigh <- macro$Accession[macro$FDR == "High"]

LBHigh <- LB$Accession[LB$FDR == "High"]



same <- intersect(MacroHigh, LBHigh)

uni2kegg <- function(x){
  out_kegg <- bitr_kegg(x, fromType = "uniprot", toType = "kegg", organism = "sey", drop = TRUE)
  return(out_kegg[2])
}

high_protein <- uni2kegg(same)
high_protein <- high_protein$kegg




library(stringr)

kk <- enrichKEGG(gene = high_protein, organism = "sey", keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
barplot(kk, showCategory = 100, title = "The KEGG Enrichment Analysis", xlab = "Counts", ylab = "Term")


a <- summary(kk)
write.csv(a, file = "D:/Documents/R/Data/MacroInCtl36/macroIN16hcontrol.csv")


library(ggplot2)

MG_result <- read.csv(file = "D:/Documents/R/Data/MacroInCtl36/macroIN16hcontrol.csv", header = TRUE)


MG_result <- MG_result[order(MG_result$p.adjust),]


MG_result <- MG_result[1:22,]


top10 <- data.frame(MG_result$Description, MG_result$Count, -log2(MG_result$p.adjust))
colnames(top10) <- c("Description","count","padj")



p <- ggplot(data=top10, aes(x=Description, y=count, fill=padj))


p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='black', size = 1.5),
                 axis.text.y=element_text(color="black",size=30))


p3 <- p2 + ylim(0,400) + scale_fill_gradient(low="#4497C4",high="#e20612")

p4 <- p3 + scale_x_discrete(limits=rev(top10[,1])) +labs(x="",y="Counts",title="Shared")+
  theme(text=element_text(size=45,  family="sans"))+
  theme(legend.text=element_text(size=20)) +theme(legend.key.size = unit(0.5, "inches"))+
  theme(axis.ticks.x=element_line(color="black",size=1.3,lineend = 10))+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size = rel(1.4)))+
  scale_y_continuous(expand=c(0.0035,0), limits = c(0,165))+
  labs(fill = expression(-log[2](padj)))


p4


png("D:/Documents/R/Data/MacroInCtl36/macro-SAME-16hcontrol.png",width=2000,height=1600)
print(p4)
dev.off()




LBdiff <- setdiff(LBHigh, MacroHigh)


high_protein <- uni2kegg(LBdiff)
high_protein <- high_protein$kegg



library(stringr)

kk <- enrichKEGG(gene = high_protein, organism = "sey", keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
barplot(kk, showCategory = 100, title = "The KEGG Enrichment Analysis", xlab = "Counts", ylab = "Term")


a <- summary(kk)
write.csv(a, file = "D:/Documents/R/Data/MacroInCtl36/LB36h-different.csv")


library(ggplot2)

MG_result <- read.csv(file = "D:/Documents/R/Data/MacroInCtl36/LB36h-different.csv", header = TRUE)


MG_result <- MG_result[order(MG_result$p.adjust),]


MG_result <- MG_result[1:4,]


top10 <- data.frame(MG_result$Description, MG_result$Count, -log2(MG_result$p.adjust))
colnames(top10) <- c("Description","count","padj")



p <- ggplot(data=top10, aes(x=Description, y=count, fill=padj))


p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='black', size = 1.5),
                 axis.text.y=element_text(color="black",size=30))


p3 <- p2 + ylim(0,400) + scale_fill_gradient(low="#00abf0",high="#e20612")

p4 <- p3 + scale_x_discrete(limits=rev(top10[,1])) +labs(x="",y="Counts",title="Control LB36h Unique")+
  theme(text=element_text(size=42,  family="sans"))+
  theme(legend.text=element_text(size=20)) +theme(legend.key.size = unit(0.5, "inches"))+
  theme(axis.ticks.x=element_line(color="black",size=1.3,lineend = 10))+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size = rel(1.4)))+
  scale_y_continuous(expand=c(0.0035,0), limits = c(0,165))+
  labs(fill = expression(-Log[2](padj)))

p4



png("D:/Documents/R/Data/MacroInCtl36/lB36h-different.png",width=1700,height=400)
print(p4)
dev.off()


MacroHigh[length(LBHigh)] <- NA
veenData <- data.frame(LBHigh, MacroHigh)


write.csv(veenData, file = "D:/Documents/R/Data/MacroInCtl36/vennData.csv")






