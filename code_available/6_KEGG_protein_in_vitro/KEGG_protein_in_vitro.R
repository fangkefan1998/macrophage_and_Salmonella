rm(list = ls())

library(clusterProfiler)
library(org.EcK12.eg.db)


three <- read.table("D:/Documents/R/Data/Macro-LB-ROS/threeOverlap.CSV", header = TRUE, sep = ",")
mac <- read.table("D:/Documents/R/Data/Macro-LB-ROS/macrophage.CSV", header = TRUE, sep = ",")

three <- three[(three$FDR == "High"), ]


hp <- three$Accession[!is.na(three$HP)]

md <- three$Accession[!is.na(three$MD)]

lb36h <- three$Accession[!is.na(three$LB36h)]

infection <- mac$Accession


overlap3 <- intersect(intersect(infection, lb36h), hp)




uni2kegg <- function(x){
  out_kegg <- bitr_kegg(x, fromType = "uniprot", toType = "kegg", organism = "sey", drop = TRUE)
  return(out_kegg[2])
}

high_protein <- uni2kegg(overlap3)
high_protein <- high_protein$kegg





library(stringr)

R.utils::setOption("clusterProfiler.download.method",'auto')
kk <- enrichKEGG(gene = high_protein, organism = "sey", keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
barplot(kk, showCategory = 100, title = "The KEGG Enrichment Analysis", xlab = "Counts", ylab = "Term")


a <- summary(kk)
write.csv(a, file = "D:/Documents/R/Data/Macro-LB-ROS/treeOverlap-kegg.csv")


library(ggplot2)

MG_result <- read.csv(file = "D:/Documents/R/Data/Macro-LB-ROS/treeOverlap-kegg.csv", header = TRUE)


MG_result <- MG_result[order(MG_result$p.adjust),]


MG_result <- MG_result[1:24,]


top10 <- data.frame(MG_result$Description, MG_result$Count, -log2(MG_result$p.adjust))
colnames(top10) <- c("Description","count","padj")



p <- ggplot(data=top10, aes(x=Description, y=count, fill=padj))


p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='black', size = 1.5),
                 axis.text.y=element_text(color="black",size=60), 
                 axis.text.x = element_text(color = "black", size = 45),
                 axis.title.x = element_text(size = 50),
                 legend.title = element_text(size = 40))


p3 <- p2 + ylim(0,400) + scale_fill_gradient(low="dodgerblue3",high="red3")

p4 <- p3 + scale_x_discrete(limits=rev(top10[,1])) +labs(x="",y="Counts",title=" ")+
  theme(text=element_text(size=35,  family="sans"))+
  theme(legend.text=element_text(size=40)) +theme(legend.key.size = unit(0.5, "inches"))+
  theme(axis.ticks.x=element_line(color="black",size=1.3,lineend = 10))+
  theme(axis.ticks.y = element_blank())+
  scale_y_continuous(expand=c(0.0035,0), limits = c(0,180))+
  labs(fill = expression(-Log[2](padj)))#



p4



png("D:/Documents/R/Data/Macro-LB-ROS/LBmacHPOverlap-kegg.png",width=2000,height=1700)
print(p4)
dev.off()




