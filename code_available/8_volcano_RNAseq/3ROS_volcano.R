rm(list = ls())

library(DESeq2)

#读入原始数据
E16_counts <- read.table("D:/Documents/R/Data/3ROSreplicates20220506/E1-6_counts.CSV", header = TRUE, sep = ",")


#构建处理组对照组表格
countData <- E16_counts[-1]
row.names(countData) <- E16_counts$Geneid

condition <- factor(c("control", "control", "control", "ROS", "ROS", "ROS"))

colData <- data.frame(row.names=colnames(countData), condition)


#准备正式构建dds矩阵
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design= ~ condition)

dds2 <- DESeq(dds)   #对原始dds进行normalize
#查看结果与信息
resultsNames(dds2)
res <- results(dds2)
summary(res)


write.csv(res, file = "D:/Documents/R/Data/3ROSreplicates20220506/deseq_res_3ROS.csv", row.names = TRUE)
#######################################################################################################################

#读取差异基因并且提取出sp1和sp2的基因，分别保存用于作图
volData <- read.table("D:/Documents/R/Data/3ROSreplicates20220506/deseq_res_3ROS.csv", sep = ",", header = TRUE)

sp1 <- c("invG","prgH","prgK","spaP","spaQ","spaR","spaS","invA","orgA","spaO","orgB","invC","invI","prgI","prgJ","invJ","sipD","sipB","invH","invE","avrA","sipA","sipB","sipC","sipD","sopA","sopB","sopE","sopE2","sptP","slrP","sopD","sspH1","steA","steB")
sp2 <- c("ssaC","ssaD","ssaJ","ssaR","ssaS","ssaT","ssaU","ssaV","ssaQ","ssaK","ssaN","ssaO","ssaG","ssaI","ssaP","sseC","sseD","sseF","sseG","pipB2","sifA","sopD2","gtgE","sseJ","sseL","steC","spvB","spvC","gtgA","sseK1","sseK2","sseK3","sspH1","sspH2","slrP","gogB","spvD","sseI","steD","pipB","sifB","srfJ","steB","steE")

sp1_deseq <- volData[(volData$X %in% sp1), ]
sp2_deseq <- volData[(volData$X %in% sp2), ]


write.csv(sp1_deseq, file = "D:/Documents/R/Data/3ROSreplicates20220506/sp1_deseq.csv")
write.csv(sp2_deseq, file = "D:/Documents/R/Data/3ROSreplicates20220506/sp2_deseq.csv")





















