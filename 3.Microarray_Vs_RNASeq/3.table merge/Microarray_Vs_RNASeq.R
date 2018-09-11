#read data in
microarray <- read.delim("microarrayTopTable.txt", header = T)
rna <- read.delim("RNAseqTopTable-Paired.txt", header = T)
blast <- read.delim("final_blast_Table.txt", header=T)

#make a table
microTable<- merge(blast, microarray, by.x ="query.name", by.y = "ProbeName")
microTableSmall<- microTable[,c("query.name", "subject","logFC")]
microRnaTable <- merge(microTableSmall, rna, by.x= "subject", by.y ="Gene" )

#plot Data
plot(microRnaTable$logFC.y~microRnaTable$logFC.x, main = "Microarray VS RNA-Seq logFC", xlab = "RNA-Seq logFC", ylab = "Microarray logFC")
plot(microRnaTable$logFC.y~microRnaTable$logFC.x, main = "Microarray VS RNA-Seq logFC", xlab = "RNA-Seq logFC", ylab = "Microarray logFC",pch=20, col=c(4,'skyblue','darkblue'))

#ks test 
web<-(ks.test(microRnaTable$logFC.y,microRnaTable$logFC.x))
# more plots
abline(lm(microRnaTable$logFC.x~microRnaTable$logFC.y),col = "red", lwd=2)
plot(lm(microRnaTable$logFC.x~microRnaTable$logFC.y))
summary(lm(microRnaTable$logFC.x~microRnaTable$logFC.y))

