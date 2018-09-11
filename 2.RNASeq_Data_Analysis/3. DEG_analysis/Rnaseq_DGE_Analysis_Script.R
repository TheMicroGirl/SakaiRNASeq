cat("\014")
#prepare libraries
library(gplots)
library(limma)
library(edgeR)

#load the data
rnaseq <- read.delim("pairedEndsCountData.txt", header = T, row.names="Gene")

#filtering
rnaseqSubset <- rnaseq[rowSums(rnaseq)>10,]

#create group labels corresponding to the experiments
groupN<- c(rep("Control (MM)",3),rep("Experiment(LL)",3))

#construct DEGList, calculate Norm factors and normalise and prepare for linear modeling it using (voom)
deglist<- DGEList(counts = rnaseqSubset, group=groupN)
normalisedRnaseqData <- calcNormFactors(deglist)
normalisedRnaseqDataV<- voom(normalisedRnaseqData)

#linear modeling
fitRNA <- lmFit(normalisedRnaseqDataV)
fitRNA <- eBayes(fitRNA)
topTableRnaseq <- toptable(fitRNA,coef=2,n=Inf,sort="p")

plot(topTableRnaseq$logFC,-log(topTableRnaseq$P.Val), main = "Rnaseq volcano plot",xlab = "log2 fold change", ylab = "-log p Value")
#clustering plot
clusterData<- t(data.frame(normalisedRnaseqDataV))
distance <- dist(clusterData)
cluster <- hclust(distance)
plot(cluster, main = "RNA-Seq Cluster Dendrogram")

#Volcano plot
plot(topTableRnaseq $logFC,-log(topTableRnaseq $P.Val), main = "RNA-Seq volcano plot",xlab = "log2 fold change", ylab = "-log p Value")
#write table
write.table(topTableRnaseq, "RNAseqTopTable-Paired.txt" , sep="\t", row.names = T, col.names = T, quote = F )
   