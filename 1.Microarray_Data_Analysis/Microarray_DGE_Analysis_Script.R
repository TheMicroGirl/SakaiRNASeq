cat("\014")
library(limma)
library(gplots)

#read data in 
summaryTable<-read.delim("summaryTable.txt",check.names=FALSE,stringsAsFactors=FALSE)
microarrayData <-read.maimages(summaryTable[,"FileName"],source="agilent",green.only=TRUE)

#normalisation
correctedMicroarrayData <-backgroundCorrect(microarrayData,method="normexp")
correctedMicroarrayData <-normalizeBetweenArrays(correctedMicroarrayData,method="quantile")

#fileting of probes
neg95 <- apply(correctedMicroarrayData$E[correctedMicroarrayData$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95)) #95 percentile of negative probes
cutoff <- matrix(1.1*neg95,nrow(correctedMicroarrayData),ncol(correctedMicroarrayData),byrow=TRUE) # keep probes that are 10% brighter than the negative probes
filtered <- rowSums(correctedMicroarrayData$E > cutoff) >= 3

#removing of negative probes
noNegData <- correctedMicroarrayData[correctedMicroarrayData$genes$ControlType==0 & filtered,]

#build a vector
Treatment <- summaryTable[,"treatment"]
levels <- c("Minimal Media","Leaf Lysate")
Treatment <- factor(Treatment,levels=levels)
design <- model.matrix(~Treatment)
design

#linear modeling of probe level values
fit <- lmFit(noNegData,design)
fit <- eBayes(fit,trend=TRUE)
plotSA(fit, main="Probe-level")
summary(decideTests(fit[,-1]))


#linear modeling including the average of the  repeat probes
noRepData <- avereps(noNegData,ID=noNegData$genes[,"SystematicName"])
fitnR <- lmFit(noRepData,design)
fitnR <- eBayes(fitnR, trend=TRUE)
plotSA(fitnR, main="Gene-level")
summary(decideTests(fitnR[,-1]))
head(fitnR)

#create a dendrogram
plot(hclust(dist(t(noRepData$E))), xlab='Sample', main = "Microarray Cluster Dendrogram", ylab="DIstance", xlab=" ")

#calculate logFC
topTable <- toptable(fit, coef=2, genelist=fit$genes, n=Inf, sort='p')

#save results
write.table(topTable, "microarrayTopTable.txt", sep="\t", row.names = T, col.names = T, quote = F )

