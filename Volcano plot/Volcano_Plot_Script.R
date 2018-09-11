library (ggplot2)

#read data
microarray <- read.table("microarrayTopTable.txt", header = T)
rna <- read.table("RNAseqTopTable-Paired.txt", header = T)

#put the treshold value
rna$threshold <- as.factor(abs(rna$logFC) > 1 & rna$P.Value < 0.05/5496)

#plot data 
a <- ggplot(data=rna, aes(x=logFC, y=-log10(P.Value), colour=threshold))
b <- a + geom_point(alpha=0.4, size=1.75)
c <- b + theme(legend.position = "none")
d <- c + xlim(c(-10, 10)) + ylim(c(0, 15)) 
e <- d + xlab("log2 fold change") + ylab("-log10 p-value")

#repeat for microarry



