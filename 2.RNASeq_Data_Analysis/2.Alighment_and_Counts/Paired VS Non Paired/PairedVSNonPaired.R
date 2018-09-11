#read data
rawPaired <- read.table("pairedEndsCountData.txt", header = T)
rawNOTPaired <- read.table("nonPairedCountData.txt",header = T)

#plot
plot((log(rawPaired$RNASeq18LL.3P))~((log(rawNOTPaired$RNASeq18LL.1.uP.2))), main = "Paired vs Unpaired Data - Sample 18LL-3",xlab= "Paired (log10 Counts)", ylab = "Unpaired (log10Counts)")                                    

#use as approperiate

















