#Prepping the bed files from 3p-Seq run so that we can feed that data into the t4 pipeline
#Need to take the bed format and move it into a 4 column format [chrom, terminationPosition, coverage, strand]
rm(list = ls())
library(dplyr)

sample <- read.delim("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/3seq/T4pipeline-master/BedIrem/316-TermSeqChlndct60.bed", header = FALSE)
sample <- sample[,-c(4,5)]

Minus <- sample[sample$V6 %in% "-",]
Plus <- sample[sample$V6 %in% "+",]

MinusCounts <- as.data.frame(table(Minus$V2))
PlusCounts <- as.data.frame(table(Plus$V3))

Minus <- distinct(Minus, V2, .keep_all = TRUE)
Minus <- arrange(Minus, V2)
Minus$coverage <- MinusCounts$Freq[match(Minus$V2, MinusCounts$Var1)]

Plus <- distinct(Plus, V3, .keep_all = TRUE)
Plus <- arrange(Plus, V3)
Plus$coverage <- PlusCounts$Freq[match(Plus$V3, PlusCounts$Var1)]

Minus <- Minus[,c(1,2,5,4)]
colnames(Minus) <- c("chrom","terminationPosition","coverage","strand")

Plus <- Plus[,c(1,3,5,4)]
colnames(Plus) <- c("chrom","terminationPosition","coverage","strand")


write.table(Minus, file = "316-TermSeqChlndct60_Minus.pooled.3p",quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)
write.table(Plus, file = "316-TermSeqChlndct60_Plus.pooled.3p",quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)



