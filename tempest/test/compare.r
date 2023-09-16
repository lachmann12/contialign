setwd("../../contialign/varconti")

kallisto = read.table("output/kallisto/abundance.tsv", stringsAsFactors=F, sep="\t", header=T)
kcounts = as.numeric(kallisto[,4])
names(kcounts) = kallisto[,1]
rust = read.table("output/counts.tsv", stringsAsFactors=F, sep="\t")
rcounts = as.numeric(rust[,2])
names(rcounts) = rust[,1]

cc = cor(log2(1+kcounts)[names(rcounts)], log2(1+rcounts))
plot(log2(1+kcounts), log2(1+rcounts)[names(rcounts)], pch=".", cex=1, ylab="rust counts", xlab="kallisto counts", main=paste(round(cc, digits=2), sum(rcounts), round(sum(kcounts))))
abline(0,1, col="blue")

