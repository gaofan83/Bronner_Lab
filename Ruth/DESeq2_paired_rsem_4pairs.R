library(tximportData)
library(tximport)
library(readr)
library(DESeq2)

directory<-'./rsem/'
files<-grep("*genes.results",list.files(directory),value=TRUE)

#2nd pair switched
KO<-files[c(1,2,4,3,7:10)]

df1<-data.frame(sampleID =factor(rep(1:4, each=2)),
                condition=factor(rep(c("Left", "right"),4)))

sampleTable<-data.frame(sampleName=KO, fileName=KO, df1)

#import files
txi <- tximport(paste(directory, KO, sep=""), type = "rsem", txIn = FALSE, txOut = FALSE)
names(txi)

# cond is a vector with conditions to be used for differential analysis
cond <- rep(c("Left", "right"),4)
sampleTable <- data.frame(sampleID =factor(rep(1:4, each=2)), condition = factor(cond))
rownames(sampleTable) <- colnames(txi$counts)

txi$length[txi$length <= 0] <- 1

ddsRSEM <- DESeqDataSetFromTximport(txi, sampleTable, ~ sampleID + condition)

#colData(dds)<-DataFrame(df1)
dds<-DESeq(ddsRSEM)

res<-results(dds)
write.table(res,file="RSEM_DEG_KO_DESEQ2_pair4.txt",quote=FALSE,sep="\t")
resOrdered<-res[order(res$padj),]
res_sig = subset(resOrdered, padj<0.1)
write.table(res_sig,file="RSEM_DEG_KO_DESEQ2_SIG_pair4.txt",quote=FALSE,sep="\t")

fpkm<-read.table("/home/fgao/Data_RNA/Bronner_Lab_RNA/Ruth_02032023/rsem/FPKM_all.txt", header=T, sep="\t", row.names=1)
fpkm_p1<-fpkm[,c(1,2)]
fpkm_p2<-fpkm[,c(4,3)]
fpkm_p3<-fpkm[,c(7,8)]
fpkm_p4<-fpkm[,c(9,10)]
fpkm_sel<-fpkm[,c(1:4,7:10)]

zscore_p1<-t(apply(t(fpkm_p1), 2, scale))
zscore_p2<-t(apply(t(fpkm_p2), 2, scale))
zscore_p3<-t(apply(t(fpkm_p3), 2, scale))
zscore_p4<-t(apply(t(fpkm_p4), 2, scale))

zscore<-t(apply(t(fpkm_sel), 2, scale))

colnames(zscore)<-colnames(fpkm_sel)
write.table(zscore, file="/home/fgao/Data_RNA/Bronner_Lab_RNA/Ruth_02032023/rsem/zscore_4pairs.txt", quote=F, sep="\t")

