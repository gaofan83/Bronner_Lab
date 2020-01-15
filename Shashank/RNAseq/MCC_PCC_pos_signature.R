#The script calculates log10fold_change using average MCC/PCC positive gene FPKM values and MCC/PCC negative gene FPKM values.

bulk_pos<-read.table("/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/cuffnorm_out/genes.fpkm_table", header=T, sep="\t", row.names=1)

mcc_pos<-bulk_pos[,c(1:192)]
pcc_pos<-bulk_pos[,c(193:384)]

mcc_neg<-read.table("/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/MCC_neg_bulk/cuffnorm_out/genes.fpkm_table", header=T, sep="\t", row.names=1)
pcc_neg<-read.table("/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/PCC_neg_bulk/cuffnorm_out/genes.fpkm_table", header=T, sep="\t", row.names=1)

mcc_ave<-data.frame(POS_AVE=rowMeans(mcc_pos), NEG_AVE=rowMeans(mcc_neg))
pcc_ave<-data.frame(POS_AVE=rowMeans(pcc_pos), NEG_AVE=rowMeans(pcc_neg))

mcc_ave$logFC<-log10(mcc_ave$POS_AVE/mcc_ave$NEG_AVE)
pcc_ave$logFC<-log10(pcc_ave$POS_AVE/pcc_ave$NEG_AVE)

write.table(mcc_ave, file="/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/MCC_FPKMave_log10FC.txt", quote=F, sep="\t")
write.table(pcc_ave, file="/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/PCC_FPKMave_log10FC.txt", quote=F, sep="\t")

gene_name<-read.table("/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/gene_id_name_uniq.txt", header=F, sep="\t", row.names=1)
colnames(gene_name)<-"Gene_Name"
mcc_ave_merge<-merge(mcc_ave, gene_name, by=0)
pcc_ave_merge<-merge(pcc_ave, gene_name, by=0)

write.table(mcc_ave_merge, file="/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/MCC_FPKMave_log10FC_name.txt", quote=F, sep="\t", row.names=F)
write.table(pcc_ave_merge, file="/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/PCC_FPKMave_log10FC_name.txt", quote=F, sep="\t", row.names=F)

mcc_ave_merge_sel<-mcc_ave_merge[mcc_ave_merge$POS_AVE>=1 & mcc_ave_merge$logFC>=1,]
pcc_ave_merge_sel<-pcc_ave_merge[pcc_ave_merge$POS_AVE>=1 & pcc_ave_merge$logFC>=1,]

write.table(mcc_ave_merge_sel, file="/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/MCC_FPKMave_log10FC_name_sel.txt", quote=F, sep="\t", row.names=F)
write.table(pcc_ave_merge_sel, file="/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/PCC_FPKMave_log10FC_name_sel.txt", quote=F, sep="\t", row.names=F)

mcc_ave_merge_sel2<-mcc_ave_merge[mcc_ave_merge$POS_AVE>=10 & mcc_ave_merge$logFC>=0.5,]
pcc_ave_merge_sel2<-pcc_ave_merge[pcc_ave_merge$POS_AVE>=10 & pcc_ave_merge$logFC>=0.5,]

write.table(mcc_ave_merge_sel2, file="/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/MCC_FPKMave_log10FC_name_sel2.txt", quote=F, sep="\t", row.names=F)
write.table(pcc_ave_merge_sel2, file="/home/fgao/Data_RNA/Bronner_Lab/SMARTseqData/PCC_FPKMave_log10FC_name_sel2.txt", quote=F, sep="\t", row.names=F)

