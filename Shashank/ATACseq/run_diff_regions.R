###############R script
data<-read.table("atac_all_bed.txt",header=T, row.names=1)
data_s<-t(apply(data,1,sort))
data_sel_sort<-data_s[data_s[,9] > data_s[,8]*2, ]
data_sel<-data[rownames(data)%in%rownames(data_sel_sort),]

group<-colnames(data_sel)[apply(data_sel,1,which.max)]

data_out<-data.frame(data_sel,group)
data_out_sort<-data_out[order(group),]

write.table(data.frame(rownames(data_out_sort), data_out_sort$group), file="ATAC_all_sel.txt",quote=F, row.names=F, col.names=F)


data_group1<-data.frame(sel1=rowMeans(data[,c(1,3)]), data[,c(2,4,5,6,7,8,9)])
data_group1s<-t(apply(data_group1,1,sort))
data_group1_sel_sort<-data_group1s[data_group1s[,8] > data_group1s[,7]*2, ]
data_group1_sel<-data.frame(data_group1[rownames(data_group1)%in%rownames(data_group1_sel_sort),], max=data_group1_sel_sort[,8])
data_group1_sig<-data_group1_sel[data_group1_sel$sel1==data_group1_sel$max & !rownames(data_group1_sel)%in%rownames(data_out_sort),]
write.table(data.frame(rownames(data_group1_sig)), file="ATAC_all_group1.txt",quote=F, row.names=F, col.names=F)

data_group2<-data.frame(sel2=rowMeans(data[,c(1,7)]), data[,c(2,3,4,5,6,8,9)])
data_group2s<-t(apply(data_group2,1,sort))
data_group2_sel_sort<-data_group2s[data_group2s[,8] > data_group2s[,7]*2, ]
data_group2_sel<-data.frame(data_group2[rownames(data_group2)%in%rownames(data_group2_sel_sort),], max=data_group2_sel_sort[,8])
data_group2_sig<-data_group2_sel[data_group2_sel$sel2==data_group2_sel$max & !rownames(data_group2_sel)%in%rownames(data_out_sort),]
write.table(data.frame(rownames(data_group2_sig)), file="ATAC_all_group2.txt",quote=F, row.names=F, col.names=F)

data_group3<-data.frame(sel3=rowMeans(data[,c(3,5)]), data[,c(1,2,4,6,7,8,9)])
data_group3s<-t(apply(data_group3,1,sort))
data_group3_sel_sort<-data_group3s[data_group3s[,8] > data_group3s[,7]*2, ]
data_group3_sel<-data.frame(data_group3[rownames(data_group3)%in%rownames(data_group3_sel_sort),], max=data_group3_sel_sort[,8])
data_group3_sig<-data_group3_sel[data_group3_sel$sel3==data_group3_sel$max & !rownames(data_group3_sel)%in%rownames(data_out_sort),]
write.table(data.frame(rownames(data_group3_sig)), file="ATAC_all_group3.txt",quote=F, row.names=F, col.names=F)

