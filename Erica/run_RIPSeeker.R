library(RIPSeeker)

bamFiles <- list.files(".", ".bam$", recursive=TRUE, full.names=TRUE)

# run Ripseek
rip_calls <- ripSeek(bamPath=grep(pattern="Cardiac", bamFiles, value=TRUE), cNAME="Sha_Nulk_RNAseq_Trunk_rep1.sort.bam", 
	genomeBuild="galgal6", uniqueHit=TRUE, assignMultihits=TRUE, 
	rerunWithDisambiguatedMultihits=TRUE, binSize=NULL, strandType=NULL, 
	paired=FALSE, biomart="ensembl", biomaRt_dataset="ggallus_gene_ensembl", 
	goAnno="org.Gg.eg.db", exportFormat="gff3", annotateFormat="txt", 
	annotateType="TSS", outDir="./rip_output", padjMethod="BH", 
	logOddCutoff = 0, pvalCutoff=1, pvalAdjCutoff = 1, eFDRCutoff = 1)

# peak counts

ripRPKM <- computeRPKM(bamFiles[1], dataset="ggallus_gene_ensembl", paired=FALSE, moreGeneInfo=FALSE,
	justRPKM=FALSE, idType="ensembl_transcript_id", featureType="exon",txDbName="biomart",
	biomart="ensembl", by="tx", saveData="/home/fgao/Data_RNA/Bronner_Lab/ripRPKM")

rulebase.results <- rulebaseRIPSeek(bamFiles=bamFiles[1], cNAME=bamFiles[2], myMin=1,
			rpkmCutoff = 0.4, fcCutoff = 3, moreRIPGeneInfo=TRUE, idType="ensembl_transcript_id", by="tx",
			biomart="ensembl",dataset="ggallus_gene_ensembl", saveData="/home/fgao/Data_RNA/Bronner_Lab/rulebase.results")

df <- rulebase.results$rpkmDF

df <- df[order(df$foldchange, decreasing=TRUE), ]

write.table(df, file="/home/yourpath/csv", append=FALSE, quote=TRUE, sep = ";", eol = "\n", na = "NA", dec = ".", row.names=TRUE, col.names=TRUE, qmethod = c("escape", "double"), fileEncoding = "")



