mt_genes <- read.csv(file="Human.MitoCarta2.0.csv",header=TRUE)
all_genes <- read.csv(file="all_human_genes.csv",header=TRUE)
check <- na.exclude(mt_genes$Symbol[mt_genes$TargetP_Score == 1])
as.vector(check)
str(check)
is.na(mt_genes$TargetP_Score)
mt_genes[is.na(mt_genes$TargetP_Score), c("Symbol","TargetP_Score","MCARTA2.0_score")]
summary(mt_genes$TargetP_Score)
non_mt_genes <- all_genes[all_genes$MCARTA2_LIST==0, ]
test <- subset(mt_genes, TargetP_Score==1, c("Symbol","TargetP_Score","MCARTA2.0_score"))
## notex that subset eliminates N/A values, while the straight command does not!
## however, subset is recommended for interactive use only, so I implicitly check for N/A
test <- mt_genes[mt_genes$TargetP_Score==1&!is.na(mt_genes$TargetP_Score), c("Symbol","TargetP_Score","MCARTA2.0_score")]
str(test)
summary(test)
table(non_mt_genes$TargetP_Score)
51+337+496 ##884 non-mitochondrial genes within cat 1-3
table(mt_genes$TargetP_Score)
193+278+128 ##599 mitochondrial genes within cat 1-3
193/(193+51)
599/(599+884) ##fraction of mitochondrial proteins in categories 1-3
599/1158 
barplot(table(non_mt_genes$TargetP_Score),ylim=c(0,800))
barplot(table(mt_genes$TargetP_Score),ylim=c(0,800))
1158/19244
mouse_genes <- read.csv(file="Mouse.MitoCarta2.0.csv",header=TRUE)
non_mt_mouse_genes <- mouse_genes[mouse_genes$MCARTA2_LIST==0, ]
mt_mouse_genes <- mouse_genes[mouse_genes$MCARTA2_LIST==1, ]
table(non_mt_mouse_genes$TargetP_Score)
table(mt_mouse_genes$TargetP_Score)

49+265+403
213+272+123 
library(ggplot2)
ggplot(mt_genes) + geom_point(aes(x=TargetP_Score, y=MCARTA2.0_score), alpha=0.1) + xlab("TargetP score") + ylab("MCarta score")
#ggplot2 works exclusively with dataframes!
#ggplot2 has many geoms (e.g., geom_line(), geom_bar(), geom_density(), geom_boxplot(), etc.)
ggplot(mt_genes) + geom_boxplot(aes(x=TargetP_Score, y=MCARTA2.0_score, group=TargetP_Score))
####Trying to get genes from Ensembl
####Playing with BioMart
##http://www.ensembl.org/info/data/biomart/biomart_r_package.html
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)
listMarts()
listEnsembl()
ensembl = useMart("ensembl")
listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)
entrez = c("673", "7157", "837")
human_mt_coding <- getSequence(id=Human.MitoCarta2.0$HumanGeneID, type="entrezgene", seqType = "coding", mart=ensembl) ### slow
human_mt_coding2 <- getSequence(id=Human.MitoCarta2.0$HumanGeneID, type="entrezgene", seqType = "coding", mart=refseq) ### slow

check <- getSequence(id=c(1537, 6390), type="entrezgene", seqType = "coding", mart=ensembl)
exportFASTA(human_mt_coding, file="human_mt_coding.fasta")
human_mt_coding <- human_mt_coding[!human_mt_coding$coding == "Sequence unavailable", ]
human_mt_coding <- human_mt_coding[nchar(human_mt_coding$coding)>=300, ]
str(human_mt_coding)
exportFASTA(human_mt_coding, file="human_mt_coding.fasta")

library("seqinr")
install.packages("seqinr")
Human.MitoCarta2.0[as.vector(Human.MitoCarta2.0$Symbol)=="CYC1", ]
as.vector(Human.MitoCarta2.0$Symbol)
exportFASTA(check, file="test.fasta")

