source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")

# extract gene id from the data
# set path
mainDir <- setwd("C:/Users/Liuhua/Dropbox/RNA_seq/GeneSetTest")
file.list <- dir(file.path(mainDir))
keep <- grepl("genes.results", file.list)
file.temp <- file.list[keep]
remove(file.list,keep,mainDir)
# read file
gene.id1 <- read.table(file.temp, header = TRUE,as.is=FALSE)$gene_id
# extract Entrez Gene ID and convert to characters
gene.id1 <- lapply(gene.id1,toString)
for (i in 1:length(gene.id1)){
  gene.id1[i] <- as.character(strsplit(gene.id1[[i]],split="|",fixed=TRUE)[[1]][2])
}
gene.id <- c(NULL,length(gene.id1))
for (i in 1:length(gene.id1)){
  gene.id[i] <- gene.id1[[i]]
}
remove(gene.id1,file.temp)
# done! get Entrez Gene ID stored in gene.id

library(org.Hs.eg.db)
mapped_genes <- mappedkeys(org.Hs.egGO)
sum(gene.id%in%mapped_genes)
# map the Gene ID to GO Term
GOMap <- as.list(org.Hs.egGO[gene.id[gene.id%in%mapped_genes]])
# get the GO Term list map to Gene ID
GO.list <- as.list(org.Hs.egGO2EG)
length(GO.list)


# function of gene set tests for simulated data using SimSeq, use gene.id, GO.list
geneSetTest.sim <- function(data.sim,gene.id,GO.list,B=100000){
  counts.sim <- data.sim$counts
  n <- ncol(counts.sim)
  n1 <- n2 <- n/2
  p <- nrow(counts.sim)
  p.de <- length(data.sim$DE.genes)
  
  de.sim <- which(data.sim$genes.subset%in%data.sim$DE.genes)
  dge.sim <- DGEList(counts=counts.sim)
  dge.sim <- calcNormFactors(dge.sim,method="TMM")
  design.sim <- cbind(c(rep(1,n1),rep(0,n2)),c(rep(0,n1),rep(1,n2)))
  dge.sim <- estimateGLMCommonDisp(dge.sim,design.sim,verbose=TRUE)
  dge.sim <- estimateGLMTrendedDisp(dge.sim,design.sim,min.n=10)
  dge.sim <- estimateGLMTagwiseDisp(dge.sim,design.sim)
  
  gene.subset <- gene.id[data.sim$genes.subset]
  gene.set.list.sim <- list()
  for (i in 1:length(GO.list)){
    GO.temp <- GO.list[i]
    set.temp <- GO.temp[[1]][GO.temp[[1]]%in%gene.subset]
    if (length(set.temp) >= 2){
      list.temp <- vector("list",2)
      list.temp[[1]] <- names(GO.temp)
      list.temp[[2]] <- set.temp
      names(list.temp) <- c("GOTerm", "GeneSet")
      gene.set.list.sim[[length(gene.set.list.sim)+1]] <- list.temp
    } 
  }
  result <- data.frame(matrix(0,nrow=length(gene.set.list.sim),ncol=4))
  colnames(result) <- c("MRPP.EUC.logtrans.pvalue","MRPP.EUC.vstrans.pvalue",
                        "MRPP.EUC.qua1trans.pvalue","MRPP.EUC.qua2trans.pvalue")
  print(length(gene.set.list.sim))
  LS <- colSums(counts.sim)
  
  for(i in 1:length(gene.set.list.sim)){
    disp <- dge.sim$tagwise.dispersion[which(gene.subset%in%gene.set.list.sim[[i]][[2]])]
    
    x <- counts.sim[which(gene.subset%in%gene.set.list.sim[[i]][[2]]),1:(dim(counts.sim)[2]/2)]
    y <- counts.sim[which(gene.subset%in%gene.set.list.sim[[i]][[2]]),
                        ((dim(counts.sim)[2]/2)+1):dim(counts.sim)[2]]

    x.logtrans <- Transform(x,disp,method="Log")
    y.logtrans <- Transform(y,disp,method="Log")
    x.vstrans <- Transform(x,disp,method="VS")
    y.vstrans <- Transform(y,disp,method="VS")
    x.qua1trans <- Transform(x,disp,method="Qua1")
    y.qua1trans <- Transform(y,disp,method="Qua1")
    x.qua2trans <- Transform(x,disp,method="Qua2")
    y.qua2trans <- Transform(y,disp,method="Qua2")
    
    disMatrix.EUC.logtrans <- DistMatrix(x.logtrans,y.logtrans,LS,method="EUC")
    disMatrix.EUC.vstrans <- DistMatrix(x.vstrans,y.vstrans,LS,method="EUC")
    disMatrix.EUC.qua1trans <- DistMatrix(x.qua1trans,y.qua1trans,LS,method="EUC")
    disMatrix.EUC.qua2trans <- DistMatrix(x.qua2trans,y.qua2trans,LS,method="EUC")
    
    result[i,1] <- MRPP_pvalue(disMatrix.EUC.logtrans,n1,B)
    result[i,2] <- MRPP_pvalue(disMatrix.EUC.vstrans,n1,B)
    result[i,3] <- MRPP_pvalue(disMatrix.EUC.qua1trans,n1,B)
    result[i,4] <- MRPP_pvalue(disMatrix.EUC.qua2trans,n1,B)
    if(i%%50 ==0) print(i)
  }
  temp <- vector("list", 4)
  temp[[1]] <- data.sim
  temp[[2]] <- gene.id
  temp[[3]] <- gene.set.list.sim
  temp[[4]] <- result
  names(temp) <- c("data.sim", "gene.id","gene.set.list", "pvalues")
  return(temp)
}


# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("edgeR")
library(edgeR)
library(SimSeq)
library(fdrtool)
library(pROC)
library(Rcpp)
library(RcppArmadillo)

# setwd("/home/allenberify/Dropbox/RNA_seq/codes")
setwd("C:/Users/Liuhua/Dropbox/RNA_seq/codes")
sourceCpp("RNA_SeqFunction.cpp")
source("RNA_SeqFunction1.R")

# load data
data(kidney)
str(kidney)
counts <- kidney$counts
dim(counts)

Cpm <- cpm(counts)
keep <- rowSums(Cpm>1)>=2
counts.keep <- counts[keep,]
gene.id.keep <- gene.id[keep]
dim(counts.keep)

replic <- kidney$replic
treatment <- kidney$treatment
nf <- apply(counts.keep, 2, quantile, 0.75)
# sample sizes to be simulated
n1 <- n2 <- 10
# dimension size to be simulated
p <- 500
# number of differential expressed genes to be simulated
p.de <- 100


data.sim <- SimData(counts = counts.keep, replic = replic, treatment = treatment,
                    sort.method = "paired", k.ind = n1, n.genes = p, n.diff = p.de,
                    norm.factors = nf)

geneSetTest.result <- geneSetTest.sim(data.sim,gene.id.keep,GO.list)

saveRDS(geneSetTest.result,file="geneSetTest.result.rds")

geneSetTest.result <- readRDS(file="geneSetTest.result.rds")

geneSetTest.pvalues <- geneSetTest.result$pvalues
pairs(geneSetTest.pvalues)
