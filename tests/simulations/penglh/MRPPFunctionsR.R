Transform <- function(y,a,method="Log"){
  if (is.vector(y)){
    if(method=="Log"){
      result <- log(y+1)
    }else if(method=="VS"){
      result <- 2/sqrt(a)*asinh(sqrt(y*a))
    }else if(method=="Qua1"){
      result <- pbeta(a*y/(1+a*y),1/3,1/3)
    }else{
      result <- pbeta( a*y/(1+a*y),1/3,1/3)*(y/a*(1+a*y)/a)^(1/6)*beta(1/3,1/3)
    }
    return(result)
  }else{
    n <- dim(y)[2]
    p <- dim(y)[1]
    if(method=="Log"){
      result <- apply(y,2,function(x) log(x+1))
    }else if(method=="VS"){
      result <- apply(y,2,function(x) 2/sqrt(a)*asinh(sqrt(x*a)))
    }else if(method=="Qua1"){
      result <- apply(y,2,function(x) pbeta(a*x/(1+a*x),1/3,1/3))
    }else{
      result <- apply(y,2,function(x) pbeta( a*x/(1+a*x),1/3,1/3)*(x/a*(1+a*x)/a)^(1/6)*beta(1/3,1/3))
    }
    return(result)
  }
}

my_line <- function(x,y,...){
  points(x,y,...)
  abline(a = 0,b = 1,...)
}

roc.pvalue <- function(pvalue,de.sim){
  p <- length(pvalue)
  p.de <- length(de.sim)
  roc1 <- matrix(0,p,2)
  colnames(roc1) <- c("FPR","TPR")
  rank.pvalue <- rank(pvalue,ties.method="random")
  for(i in 1:p){
    cutOff <- i
    set1 <- which(rank.pvalue<=cutOff)
    set2 <- which(rank.pvalue>cutOff)
    roc1[i,1] <- 1-sum(set2%in%c(1:p)[-de.sim])/(p-p.de)
    roc1[i,2] <- sum(set1%in%de.sim)/p.de
  }
  return(roc1)
}

roc.pvalue.ave <- function(pvalue.matrix,de.sim.matrix){
  S <- nrow(pvalue.matrix)
  p <- ncol(pvalue.matrix)
  p.de <- ncol(de.sim.matrix)
  roc1 <- matrix(0,p,2)
  colnames(roc1) <- c("FPR","TPR")
  rank.pvalue.matrix <- pvalue.matrix
  for(s in 1:S){
    rank.pvalue.matrix[s,] <- rank(pvalue.matrix[s,],ties.method="random")
  }
  for(i in 1:p){
    cutOff <- i
    fpr1 <- rep(0,S)
    tpr1 <- rep(0,S)
    for(s in 1:S){
      set1 <- which(rank.pvalue.matrix[s,]<=cutOff)
      set2 <- which(rank.pvalue.matrix[s,]>cutOff)
      fpr1[s] <- 1-sum(set2%in%c(1:p)[-de.sim.matrix[s,]])/(p-p.de)
      tpr1[s] <- sum(set1%in%de.sim.matrix[s,])/p.de
    }
    roc1[i,1] <- mean(fpr1)
    roc1[i,2] <- mean(tpr1)
  }
  return(roc1[order(roc1[,1]),])
}

auc.partial <- function(pvalue.matrix,de.sim.matrix,a,b){
  S <- nrow(pvalue.matrix)
  p <- ncol(pvalue.matrix)
  result <- rep(0,S)
  for(s in 1:S){
    temp <- (c(1:p)%in%de.sim.matrix[s,])*1
    result[s] <- auc(temp,rank(pvalue.matrix[s,]),partial.auc=c(1-a,1-b))
  }
  return(result)
}

fdr.sim <- function(S=100, B=10000, n1=10, n2=10, p=500, p.de=100){
  n <- n1+n2
  counts.list <- vector("list",S)
  de.matrix <- matrix(0,S,p.de)
  edgeR.fdr <- matrix(0,S,p)
  t.test.norm.fdr <- matrix(0,S,p)
  MRPP.PLLR.fdr <- matrix(0,S,p)
  MRPP.EUC.logtrans.fdr <- matrix(0,S,p)
  MRPP.EUC.vstrans.fdr <- matrix(0,S,p)
  MRPP.EUC.qua1trans.fdr <- matrix(0,S,p)
  MRPP.EUC.qua2trans.fdr <- matrix(0,S,p)
  s= 1
  while(s<=S){
    data.sim <- SimData(counts = counts.keep, replic = replic, treatment = treatment,
                        sort.method = "paired", k.ind = n1, n.genes = p, n.diff = p.de,
                        norm.factors = nf)
    de.sim <- which(data.sim$genes.subset%in%data.sim$DE.genes)
    de.matrix[s,] <- de.sim
    
    counts.sim <- data.sim$counts
    counts.list[[s]] <- counts.sim
    
    if(min(rowSums(counts.sim))==0){
      cat("there is a row with all zero counts!\n")
    }else{
      dge.sim <- DGEList(counts=counts.sim)
      dge.sim <- calcNormFactors(dge.sim,method="TMM")
      design.sim <- cbind(c(rep(1,n1),rep(0,n2)),c(rep(0,n1),rep(1,n2)))
      dge.sim <- estimateGLMCommonDisp(dge.sim,design.sim,verbose=TRUE)
      dge.sim <- estimateGLMTrendedDisp(dge.sim,design.sim,min.n=10)
      dge.sim <- estimateGLMTagwiseDisp(dge.sim,design.sim)
  
      dge.sim.fit <- glmFit(dge.sim,design.sim)
      dge.sim.de <- glmLRT(dge.sim.fit,contrast=c(1,-1))
      tp.sim <- topTags(dge.sim.de,sort="none",n=Inf)
      edgeR.fdr[s,] <- tp.sim$table$FDR
#      cat("edgeR",sim.edgeR.fdr[s,],"\n")
      counts.sim.norm <- counts.sim
      for (j in 1:n){
        counts.sim.norm[,j] <- counts.sim.norm[,j]*mean(dge.sim$sample$lib.size)/(dge.sim$sample$norm.factors[j]*dge.sim$sample$lib.size[j])
      }
      
      t.test.norm.pvalue <- apply(counts.sim.norm,1,function(x) 
        t.test(log(x[1:(length(x)/2)]+1),log(x[(length(x)/2+1):length(x)]+1),var.equal=FALSE)$p.value)
#      cat("t.test",t.test.norm.pvalue,"\n")
      t.test.norm.fdr[s,] <- fdrtool(t.test.norm.pvalue,statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
      
      # use the tagwise estimated dispersion from edgeR
      disp <- dge.sim$tagwise.dispersion
      x <- counts.sim.norm[,1:n1]
      y <- counts.sim.norm[,(n1+1):n]
      x.PLLR <- counts.sim[,1:n1]
      y.PLLR <- counts.sim[,(n1+1):n]
      x.logtrans <- Transform(x,disp,method="Log")
      y.logtrans <- Transform(y,disp,method="Log")
      x.vstrans <- Transform(x,disp,method="VS")
      y.vstrans <- Transform(y,disp,method="VS")
      x.qua1trans <- Transform(x,disp,method="Qua1")
      y.qua1trans <- Transform(y,disp,method="Qua1")
      x.qua2trans <- Transform(x,disp,method="Qua2")
      y.qua2trans <- Transform(y,disp,method="Qua2")
      
      MRPP.EUC.logtrans.pvalue <- MRPP.PLLR.pvalue <- rep(0,p)
      MRPP.EUC.vstrans.pvalue <- MRPP.EUC.qua1trans.pvalue <- MRPP.EUC.qua2trans.pvalue <- rep(0,p)
      LS <- colSums(counts.sim)
      for(i in 1:p){
        disMatrix.PLLR <- DistMatrix(x.PLLR[i,],y.PLLR[i,],LS,method="PLLR")
        disMatrix.EUC.logtrans <- DistMatrix(x.logtrans[i,],y.logtrans[i,],LS,method="EUC")
        disMatrix.EUC.vstrans <- DistMatrix(x.vstrans[i,],y.vstrans[i,],LS,method="EUC")
        disMatrix.EUC.qua1trans <- DistMatrix(x.qua1trans[i,],y.qua1trans[i,],LS,method="EUC")
        disMatrix.EUC.qua2trans <- DistMatrix(x.qua2trans[i,],y.qua2trans[i,],LS,method="EUC")
        
        MRPP.PLLR.pvalue[i] <- MRPP_pvalue(disMatrix.PLLR,n1,B)
        MRPP.EUC.logtrans.pvalue[i] <- MRPP_pvalue(disMatrix.EUC.logtrans,n1,B)
        MRPP.EUC.vstrans.pvalue[i] <- MRPP_pvalue(disMatrix.EUC.vstrans,n1,B)
        MRPP.EUC.qua1trans.pvalue[i] <- MRPP_pvalue(disMatrix.EUC.qua1trans,n1,B)
        MRPP.EUC.qua2trans.pvalue[i] <- MRPP_pvalue(disMatrix.EUC.qua2trans,n1,B)
      }
      MRPP.PLLR.fdr[s,] <- fdrtool(MRPP.PLLR.pvalue,statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
      MRPP.EUC.logtrans.fdr[s,] <- 
        fdrtool(MRPP.EUC.logtrans.pvalue,statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
      MRPP.EUC.vstrans.fdr[s,] <- 
        fdrtool(MRPP.EUC.vstrans.pvalue,statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
      MRPP.EUC.qua1trans.fdr[s,] <- 
        fdrtool(MRPP.EUC.qua1trans.pvalue,statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
      MRPP.EUC.qua2trans.fdr[s,] <- 
        fdrtool(MRPP.EUC.qua2trans.pvalue,statistic="pvalue",plot=FALSE,verbose=FALSE)$lfdr
      
      s <- s+1
    }
    cat(s,"\n")
  }
  return(list("counts.list"=counts.list,"de.matrix"=de.matrix,"edgeR.fdr"=edgeR.fdr,
              "t.test.norm.fdr"=t.test.norm.fdr,
              "MRPP.PLLR.fdr"=MRPP.PLLR.fdr,"MRPP.EUC.logtrans.fdr"=MRPP.EUC.logtrans.fdr,
              "MRPP.EUC.vstrans.fdr"=MRPP.EUC.vstrans.fdr,
              "MRPP.EUC.qua1trans.fdr"=MRPP.EUC.qua1trans.fdr,
              "MRPP.EUC.qua2trans.fdr"=MRPP.EUC.qua2trans.fdr))
}

pvalue.sim <- function(S=100, B=10000, n1=10, n2=10, p=500, p.de=100){
  n <- n1+n2
  counts.list <- vector("list",S)
  de.matrix <- matrix(0,S,p.de)
  limma.pvalue <- matrix(0,S,p)
  edgeR.pvalue <- matrix(0,S,p)
  DESeq.pvalue <- matrix(0,S,p)
  t.test.pvalue <- matrix(0,S,p)
  t.test.norm.pvalue <- matrix(0,S,p)
  MRPP.PLLR.pvalue <- matrix(0,S,p)
  MRPP.EUC.logtrans.pvalue <- matrix(0,S,p)
  MRPP.EUC.vstrans.pvalue <- matrix(0,S,p)
  MRPP.EUC.qua1trans.pvalue <- matrix(0,S,p)
  MRPP.EUC.qua2trans.pvalue <- matrix(0,S,p)
  MRPP.EUC.voomtrans.pvalue <- matrix(0,S,p)
  t.test.vstrans.pvalue <- matrix(0,S,p)
  t.test.qua1trans.pvalue <- matrix(0,S,p)
  t.test.qua2trans.pvalue <- matrix(0,S,p)
  t.test.voomtrans.pvalue <- matrix(0,S,p)
  s= 1
  while(s<=S){
    data.sim <- SimData(counts = counts.keep, replic = replic, treatment = treatment,
                        sort.method = "paired", k.ind = n1, n.genes = p, n.diff = p.de,
                        norm.factors = nf)
    de.sim <- which(data.sim$genes.subset%in%data.sim$DE.genes)
    de.matrix[s,] <- de.sim
    
    counts.sim <- data.sim$counts
    counts.list[[s]] <- counts.sim
    
    if(min(rowSums(counts.sim))==0){
      cat("there is a row with all zero counts!\n")
    }else{
      dge.sim <- DGEList(counts=counts.sim)
      dge.sim <- calcNormFactors(dge.sim,method="TMM")
      
      design.limma.sim <- cbind(rep(1,n),c(rep(0,n1),rep(1,n2)))
      voom.trans <- voom(dge.sim,design.limma.sim,plot=FALSE)
      counts.voomtrans <- voom.trans$E
      limma.fit <- lmFit(voom.trans,design.limma.sim)
      limma.fit <- eBayes(limma.fit,trend=TRUE)
      limma.tp.sim <- topTable(limma.fit,sort="none",coef=ncol(design.limma.sim),n=Inf)
      limma.pvalue[s,] <- limma.tp.sim$P.Value
      
      design.sim <- cbind(c(rep(1,n1),rep(0,n2)),c(rep(0,n1),rep(1,n2)))
      dge.sim <- estimateGLMCommonDisp(dge.sim,design.sim,verbose=TRUE)
      dge.sim <- estimateGLMTrendedDisp(dge.sim,design.sim,min.n=10)
      dge.sim <- estimateGLMTagwiseDisp(dge.sim,design.sim)
      
      dge.sim.fit <- glmFit(dge.sim,design.sim)
      dge.sim.de <- glmLRT(dge.sim.fit,contrast=c(1,-1))
      tp.sim <- topTags(dge.sim.de,sort="none",n=Inf)
      edgeR.pvalue[s,] <- tp.sim$table$PValue
      
      condition <- c(rep(0,n1),rep(1,n2))
      cds <- newCountDataSet(counts.sim,condition)
      cds <- estimateSizeFactors(cds)
      # sizeFactors(cds)
      # head(counts(cds,normalized=TRUE))
      cds <- estimateDispersions(cds)
      # str(fitInfo(cds))
      # plotDispEsts(cds)
      # head(fData(cds))
      res <- nbinomTest(cds,"0","1")
      # head(res)
      DESeq.pvalue[s,] <- res$pval
      
      t.test.pvalue[s,] <- apply(counts.sim,1,function(x) 
        t.test(log(x[1:(length(x)/2)]+1),log(x[(length(x)/2+1):length(x)]+1),var.equal=FALSE)$p.value)

      counts.sim.norm <- counts.sim
      for (j in 1:n){
        counts.sim.norm[,j] <- counts.sim.norm[,j]*mean(dge.sim$sample$lib.size)/(dge.sim$sample$norm.factors[j]*dge.sim$sample$lib.size[j])
      }
      
      t.test.norm.pvalue[s,] <- apply(counts.sim.norm,1,function(x) 
        t.test(log(x[1:(length(x)/2)]+1),log(x[(length(x)/2+1):length(x)]+1),var.equal=FALSE)$p.value)
      
      # use the tagwise estimated dispersion from edgeR
      disp <- dge.sim$tagwise.dispersion
      x <- counts.sim.norm[,1:n1]
      y <- counts.sim.norm[,(n1+1):n]
      x.PLLR <- counts.sim[,1:n1]
      y.PLLR <- counts.sim[,(n1+1):n]
      x.logtrans <- Transform(x,disp,method="Log")
      y.logtrans <- Transform(y,disp,method="Log")
      x.vstrans <- Transform(x,disp,method="VS")
      y.vstrans <- Transform(y,disp,method="VS")
      x.qua1trans <- Transform(x,disp,method="Qua1")
      y.qua1trans <- Transform(y,disp,method="Qua1")
      x.qua2trans <- Transform(x,disp,method="Qua2")
      y.qua2trans <- Transform(y,disp,method="Qua2")
      x.voomtrans <- counts.voomtrans[,1:n1]
      y.voomtrans <- counts.voomtrans[,(n1+1):n]
      
      t.test.vstrans.pvalue[s,] <- apply(cbind(x.vstrans,y.vstrans),1,function(x) 
        t.test(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)],var.equal=FALSE)$p.value)  
      t.test.qua1trans.pvalue[s,] <- apply(cbind(x.qua1trans,y.qua1trans),1,function(x) 
        t.test(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)],var.equal=FALSE)$p.value)
      t.test.qua2trans.pvalue[s,] <- apply(cbind(x.qua2trans,y.qua2trans),1,function(x) 
        t.test(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)],var.equal=FALSE)$p.value)
      t.test.voomtrans.pvalue[s,] <- apply(counts.voomtrans,1,function(x)
        t.test(x[1:(length(x)/2)],x[(length(x)/2+1):length(x)],var.equal=FALSE)$p.value)
            
      LS <- colSums(counts.sim)
      for(i in 1:p){
        disMatrix.PLLR <- DistMatrix(x.PLLR[i,],y.PLLR[i,],LS,method="PLLR")
        disMatrix.EUC.logtrans <- DistMatrix(x.logtrans[i,],y.logtrans[i,],LS,method="EUC")
        disMatrix.EUC.vstrans <- DistMatrix(x.vstrans[i,],y.vstrans[i,],LS,method="EUC")
        disMatrix.EUC.qua1trans <- DistMatrix(x.qua1trans[i,],y.qua1trans[i,],LS,method="EUC")
        disMatrix.EUC.qua2trans <- DistMatrix(x.qua2trans[i,],y.qua2trans[i,],LS,method="EUC")
        disMatrix.EUC.voomtrans <- DistMatrix(x.voomtrans[i,],y.voomtrans[i,],LS,method="EUC")
        
        MRPP.PLLR.pvalue[s,i] <- MRPP_pvalue(disMatrix.PLLR,n1,B)
        MRPP.EUC.logtrans.pvalue[s,i] <- MRPP_pvalue(disMatrix.EUC.logtrans,n1,B)
        MRPP.EUC.vstrans.pvalue[s,i] <- MRPP_pvalue(disMatrix.EUC.vstrans,n1,B)
        MRPP.EUC.qua1trans.pvalue[s,i] <- MRPP_pvalue(disMatrix.EUC.qua1trans,n1,B)
        MRPP.EUC.qua2trans.pvalue[s,i] <- MRPP_pvalue(disMatrix.EUC.qua2trans,n1,B)
        MRPP.EUC.voomtrans.pvalue[s,i] <- MRPP_pvalue(disMatrix.EUC.voomtrans,n1,B)
      }
      
      s <- s+1
    }
    cat(s,"\n")
  }
  return(list("counts.list"=counts.list,"de.matrix"=de.matrix,
              "limma.pvalue"=limma.pvalue,
              "edgeR.pvalue"=edgeR.pvalue,
              "DESeq.pvalue"=DESeq.pvalue,
              "t.test.pvalue"=t.test.pvalue,
              "t.test.norm.pvalue"=t.test.norm.pvalue,
              "t.test.vstrans.pvalue"=t.test.vstrans.pvalue,
              "t.test.qua1trans.pvalue"=t.test.qua1trans.pvalue,
              "t.test.qua2trans.pvalue"=t.test.qua2trans.pvalue,
              "t.test.voomtrans.pvalue"=t.test.voomtrans.pvalue,
              "MRPP.PLLR.pvalue"=MRPP.PLLR.pvalue,
              "MRPP.EUC.logtrans.pvalue"=MRPP.EUC.logtrans.pvalue,
              "MRPP.EUC.vstrans.pvalue"=MRPP.EUC.vstrans.pvalue,
              "MRPP.EUC.qua1trans.pvalue"=MRPP.EUC.qua1trans.pvalue,
              "MRPP.EUC.qua2trans.pvalue"=MRPP.EUC.qua2trans.pvalue,
              "MRPP.EUC.voomtrans.pvalue"=MRPP.EUC.voomtrans.pvalue))
}

fdr.adj <- function(pvalue.matrix){
  p <- ncol(pvalue.matrix)
  S <- nrow(pvalue.matrix)
  fdr.matrix <- pvalue.matrix
  for(s in 1:S){
    temp <- p*pvalue.matrix[s,]/rank(pvalue.matrix[s,],ties.method="random")
    temp[which(temp>1)] <- 1
    fdr.matrix[s,] <- temp
  }
  return(fdr.matrix)
}