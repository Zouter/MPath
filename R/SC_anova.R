#' SC_anova run ANOVA test to identify differentially expressed genes
#' @param inputfile a tab delimied file containing expression values (TPM, CPM, FPKM, etc). Columns are cells and rows are genes.
#' @param targetfile a tab delimited file indicating cell annoation with two columns. The first column indicates cell ID, the rest columns indicates cell annotations.
#' @param iflog2 a boolean value to indicate whether the expression values will be log2 transformed.
#' @param p_threshold the cutoff of p values for DEGs.
#' @param factor column name of targetfile, indicating which column will be used as cell annotation for comparison.
#' @param baseName prefix name of resulting files.
#' @importFrom ez ezANOVA
#' @importFrom plyr ldply
#' @export
#' @examples
#' \dontrun{
#' inputfile="GSE52529_fpkm_matrix_nooutliers_geneQC0.05anyGroup.txt";
#' targetfile="GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG_landmark_cluster.txt";
#' baseName="GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG_landmark_cluster";
#' SC_anova(inputfile = inputfile,
#'          targetfile = targetfile,
#'          iflog2 = TRUE,p_threshold=0.05,factor="landmark_cluster",
#'          baseName = baseName)
#' }
SC_anova <- function(inputfile,targetfile,iflog2,p_threshold,factor,baseName){
  target <- read.table(targetfile,header=T,sep="\t")
  data <- read.table(inputfile,row.names=1,check.names=FALSE,header=T,sep="\t")
  
  data <- data[,as.character(target[,1])]
  
  nonlog_data <- data
  if(iflog2){
    data <- apply(data,c(1,2),function(x) if(x>1) log2(x) else 0)
  }
  
  rowmax <- apply(data,1,max)
  data <- data[which(rowmax>0),]
  
  iv <- as.factor(target[,factor])
  
  onewayanova.single <- function(exprs){
    m <- data.frame(iv,exprs)
    m$id <- as.factor(rownames(m))
    ezANOVA(data=m,dv=exprs,wid=id,between=iv,return_aov=T)
  }
  
  anovares <- apply(data,1,onewayanova.single)
  
  ########### factor deg
  pval <- data.frame(lapply(anovares,function(x){x$ANOVA[,"p"]}),check.names=FALSE)
  pval <- t(pval)
  BH.adjusted.ANOVA.pval <- p.adjust(pval, method = "BH")
  names(BH.adjusted.ANOVA.pval) <- rownames(pval)
  deg <- data.frame(BH.adjusted.ANOVA.pval[BH.adjusted.ANOVA.pval<p_threshold])
  colnames(deg) <- "BH adjusted ANOVA pval"
  write.table(deg,file=paste(baseName,"anova_",factor,".txt",sep=""),sep="\t",col.names=NA)
  write.table(nonlog_data[as.character(rownames(deg)),],paste(baseName,"_anova_sigGenes.txt",sep=""),sep="\t",col.names=NA) 
  
  ########### post hoc tukey test
  tukeyres <- ldply(anovares,function(x){
    o <- data.frame()
    tr <- TukeyHSD(x$aov)
    for(term in names(tr)){
      o1 <- as.data.frame(tr[[term]])
      o1$comparison <- rownames(o1)
      o <- rbind(o, o1)
    }
    return(o)
  })
  write.table(tukeyres,file=paste(baseName,"anova_tukey.txt",sep=""),sep="\t",row.names=FALSE)
  
  ########## merging files
  tukey <- read.table(paste(baseName,"anova_tukey.txt",sep=""),header=T,sep="\t")
  colnames(tukey)[1] <- "gene"
  
  anova_factor <- read.table(paste(baseName,"anova_",factor,".txt",sep=""),header=T,sep="\t")
  colnames(anova_factor) <- c("gene",paste("anova_",factor,"_pVal",sep=""))
  
  anova_tukey <- merge(tukey,anova_factor,by="gene",all=TRUE)
  
  anova_sig <- !is.na(anova_tukey[,paste("anova_",factor,"_pVal",sep="")])
  anova_tukey_sig <- anova_tukey[anova_tukey$p.adj<0.05 & anova_sig,]
  write.table(anova_tukey_sig,paste(baseName,"anova_tukey_sig.txt",sep=""),sep="\t",row.names=FALSE)  
}
