#' QC_gene removes genes that have TPM values < 1 in more than 95 percent of cells in each group
#' 
#' @param rpkmFile a tab delimited file containing expression data (TPM, CPM, FPKM, etc), columns are cells and rows are genes.
#' @param rpkmQCFile resulting file after QC
#' @param sampleFile a tab delimited file containing sample annotations with two columns. The first column indicates SampleID, the second column indicates GroupID.
#' @param threshold the cutoff of percentage of cells in which the given gene is not expressed.Default is 0.05.
#' @param method keep genes whose rpkm values are not less than 1 in at least 5 percent of the cells in every group, method="any": keep genes whose rpkm values are not less than 1 in at least 5 percent of the cells in any group. Default is 'any'.
#' @export
#' @examples
#' \dontrun{
#' QC_gene(rpkmFile = "TPM_MDP_CDP_preDC_Mar2015_noOutlier.txt",
#'         rpkmQCFile = "TPM_MDP_CDP_preDC_Mar2015_noOutlier_geneQC0.05anyGroup.txt",
#'         sampleFile = "sample_MDP_CDP_preDC_Mar2015_noOutlier.txt",
#'         threshold = 0.05,method = "any")
#' }
QC_gene <- function(rpkmFile="TPM_monocyte_Mar2015_noOutlier.txt",
                    rpkmQCFile="TPM_monocyte_Mar2015_noOutlier_geneQC0.05perGroup.txt",
                    sampleFile="monocyte_sample_Mar2015_noOutlier.txt",threshold=0.05,method="any"){
  
  rpkm <- read.table(rpkmFile,sep="\t",header=T,check.names=FALSE,row.names=1)
  

  if(is.null(sampleFile)){
     rpkm_1 <- rpkm > 1
     n_exprs <- apply(rpkm_1,1,sum)  
     pass <- n_exprs > ceiling(ncol(rpkm)*threshold)
     rpkm_qc <- rpkm[pass,]     
  }else{     
    sample <- read.table(sampleFile,sep="\t",header=T)  
    groups <- as.character(unique(sample[,2]))
    pass <- matrix(ncol=length(groups),nrow=nrow(rpkm))
    rownames(pass) <- row.names(rpkm)
    colnames(pass) <- groups
    for(i in 1:length(groups)){
      group_i <- groups[i]
      sample_i <- sample[sample[,2]==group_i,]
      rpkm_i <- rpkm[,as.character(sample_i[,1])]
      rpkm_i_1 <- rpkm_i > 1
      n_exprs_i <- apply(rpkm_i_1,1,sum)
      pass[,i] <- n_exprs_i > ceiling(ncol(rpkm_i)*threshold)     
    }
    pass_sum <- apply(pass,1,sum)
    if(method=="every"){
      rpkm_qc <- rpkm[pass_sum==length(groups),]
    }else{
      rpkm_qc <- rpkm[pass_sum>=1,]
    }
  }
  write.table(rpkm_qc,rpkmQCFile,sep="\t",col.names = NA)
}
