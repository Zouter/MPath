#' Mpath: an analysis algorithm that  maps multi-branching single-cell trajectories from single-cell RNA-sequencing data
#'
#' This package provides a new algorithm that maps multi-branching cell developmental pathways 
#' and aligns individual cells along the continuum of developmental trajectories. 
#' Mpath computationally reconstructs cell developmental pathways as a multi-destination journey on a map of connected landmarks 
#' wherein individual cells are placed in order along the paths connecting the landmarks. 
#' To achieve that, it first identifies clusters of cells and designates landmark clusters each defines a discrete cellular state. 
#' Subsequently it identifies and counts cells that are potentially transitioning from one landmark state to the next based on transcriptional similarities. 
#' It then uses the cell counts to infer putative transitions between landmark states giving rise to a state transition network. 
#' After that, Mpath sorts individual cells according to their various stages during transition to resemble the landmark-to-landmark continuum. 
#' Lastly, Mpath detects genes that were differentially expressed along the single-cell trajectories and identifies candidate regulatory markers.
#' 
#' 
#' @examples
#' \dontrun{
#' #### Install and load Mpath package
#' 
#' install.packages("Mpath_1.0.tar.gz",repos = NULL, type="source")
#' library(Mpath)
#' 
#' #################################################
#' ###### Analysis of mouse DC dataset GSE60783 ####
#' #################################################
#' 
#' path <- getwd()
#' setwd(paste(path,"/GSE60783",sep=""))
#' 
#' ##### remove low detection rate genes
#' 
#' rpkmFile="TPM_GSE60783_noOutlier.txt";
#' rpkmQCFile="TPM_GSE60783_noOutlier_geneQC0.05anyGroup.txt";
#' sampleFile="sample_GSE60783_noOutlier.txt";
#' QC_gene(rpkmFile=rpkmFile,
#'         rpkmQCFile=rpkmQCFile,
#'         sampleFile=sampleFile,threshold=0.05,method="any")
#' 
#' ### Mpath using spleenic CD4 vs CD8 DEGs 
#' 
#' rpkmFile = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
#' sampleFile = "sample_GSE60783_noOutlier.txt";
#' find_optimal_cluster_number(rpkmFile = rpkmFile,
#'                             sampleFile = sampleFile,
#'                             min_cluster_num = 7, max_cluster_num = 15,
#'                             diversity_cut = 0.6, size_cut = 0.05)
#'  
#' ### Landmark designation
#' 
#' rpkmFile = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
#' baseName = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG";
#' sampleFile = "sample_GSE60783_noOutlier.txt";
#' landmark_cluster <- landmark_designation(rpkmFile = rpkmFile,
#'                                          baseName = baseName,
#'                                          sampleFile = sampleFile,
#'                                          method = "diversity_size",
#'                                          numcluster = 11, diversity_cut=0.6,
#'                                          size_cut=0.05)
#' 
#' ### Plot hierachical clustering 
#' 
#' dataFile = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
#' SC_hc_colorCode(dataFile = dataFile,
#'                 cuttree_k = 11,
#'                 sampleFile= "sample_GSE60783_noOutlier.txt",
#'                 width = 22, height = 10, iflog2 = TRUE,
#'                 colorPalette = c("red","green","blue"))
#' 
#' ### Construct weighted neighborhood network
#' 
#' exprs = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
#' baseName = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG";
#' neighbor_network <- build_network(exprs = exprs,
#'                                   landmark_cluster = landmark_cluster,
#'                                   baseName = baseName)
#' 
#' ### TrimNet: trim edges of lower weights
#' 
#' baseName = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG";
#' trimmed_net <- trim_net(neighbor_network,textSize=30,
#'                         baseName = baseName,
#'                         method = "mst")
#' 
#' ### plot trimmed net and color-code the nodes by gene expression
#' 
#' rpkmFile="TPM_GSE60783_noOutlier.txt";
#' lmFile="TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG_landmark_cluster.txt";
#' color_code_node_2(networkFile=trimmed_net,
#'                   rpkmFile=rpkmFile,
#'                   lmFile=lmFile,
#'                   geneName=c("Irf8","Id2","Batf3"),
#'                   baseName="cDC1_marker",
#'                   seed=NULL)
#' 
#' ### Re-order the cells on the path connecting landmark 
#' ### "CDP_2","CDP_1","PreDC_9","PreDC_3"
#' 
#' exprs = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
#' ccFile = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG_landmark_cluster.txt";
#' order <- nbor_order(exprs = exprs,
#'                     ccFile = ccFile, 
#'                     lm_order = c("CDP_2","CDP_1","preDC_9","preDC_3"),
#'                     if_bb_only=FALSE,
#'                     method=1)
#' 
#' ### identify genes that changed along the cell re-ordering
#' 
#' deg <- vgam_deg(exprs = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup.txt",
#'                 order = order,
#'                 lm_order = c("CDP_2","CDP_1","preDC_9","preDC_3"),
#'                 min_expr=1,
#'                 p_threshold=0.05)
#' 
#' ### plot heatmap of genes that changed along the cell re-ordering
#' 
#' heatmap_nbor(exprs = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup.txt",                        
#'              cell_order = "CDP_2_CDP_1_preDC_9_preDC_3_order.txt", 
#'              plot_genes = "CDP_2_CDP_1_preDC_9_preDC_3_vgam_deg0.05.txt",
#'              cell_annotation = "sample_GSE60783_noOutlier.txt", 
#'              num_gene_cluster = 6,
#'              hm_height = 15, hm_width = 10,
#'              baseName = "CDP_2_CDP_1_preDC_9_preDC_3_order_vgam_deg0.05")
#'              
#' ###################################################
#' ### Analysis of human myoblast dataset GSE52529 ###
#' ###################################################
#' 
#' setwd(paste(path,"/GSE52529",sep=""))
#' 
#' ### remove low detection rate genes
#' 
#' QC_gene(rpkmFile = "GSE52529_fpkm_matrix_nooutliers.txt",
#'         rpkmQCFile = "GSE52529_fpkm_matrix_nooutliers_geneQC0.05anyGroup.txt",
#'         sampleFile = "sample_nooutlier.txt",threshold=0.05,method="any")
#'
#' rpkmFile = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG.txt";
#' sampleFile = "sample_nooutlier.txt";
#' find_optimal_cluster_number(rpkmFile = rpkmFile,
#'                             sampleFile = sampleFile,
#'                             min_cluster_num = 7, max_cluster_num = 18,
#'                             diversity_cut = 0.9, size_cut = 0.05)
#'
#' rpkmFile = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG.txt";
#' baseName = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG";
#' landmark_cluster <- landmark_designation(rpkmFile = rpkmFile,
#'                                          baseName = baseName,
#'                                          sampleFile = "sample_nooutlier.txt",
#'                                          method = "diversity_size",
#'                                          numcluster = 14, diversity_cut=0.9, 
#'                                          size_cut=0.05)
#'
#' SC_hc_colorCode(dataFile = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG.txt",
#'                 cuttree_k = 14,
#'                 sampleFile= "sample_nooutlier.txt",
#'                 width = 22, height = 10, iflog2 = TRUE)
#'
#' ### Construct weighted neighborhood network
#' 
#' exprs = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG.txt";
#' neighbor_network <- build_network(exprs = exprs,
#'                                   landmark_cluster = landmark_cluster)
#'
#'### TrimNet: trim edges of lower weights
#'
#' trimmed_net <- trim_net(neighbor_network,textSize=30,
#'                        baseName = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG",
#'                        method = "mst")
#'
#' ### Color code the landmarks with marker expression
#' 
#' networkFile="GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG_state_transition_mst.txt";
#' lmFile = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG_landmark_cluster.txt";
#' rpkmFile = "GSE52529_fpkm_matrix_nooutliers_geneSymbol.txt";
#' geneName=c("SPHK1","PBX1","XBP1","ZIC1","MZF1","CUX1","ARID5B","POU2F1","CDK1","MYOG");
#' color_code_node_2(networkFile = networkFile,
#'                   rpkmFile = rpkmFile,
#'                   lmFile = lmFile,
#'                   geneName = geneName,
#'                   baseName = "Marker",
#'                   seed=3)
#'
#'
#' ### path 1: "T0_2","T0_1","T24_8","T48_10","T72_13"
#' 
#' ccFile = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG_landmark_cluster.txt";
#' order <- nbor_order(exprs = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG.txt",
#'                     ccFile = ccFile, 
#'                     lm_order = c("T0_2","T0_1","T24_8","T48_10","T72_13"),
#'                     if_bb_only=TRUE,
#'                     method=1)
#'
#' deg <- vgam_deg(exprs = "GSE52529_fpkm_matrix_nooutliers_geneQC0.05anyGroup.txt",
#'                order = order,
#'                lm_order = c("T0_2","T0_1","T24_8","T48_10","T72_13"),
#'                min_expr=1,
#'                p_threshold=0.05)
#'
#' heatmap_nbor(exprs = "GSE52529_fpkm_matrix_nooutliers_geneQC0.05anyGroup.txt",                        
#'             cell_order = order,
#'             plot_genes = row.names(deg),
#'             cell_annotation = "sample_nooutlier.txt", 
#'             num_gene_cluster = 6,
#'             hm_height = 15, hm_width = 10,
#'             baseName = "GSE52529_path1_order_vgam_deg0.05")
#'
#' ### path 2: "T0_2","T0_1","T24_7","T48_4","T72_14"
#' ccFile = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG_landmark_cluster.txt";
#' order <- nbor_order(exprs = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG.txt",
#'                     ccFile = ccFile, 
#'                     lm_order = c("T0_2","T0_1","T24_7","T48_4","T72_14"),
#'                     if_bb_only = TRUE,
#'                     method=1)
#'
#' deg <- vgam_deg(exprs = "GSE52529_fpkm_matrix_nooutliers_geneQC0.05anyGroup.txt",
#'                order = order,
#'                lm_order = c("T0_2","T0_1","T24_7","T48_4","T72_14"),
#'                min_expr=1,
#'                p_threshold=0.05)
#'
#' heatmap_nbor(exprs = "GSE52529_fpkm_matrix_nooutliers_geneQC0.05anyGroup.txt",                        
#'             cell_order = order,
#'             plot_genes = row.names(deg),
#'             cell_annotation = "sample_nooutlier.txt", 
#'             num_gene_cluster = 6,
#'             hm_height = 15, hm_width = 10,
#'             baseName = "GSE52529_path2_order_vgam_deg0.05")
#'
#' ### Generate Figure 5
#' 
#' ccFile = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG_landmark_cluster.txt";
#' order1 <- nbor_order(exprs = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG.txt",
#'                      ccFile = ccFile, 
#'                      lm_order = c("T0_2","T0_1","T24_8","T48_10","T72_13"),
#'                      if_bb_only=TRUE,
#'                      method=1)
#' 
#' ccFile = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG_landmark_cluster.txt";
#' order2 <- nbor_order(exprs = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG.txt",
#'                      ccFile = ccFile, 
#'                      lm_order = c("T0_2","T0_1","T24_7","T48_4","T72_14"),
#'                      if_bb_only=TRUE,
#'                      method=1)
#'
#' deg1 <- read.table("T0_2_T0_1_T24_8_T48_10_T72_13_vgam_deg0.05.txt",sep="\t",header=T)
#' deg2 <- read.table("T0_2_T0_1_T24_7_T48_4_T72_14_vgam_deg0.05.txt",sep="\t",header=T)
#' deg <- unique(c(as.character(deg1[,1]),as.character(deg2[,1])))
#'
#' heatmap_nbor(exprs = "GSE52529_fpkm_matrix_nooutliers_ANOVA_p0.05_DEG.txt",                        
#'             cell_order = c(order1,order2), 
#'             plot_genes = deg,
#'             cell_annotation = "sample_nooutlier.txt", 
#'             num_gene_cluster = 7,
#'             hm_height = 15, hm_width = 10,
#'             baseName = "Path12_method1orderedbackbone_progression_heatmap", 
#'             n_linechart = list(order1,order2))
#' }
#'              
#' @docType package
#' @name Mpath-package
NULL


#' heatmap_nbor plot heatmap of gene expression
#' 
#' @param exprs a data frame or matrix of expression data(ie. rpkm, TPM, fpkm) containing cells in columns and genes in rows
#' @param cell_order a vector storing the order of cells with cell ID or name, same as appeared in column names of \code{exprs}
#' @param plot_genes a vector storing the genes selected for plot, same as appeared in the row names of \code{exprs}
#' @param cell_annotation a two column data frame or matrix annotating cells(cell ID or name) with cell types
#' @param num_gene_cluster a integer indicating the number of gene clusters to generate by cutting the dendrogram tree, if num_gene_cluster = NULL, no gene clusters will be generated
#' @param hm_height an integeter to specify the heatmap height
#' @param hm_width an integeter to specify heatmap width
#' @param baseName a character string to specify the prefix of name of result files
#' @param colorPalette a character vector to specify the color scheme of column side bar that labels cell type of individual cells. Default value is NULL, a default color scheme will be deployed. Alternatively users can specfiy their desired color paletter, for examples, colorPalette=c("red","green","blue","black","orange","purple","burlywood4") 
#' @param n_linechart a string vector to specify for which cell ordering to plot the line chart. Default value is NULL, all the cells in cell_order will be plot. Alternatively, users could specify, for examples two line charts, by n_linechart = list(c("cell1","cell2"),c("cell3","cell4"))
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples 
#' \dontrun{
#' heatmap_nbor(exprs = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup.txt",                        
#'              cell_order = "CDP_2_CDP_1_preDC_9_preDC_3_order.txt", 
#'              plot_genes = "CDP_2_CDP_1_preDC_9_preDC_3_vgam_deg0.05.txt",
#'              cell_annotation = "sample_GSE60783_noOutlier.txt", 
#'              num_gene_cluster = 6,
#'              hm_height = 15, hm_width = 10,
#'              baseName = "CDP_2_CDP_1_preDC_9_preDC_3_order_vgam_deg0.05") 
#' }
heatmap_nbor <- function(exprs, cell_order, plot_genes, 
                         cell_annotation, num_gene_cluster = 4,
                         hm_height = 10, hm_width = 10, baseName, 
                         colorPalette = NULL, n_linechart=NULL){
  if(is.character(exprs)){
    if(file.exists(exprs))
      exprs <- read.table(exprs, sep = "\t", header = T, row.names = 1, check.names = FALSE)
  }  
  
  if(is.character(cell_order)){
    if(file.exists(cell_order))
      cell_order <- as.character(read.table(cell_order)[,1])
  }
  
  if(is.character(plot_genes)){
    if(file.exists(plot_genes))
      plot_genes <- as.character(read.table(plot_genes,sep="\t",header=T)[,1])
  }
  plot_genes <- unique(plot_genes)
  
  if(is.character(cell_annotation)){
    if(file.exists(cell_annotation))
      cell_annotation <- read.table(cell_annotation, header = TRUE, row.names = 1)
  }
  
  if(is.null(colorPalette)){
    colorPalette <- brewer.pal(length(unique(cell_annotation[,1])),"Set2")
  }
  
  group_color <- colorPalette[as.factor(cell_annotation[cell_order,1])]
  

  log2exprs <- apply(exprs,c(1,2),function(x) if(x>1) log2(x) else 0)
  if(!is.null(plot_genes)){
     log2exprs <- log2exprs[plot_genes[plot_genes %in% row.names(log2exprs)],]
  }
  if(!is.null(cell_order)){
    log2exprs <- log2exprs[,cell_order]    
  }
  
  log2exprs <- log2exprs[rowSums(log2exprs)>0,]
  log2exprs_z <- apply(log2exprs,1,scale)
  rownames(log2exprs_z) <- colnames(log2exprs)
  log2exprs_z <- t(log2exprs_z)

  hclust2 <- function(x, method="ward")  hclust(x, method="ward")
  dist2 <- function(x)  as.dist(1-cor(t(x), method="pearson"))

  if(is.null(num_gene_cluster)){
    #pdf(paste(baseName,"_heatmap.pdf",sep=""),height=hm_height,width=hm_width)
    par(oma = c(5, 1, 1, 1))
    hm <- heatmap.2(as.matrix(log2exprs_z),col=bluered,breaks=seq(-2,2,0.1),trace="none",density="none",scale="none",margins = c(8,8),ColSideColors=group_color,keysize=1,dendrogram="row",Colv=FALSE,distfun=dist2,hclustfun=hclust2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom",legend=unique(cell_annotation[cell_order,1]),pch=15,col=unique(group_color), bty = "n")
    #dev.off()
  }else{
    hm <- heatmap.2(as.matrix(log2exprs_z),col=bluered,breaks=seq(-2,2,0.1),trace="none",density="none",scale="none",margins = c(8,8),ColSideColors=group_color,keysize=1,dendrogram="row",Colv=FALSE,distfun=dist2,hclustfun=hclust2)
    gene_cluster <- cutree(as.hclust(hm$rowDendrogram),k=num_gene_cluster)
    my_palette <- brewer.pal(num_gene_cluster,"Set2")
    gene_color <- my_palette[gene_cluster]
    #pdf(paste(baseName,"_heatmap.pdf",sep=""),height=hm_height,width=hm_width) 
    par(oma = c(5, 1, 1, 1))
    heatmap.2(as.matrix(log2exprs_z),col=bluered,breaks=seq(-2,2,0.1),trace="none",density="none",scale="none",margins = c(8,8),RowSideColors=gene_color,ColSideColors=group_color,keysize=1,dendrogram="row",Colv=FALSE,distfun=dist2,hclustfun=hclust2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottomleft",legend=unique(cell_annotation[cell_order,1]),pch=15,col=unique(group_color),bty = "n")
    legend("bottomright",legend=unique(gene_cluster),pch=15,col=unique(gene_color), bty = "n")
    #dev.off()
    
    log2exprs_gene_cluster <- aggregate(log2exprs,list(gene_cluster=gene_cluster),mean)
    if(is.null(n_linechart)){
       #pdf(paste(baseName,"_linechart.#pdf",sep=""))
       par(oma = c(2,2,2,2))
       plot(lowess(t(log2exprs_gene_cluster[log2exprs_gene_cluster$gene_cluster==1,])~c(1:ncol(log2exprs_gene_cluster)),f=1/5), col=my_palette[1], lwd=3, type = "l", bty = "n", ylim=c(min(log2exprs_gene_cluster[,cell_order]),max(log2exprs_gene_cluster[,cell_order])), xlab="Cell no in order", ylab="log2TPM")
       for(i in 2:num_gene_cluster){
           lines(lowess(t(log2exprs_gene_cluster[log2exprs_gene_cluster$gene_cluster==i,])~c(1:ncol(log2exprs_gene_cluster)),f=1/5),col=my_palette[i],lwd=3)
       } 
       par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
       plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
       legend("topright",legend=unique(gene_cluster),pch=15,col=unique(gene_color), bty = "n")    
       #dev.off()
    }else{
       for(k in 1:length(n_linechart)){
         ##pdf(paste(baseName,"_linechart","_order",k,".#pdf",sep=""))
         par(oma = c(2,2,2,2))
         plot(lowess(t(log2exprs_gene_cluster[log2exprs_gene_cluster$gene_cluster==1,n_linechart[[k]]])~c(1:length(n_linechart[[k]])),f=1/5), col=my_palette[1], lwd=3, type = "l", bty = "n", ylim=c(min(log2exprs_gene_cluster[,cell_order]),max(log2exprs_gene_cluster[,cell_order])), xlab="Cell no in order", ylab="log2TPM")
         for(i in 2:num_gene_cluster){
           lines(lowess(t(log2exprs_gene_cluster[log2exprs_gene_cluster$gene_cluster==i,n_linechart[[k]]])~c(1:length(n_linechart[[k]])),f=1/5),col=my_palette[i],lwd=3)
         } 
         par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
         plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
         legend("topright",legend=unique(gene_cluster),pch=15,col=unique(gene_color), bty = "n")    
         #dev.off()
       }
    }
    
    write.table(data.frame(gene=names(gene_cluster),cluster=gene_cluster),paste(baseName,"_gene_cluster.txt",sep=""),sep="\t",row.names=F)
  }
}

#' nbor_order sorts individual cells according to their various stages during transition to resemble the landmark-to-landmark continuum
#' @param exprs a data frame or matrix of expression data(ie. rpkm, TPM, fpkm) or a tab delimited txt file of expression data, containing cells in columns and genes in rows
#' @param ccFile a data frame or matrix of two columns or a tab delimited file of landmark cluster assignment of individual cells. The first column indicates cell ID, the second column indicates the landmark cluster which the cell was assigned to.  
#' @param lm_order a vector of landmark IDs indicating along which path the cells are to be sorted
#' @param if_bb_only a boolean to indicate if only cells on backbone will be sorted. Default is FALSE
#' @param method 1 or 2 to indicate which method to be used for sorting. Default is 1
#' @param writeRes a boolean to indicate whether to save result files
#' @return a vector of  re-orderd cell IDs
#' @export
#' @examples
#' \dontrun{
#' exprs = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
#' ccFile = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG_landmark_cluster.txt";
#' order <- nbor_order(exprs = exprs,
#'                     ccFile = ccFile, 
#'                     lm_order = c("CDP_2","CDP_1","preDC_9","preDC_3"),
#'                     if_bb_only=FALSE,
#'                     method=1)
#' }
nbor_order <- function(exprs, ccFile,
                       lm_order = c("CD115+CDP_1","PreDC_9","PreDC_10","PreDC_11"),
                       if_bb_only = FALSE, method = 1, writeRes = TRUE){
  
  if(is.character(exprs)){
    if(file.exists(exprs))
      exprs <- read.table(exprs, sep = "\t", header = T, row.names = 1, check.names = FALSE)
  }
  log2exprs <- apply(exprs,c(1,2),function(x) if(x>1) log2(x) else 0)
  
  lm <- generate_lm(ccFile,log2exprs)
  nbor_res <- nbor(log2exprs, lm = lm, lm_order = lm_order, if_bb_only = if_bb_only, method = method) 
  
  order <- vector(length=0)
  for(i in 1:length(nbor_res)){
    order <- c(order,nbor_res[[i]]$order)
  } 
  if (writeRes) {
    write.table(order,paste(paste(lm_order,collapse="_"),"_order.txt",sep=""),sep="\t",col.names=F,row.names=F)
  }
  return(order)
}

#' vgam_deg identifies genes that were differentially expressed along the re-ordered single-cell trajectories using vgam
#' @param exprs a data frame or matrix of expression data(ie. rpkm, TPM, fpkm) or a tab delimited txt file of expression data, containing cells in columns and genes in rows
#' @param order a vector of re-ordered cell IDs
#' @param lm_order a vector of landmark IDs indicating along which path the cells are to be sorted
#' @param min_expr a numeric value indicating the minimum TPM value for a gene to be considered as expressed. Default is 1.
#' @param p_threshold p value cutoff for selecting differentially expressed genes.
#' @return deg: a list of differentially expressed genes
#' @export 
#' @examples 
#' \dontrun{
#' deg <- vgam_deg(exprs = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup.txt",
#'                 order = order,
#'                 lm_order = c("CDP_2","CDP_1","preDC_9","preDC_3"),
#'                 min_expr=1,
#'                 p_threshold=0.05)
#' }

vgam_deg <- function(exprs="TPM_GSE60783_noOutlier_geneQC0.05anyGroup.txt",
                     order,lm_order = c("CD115+CDP_2","CD115+CDP_1","PreDC_9","PreDC_3"),
                     min_expr = 1, p_threshold = 0.05){
  
  if(is.character(exprs)){
    if(file.exists(exprs))
      exprs <- read.table(exprs, sep = "\t", header = T, row.names = 1, check.names = FALSE)
  }
  if(length(grep("^Gm[0-9][0-9]*$",row.names(exprs)))>0){
     exprs <- exprs[-grep("^Gm[0-9][0-9]*$",row.names(exprs)),]
  }
  if(length(grep("^n-R5",row.names(exprs)))>0){
     exprs <- exprs[-grep("^n-R5",row.names(exprs)),]
  }
  
  log2exprs <- apply(exprs, c(1,2), function(x) if(x>1) log2(x) else 0)
    
  log2exprs_order <- log2exprs[,order]
  
  p <- apply(log2exprs_order, 1, FUN=function(x) vgam_perGene(x, c(1:length(order)), min_expr))
  p <- p[which(!is.na(p))]
  p.adj<- p.adjust(p,method="BH")
  res <- data.frame(p=p,p.adj=p.adj)
  deg <- res[res$p.adj<0.05,]
  
  write.table(deg,paste(paste(lm_order,collapse="_"),"_vgam_deg",p_threshold,".txt",sep=""),sep="\t",col.names=NA)
  return(deg)
}

#' vgam_perGene determines if a gene is differentially expressed along the re-ordered single-cell trajectories using vgam 
#' @param expr a vector of one gene's expression in different cells (ie. rpkm, TPM, fpkm)
#' @param order a vector of re-ordered cell IDs
#' @param min_expr a numeric value indicating the minimum TPM value for a gene to be considered as expressed. Default is 1.
#' @return pval: p value of significance of the gene being differentially expressed
#' @importFrom VGAM vgam
#' @export
#' @examples
#' \dontrun{
#' p_val <- vgam_perGene(expr,order,min_expr=1)
#' }

vgam_perGene <- function(expr,order,min_expr){
  expr_order <- data.frame(expr=expr,order=order)
  full_model <- fit_fullmodel(expr_order,min_expr)
  reduced_model <- fit_reducedmodel(expr_order,min_expr)
  
  if(is.null(full_model)) return(NA)
  lrt <- lrtest(full_model,reduced_model)
  pval=lrt@Body["Pr(>Chisq)"][2,]
  return(pval)
}

#' @importFrom VGAM vgam
fit_fullmodel <- function(expr_order,min_expr){
  tryCatch({
    full_model <- vgam(expr~s(order,3),data=expr_order,family=tobit(Lower=log2(min_expr), Upper=Inf))
    full_model
  },
  error=function(e){NULL}
  )  
}

#' @importFrom VGAM vgam
fit_reducedmodel <- function(expr_order,min_expr){
  tryCatch({
    reduced_model <- vgam(expr~1,data=expr_order,family=tobit(Lower=log2(min_expr), Upper=Inf))
    reduced_model
  },
  error=function(e){NULL}
  )  
}

#' find_optimal_cluster_number identifies the optimal number of initial cluster number by searching from min_cluster_num to max_cluster_num
#' @param rpkmFile a tab delimited txt file of expression data, containing cells in columns and genes in rows 
#' @param sampleFile a tab delimited txt file of sample annotation with two columns, the first column is cell ID, the second column is group ID
#' @param min_cluster_num minimum number of initial clusters
#' @param max_cluster_num maximum number of initial clusters
#' @param diversity_cut the cutoff value of diversity for differentiating landmark clusters from non-landmark clusters. The diversity of a landmark cluster must be below this cutoff.
#' @param size_cut the cutoff value of size i.e. number of cells for differentiating landmark clusters from non-landmark clusters. The number of cells in a landmark cluster must be greater than this cutoff.
#' @import ggplot2
#' @export
#' @examples 
#' \dontrun{
#' rpkmFile = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
#' sampleFile = "sample_GSE60783_noOutlier.txt";
#' find_optimal_cluster_number(rpkmFile = rpkmFile,
#'                             sampleFile = sampleFile,
#'                             min_cluster_num = 7, max_cluster_num = 15,
#'                             diversity_cut = 0.6, size_cut = 0.05)
#' }
find_optimal_cluster_number <- function(rpkmFile, sampleFile,
                                        min_cluster_num = 7, max_cluster_num = 13,
                                        diversity_cut = 0.6, size_cut = 0.05){
  name <- sub(".txt","",rpkmFile)  
  
  len <- max_cluster_num-min_cluster_num+1
  ncluster_nlm <- data.frame(ncluster=c(min_cluster_num:max_cluster_num),nlm=rep(0,len))
  
  for(i in 1:len){
    numcluster <- ncluster_nlm$ncluster[i]
    baseName=paste(name,"_clusternum",numcluster,"diversity",diversity_cut,"size",size_cut,sep="") 
    res <- landmark_designation(rpkmFile=rpkmFile,baseName=baseName,sampleFile=sampleFile,method="diversity_size",numcluster=numcluster,diversity_cut=diversity_cut,size_cut=size_cut,saveRes=F)  
    ncluster_nlm$nlm[i] <- length(unique(res$landmark_cluster))
  }
  
  myplot <- ggplot(ncluster_nlm, aes(x = ncluster, y = nlm)) + geom_point(size=8) + geom_line(linetype="longdash",size=1) + xlab("Number of clusters") + ylab("Number of landmark clusters") 
  myplot + scale_x_continuous(breaks=ncluster_nlm$ncluster) + scale_y_continuous(breaks=ncluster_nlm$nlm) + theme_bw() + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.text = element_text(size=20), axis.title=element_text(size=25))
  ggsave(paste(name,"_ncluster_vs_nlm.#pdf",sep=""),height=8,width=8)
}

#' landmark_designation clusters cells and determines landmark clusters
#' @param rpkmFile a tab delimited txt file of expression data, containing cells in columns and genes in rows 
#' @param baseName a character string indicating of prefix name of resulting files
#' @param sampleFile a tab delimited txt file of sample annotation with two columns, the first column is cell ID, the second column is group ID
#' @param distMethod the method for calculating dissimilarity between cells. distMethod can be one of "pearson", "kendall", "spearman" or "euclidean". Default is "euclidean".
#' @param method method for distinguishing landmark clusters from non-landmark clusters.method can be "kmeans" or "diversity" or "size" or "diversity_size". When method="diversity", numlm needs to be specified. Default is "diversity_size".
#' @param numcluster number of initial clusters
#' @param diversity_cut the cutoff value of diversity for differentiating landmark clusters from non-landmark clusters. The diversity of a landmark cluster must be below this cutoff.
#' @param size_cut the cutoff value of size i.e. number of cells for differentiating landmark clusters from non-landmark clusters. The number of cells in a landmark cluster must be greater than this cutoff.
#' @param saveRes a boolean to indicate whether to save result files
#' @return a dataframe of two columns, the first column is cell ID, the second column is the landmark cluster the cell belongs to
#' @importFrom stats hclust as.dist cor dist cutree kmeans
#' @importFrom graphics plot text abline legend 
#' @importFrom vegan diversity
#' @importFrom igraph compare
#' @export
#' @examples 
#' \dontrun{
#' rpkmFile = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
#' baseName = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG";
#' sampleFile = "sample_GSE60783_noOutlier.txt";
#' landmark_cluster <- landmark_designation(rpkmFile = rpkmFile,
#'                                          baseName = baseName,
#'                                          sampleFile = sampleFile,
#'                                          method = "diversity_size",
#'                                          numcluster = 11, diversity_cut=0.6,
#'                                          size_cut=0.05)
#' }
landmark_designation <- function(rpkmFile, baseName, sampleFile, distMethod="euclidean",
                                 method="kmeans", numcluster=NULL, diversity_cut=0.6,
                                 size_cut=0.05, saveRes=TRUE){
  if (is.character(rpkmFile)) {
    rpkm <- read.table(rpkmFile,sep="\t",header=T,row.names=1)
  } else {
    rpkm <- rpkmFile
  }
  log2rpkm <- apply(rpkm,c(1,2),function(x) if(x>1) log2(x) else 0)
  
  if (is.character(sampleFile)) {
    sample <- read.table(sampleFile,sep="\t",header=T)
  } else {
    sample <- sampleFile
  }
  row.names(sample) <- sample[,1]
  sample <- sample[colnames(log2rpkm),]
  
  if(distMethod=="euclidean"){
    hc <- stats::hclust(stats::dist(t(log2rpkm)),method="ward.D")
  } else {
    data.cor <- stats::cor(log2rpkm,method=distMethod);
    data.cor <- stats::as.dist(1-data.cor); #since the algorithm wants a distance measure, the resulting correlation is trasformed in a distance
    hc <- stats::hclust(data.cor,method="ward.D")
  }
  
  nmi_res <- data.frame(cluster_num=c(1:ncol(log2rpkm)),nmi=vector(length=ncol(log2rpkm)))
  for(k in 1:ncol(log2rpkm)){
    ct <- stats::cutree(hc,k=k)
    
    if(sum(names(ct)!=row.names(sample))==0){
      nmi_res[k,"nmi"] <- igraph::compare(as.factor(ct), as.factor(sample[,2]),method="nmi")
    }else{
      print("Error: the order of cells in sample.txt is different from rpkm file.\n")
    }
  }
  
  nmi_res_sort <- nmi_res[order(nmi_res$nmi,decreasing=TRUE),]  
  optimal_cluster_num <- nmi_res_sort[1,"cluster_num"]
  optimal_nmi <- nmi_res_sort[1,"nmi"]
  
  #pdf(paste(baseName,"_nmi.pdf",sep=""),height=8,width=8)
  #plot(nmi_res,xlab="Number of clusters",ylab="Normalized mutual information (NMI)")
  #points(nmi_res_sort[1,],col="blue",pch=1)
  #text(optimal_cluster_num+50,optimal_nmi+0.01,labels=paste("number of clusters = ",optimal_cluster_num,", nmi=",round(optimal_nmi,digit=2),sep=""),col="blue")
  ##dev.off()
  
  if(is.null(numcluster)){
    ct <- cutree(hc,k=optimal_cluster_num)
  }else{
    ct <- cutree(hc,k=numcluster)
  }
  
  clusters <- data.frame(table(ct))  
  clusters$diversity <- vegan::diversity(t(table(sample[,2],ct)), index = "shannon")
  if(method=="kmeans"){
    km <- kmeans(clusters[,2:3],2)
    km_c <- km$cluster
    km_center <- km$center
    if(km_center[1,1]>km_center[2,1] & km_center[1,2]>km_center[2,2]){
      km_c[km_c==1] <- "Landmark clusters"
      km_c[km_c==2] <- "Nonlandmark clusters"
    }else{
      km_c[km_c==2] <- "Landmark clusters"
      km_c[km_c==1] <- "Nonlandmark clusters"
    }
  }else if(method=="diversity"){
    km_c <- vector(length=nrow(clusters))
    km_c[clusters$diversity<=diversity_cut] <- "Landmark clusters"
    km_c[clusters$diversity>diversity_cut] <- "Nonlandmark clusters"
  }else if(method=="size"){
    km_c <- vector(length=nrow(clusters))
    km_c[clusters$Freq>=ceiling(size_cut*ncol(rpkm))+1] <- "Landmark clusters"
    km_c[clusters$Freq<ceiling(size_cut*ncol(rpkm))+1] <- "Nonlandmark clusters"
  }else if(method=="diversity_size"){
    km_c <- vector(length=nrow(clusters))
    km_c[clusters$diversity<=diversity_cut&clusters$Freq>=ceiling(size_cut*ncol(rpkm))+1] <- "Landmark clusters"
    km_c[clusters$diversity>diversity_cut|clusters$Freq<ceiling(size_cut*ncol(rpkm))+1] <- "Nonlandmark clusters"
  }
  
  cc_num <- as.character(clusters[km_c=="Landmark clusters","ct"])
  ct_cc <- ct[ct %in% cc_num]
  cc <- data.frame(cell=names(ct_cc),cluster=ct_cc)
  cc$true_group <- sample[as.character(cc$cell),"GroupID"]
  count <- table(cc$cluster,cc$true_group)
  num_name <- data.frame(num=vector(length=length(cc_num)),name=vector(length=length(cc_num)))
  
  if(nrow(count) == 0) {
    warning("No landmarks selected!! (added by Wouter as extra check)")
  }
  
  for(i in seq_len(nrow(count))){
    num_name$num[i] <- rownames(count)[i]
    num_name$name[i] <- colnames(count)[order(count[i,],decreasing=TRUE)[1]]
  }
  row.names(num_name) <- num_name[,"num"]
  cc$name <- num_name[as.character(cc$cluster),"name"]
  cc$landmark_cluster <- paste(cc$name,cc$cluster,sep="_")
  
  if(saveRes){
    write.table(cc[,c("cell","landmark_cluster")],paste(baseName,"_landmark_cluster.txt",sep=""),sep="\t",row.names=F)
  }
  
  
  cc_all <- data.frame(cell=names(ct),cluster=ct)
  cc_all$true_group <- sample[as.character(cc_all$cell),"GroupID"]
  count_all <- table(cc_all$cluster,cc_all$true_group)
  num_name <- data.frame(num=vector(length=nrow(count_all)),name=vector(length=nrow(count_all)))
  for(i in 1:nrow(count_all)){
    num_name$num[i] <- rownames(count_all)[i]
    num_name$name[i] <- colnames(count_all)[order(count_all[i,],decreasing=TRUE)[1]]
  }
  num_name$name <- paste(num_name$name,"_",num_name$num,sep="")
  row.names(num_name) <- num_name$num
  
  clusters$name <- num_name[as.character(clusters$ct),"name"]
  col_palatte <- c("red","black",height=2,width=2)
  if(saveRes){
    #pdf(paste(baseName,"_landmark_cluster.pdf",sep=""))
    plot(clusters$Freq,clusters$diversity,xlab="Cluster size",ylab="Cluster diversity",col= col_palatte[as.factor(km_c)],pch=16,xlim=c(min(clusters$Freq)-1,max(clusters$Freq)+1),ylim=c(min(clusters$diversity)-0.1,max(clusters$diversity)+0.1))
    text(clusters$Freq,clusters$diversity-0.02, labels=clusters$name, cex= 1, adj = c(0.5,0.9))
    if(method=="diversity"){
      abline(h=diversity_cut,lty=2)
    }else if(method=="size"){
      abline(v=ceiling(size_cut*ncol(rpkm))+1,lty=2)
    }else if(method=="diversity_size"){
      abline(h=diversity_cut,lty=2)
      abline(v=ceiling(size_cut*ncol(rpkm))+1,lty=2)
    }
    legend("topleft",legend=unique(km_c),pch=16,col= col_palatte[as.factor(unique(km_c))],bty="n")  
    ##dev.off()
  }
  #return(data.frame(nlm=sum(km_c=="Landmark clusters"),median_size=mean(clusters[km_c=="Landmark clusters","Freq"]),median_diversity=mean(clusters[km_c=="Landmark clusters","diversity"])))
  return(cc[,c("cell","landmark_cluster")])
}

#' build_network constructs weighted neighborhood network
#' @param exprs a data frame or matrix of expression data(ie. rpkm, TPM, fpkm) containing cells in columns and genes in rows
#' @param landmark_cluster a data frame or matrix of two columns or a tab delimited file of landmark cluster assignment of individual cells. The first column indicates cell ID, the second column indicates the landmark cluster which the cell was assigned to.  
#' @param distMethod the method for calculating dissimilarity between cells. distMethod can be one of "pearson", "kendall", "spearman" or "euclidean". Default is "euclidean".
#' @param baseName output directory
#' @param writeRes a boolean to indicate whether to save result files
#' @return a matrix of weighted neighborhood network, column and row names are landmarks, the values represent the weights of the edges connecting two landmarks
#' @export
#' @examples
#' \dontrun{
#' exprs = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
#' baseName = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG";
#' landmark = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG_landmark_cluster.txt"; 
#' # or landmark can be the return value of landmark_designation function 
#' neighbor_network <- build_network(exprs = exprs,
#'                                   landmark_cluster = landmark,
#'                                   baseName = baseName)
#' }
build_network <- function(exprs, landmark_cluster, distMethod = "euclidean", baseName = NULL, writeRes = TRUE){
  if(is.null(baseName)){
    if(is.character(exprs)){
       baseName <- sub(".txt","",exprs)
    }
  }
  
  if(is.character(exprs)){
    if(file.exists(exprs))
      exprs <- read.table(exprs, sep = "\t", header = T, row.names = 1, check.names = FALSE)
  }  
  log2exprs <- apply(exprs,c(1,2),function(x) if(x>1) log2(x) else 0)
  
  lm <- generate_lm(landmark_cluster=landmark_cluster,log2exprs)
  neighbor_network <- transition(log2exprs, lm, writeRes=writeRes, ifPlot=TRUE, textSize=30, baseName=baseName, distMethod = distMethod)
  return(neighbor_network) 
}

generate_lm <- function(landmark_cluster = "canonical_cluster.txt",log2exprs){
  if(is.character(landmark_cluster)){
    if(file.exists(landmark_cluster)){
       landmark_cluster <- read.table(landmark_cluster,sep="\t",header=T)
    }
  }
  
  log2exprs_landmark_cluster <- log2exprs[,as.character(landmark_cluster[,1])]
  lm <- aggregate(x = t(log2exprs_landmark_cluster), by = list(landmark_cluster[,2]), FUN = "mean")
  row.names(lm) <- lm[,1]
  lm <- lm[-1]
}

####### infer state transition between land marks
# distMethod can be one of "pearson", "kendall", "spearman" or "euclidean"

#' @importFrom igraph graph.adjacency E plot.igraph
transition <- function(log2exprs, lm, writeRes = TRUE, ifPlot = TRUE, textSize = 30, baseName, distMethod = "euclidean"){
  nb <- apply(log2exprs,2,function(x) neighbor(x,lm,distMethod=distMethod))
  if (writeRes) {
    write.table(t(nb),paste(baseName,"_neighbors_in_order.txt",sep=""),sep="\t")
  }
  row.names(nb) <- paste("neighbor",1:nrow(nb),sep="")
  nb12 <- table(nb["neighbor1",],nb["neighbor2",])
  
  absent <- rownames(lm)[!rownames(lm) %in% colnames(nb12)]
  if(length(absent)>0){
    temp <- matrix(0,nrow=nrow(nb12),ncol=length(absent))
    colnames(temp) <- absent
    rownames(temp) <- rownames(nb12)
    nb12 <- cbind(nb12,temp)
  }
  absent <- rownames(lm)[!rownames(lm) %in% rownames(nb12)]
  if(length(absent)>0){
    temp <- matrix(0,ncol=ncol(nb12),nrow=length(absent))
    rownames(temp) <- absent
    colnames(temp) <- colnames(nb12)
    nb12 <- rbind(nb12,temp)    
  }
  nb12 <- nb12[rownames(lm),rownames(lm)]
  
  if(writeRes){
    write.table(nb12,paste(baseName,"_NearestNeighbor_landMarks.txt",sep=""),sep="\t",col.names=NA)
    if(ifPlot){
      nn <- read.table(paste(baseName,"_NearestNeighbor_landMarks.txt",sep=""), row.names = 1, header = TRUE,check.names=F)
      network <- graph.adjacency(as.matrix(nn), weighted = TRUE)
      #pdf(paste(baseName,"_state_transition.pdf",sep=""),width=10,height=10)
      set.seed(1234)
      plot.igraph(network, edge.label=round(E(network)$weight),vertex.size=textSize,vertex.color="lightgrey",vertex.label.color="black",vertex.label.cex=1.3,edge.label.color="black",edge.label.cex=1.3)
      ##dev.off()
    }
  }  
  return(nb12)
}

#' trim_net trimms the weighted neighborhood network by removing edges of lower weights
#' @param nb12 a matrix of weighted neighborhood network, column and row names are landmarks, the values represent the weights of the edges connecting two landmarks
#' @param baseName a character string indicating the prefix name of resulting files
#' @param method trimming method, method can be one of "TrimNet" or "mst". When method="TrimNet" the initial node needs to be specified. Default is "mst"
#' @param start starting landmark, needs to be specified when method="TrimNet".
#' @param textSize the size of text
#' @param writeRes a boolean to indicate whether to save result files
#' @return a matrix of trimmed state transition network, column and row names are landmarks, the values are 0 or 1 indicating whether the two landmarks are connected.
#' @importFrom igraph graph.adjacency plot.igraph
#' @importFrom ape mst
#' @export
#' @examples
#' \dontrun{
#' baseName = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG";
#' trimmed_net <- trim_net(neighbor_network,textSize=30,
#'                         baseName = baseName,
#'                         method = "mst")
#' }
trim_net <- function(nb12,textSize=20,baseName=NULL,method="mst",start="MDP_6", writeRes = TRUE){
  if(method=="mst"){
    bb <- mst(-nb12)
    
    if (writeRes) {
      network <- graph.adjacency(as.matrix(bb))
      
      #pdf(paste(baseName,"_state_transition_mst.pdf",sep=""))
      plot.igraph(network, vertex.size=textSize, vertex.color="lightgrey")
      ##dev.off()  
      
      write.table(bb,paste(baseName,"_state_transition_mst.txt",sep=""),sep="\t",col.names=NA)
    }
    return(bb)
    
  }else if(method=="TrimNet"){
    bb <- nb12-nb12
    total <- unique(c(colnames(nb12),rownames(nb12)))
    if(is.null(start)){
      max_value <- which(nb12 == max(nb12), arr.ind = TRUE)
      added <- c(rownames(nb12)[max_value[,"row"]],colnames(nb12)[max_value[,"col"]])
    }else{
      added <- c(start)
    }
    not_added <- total[which(!total%in%added)]
    
    if(length(added)>0){
      if(length(added)==2){
        bb[added[1],added[2]] <- nb12[added[1],added[2]]
        bb[added[2],added[1]] <- nb12[added[2],added[1]]
      }else{
        for(m in 1:length(added)){
          for(n in 1:length(added)){
            bb[added[m],added[n]] <- nb12[added[m],added[n]]
          }
        }
      }
    }
    while (length(not_added)>0) {
      max <- 0
      max_i <- ""
      max_j <- ""
      for(i in 1:length(added)){
        for(j in 1:length(not_added)){
          lm_i <- added[i]
          lm_j <- not_added[j]
          sum <- 0
          if(lm_i %in% rownames(nb12) & lm_j %in% colnames(nb12)){
            sum <- sum + nb12[lm_i,lm_j]
          }
          if(lm_i %in% colnames(nb12) & lm_j %in% rownames(nb12)){
            sum <- sum + nb12[lm_j,lm_i]
          }
          if(sum>max){
            max <- sum
            max_i <- lm_i
            max_j <- lm_j
          }
        }        
      }
      added <- c(added,max_j)
      not_added <- total[which(!total%in%added)]
      bb[max_i,max_j] <- nb12[max_i,max_j]
      bb[max_j,max_i] <- nb12[max_j,max_i]
    }  
    network <- graph.adjacency(as.matrix(bb))
    
    #pdf(paste(baseName,"_state_transition_backbone_draf_TrimNet_fat.pdf",sep=""))
    #plot(network,vertex.size=textSize)
    ###dev.off()  
    #write.table(bb,paste(baseName,"_state_transition_backbone_draf_TrimNet_fat.txt",sep=""),sep="\t",col.names=NA)
    
    bb_thin <- apply(bb>0,c(1,2),sum)
    
    if (writeRes) {
      network <- graph.adjacency(as.matrix(bb_thin))
      #pdf(paste(baseName,"_state_transition_TrimNet.pdf",sep=""))
      #set.seed(3959)
      
      set.seed(3);    
      plot.igraph(network, vertex.size=textSize, vertex.color="lightgrey", vertex.label.color="black", vertex.label.cex=0.8,edge.arrow.size = 0.5)
      ##dev.off()
      write.table(bb_thin,paste(baseName,"_state_transition_TrimNet.txt",sep=""),sep="\t",col.names=NA)
    }
    
    return(bb_thin)
  } 
}

###### Find neighboring landmarks ######
# distMethod can be one of "pearson", "kendall", "spearman" or "euclidean"

neighbor <- function(log2rpkm_percell,lm,distMethod="euclidean"){
  if(sum(names(log2rpkm_percell)!=colnames(lm))==0){
    data4dist <- data.frame(t(lm),cell=log2rpkm_percell,check.names=FALSE)
  }
  ####
  if(distMethod=="spearman" | distMethod=="pearson"){
    data.cor<-cor(data4dist,method=distMethod);
    d <- as.matrix(as.dist(1-data.cor)); 
  }else{
    d <- as.matrix(dist(t(data4dist),distMethod))
  }
  d_cell <- d["cell",]
  ord <- order(d_cell)
  names(d_cell)[ord[2:length(ord)]]
}

#### place cells in linear order along backbone paths
#### if method=1, the cells are ordered based on coordinates calcualted from distance to landmarks
#### if method=2, the cells are ordered based on the distance to the nearest neighbor landmark

nbor <- function(log2rpkm,lm,lm_order=c("T0_1","T24_8","T48_10","T72_13"),if_bb_only=FALSE,dist_method="euclidean",method=1){
  nb <- apply(log2rpkm,2,function(x) neighbor(x,lm))
  row.names(nb) <- paste("neighbor",1:nrow(nb),sep="")
  
  if(nrow(nb)<ncol(nb)){
    nb <- t(nb)
  }
  
  k <- 1
  res <- list()
  for(i in 1:length(lm_order)){
    for(j in 1:length(lm_order)){
      nb1 <- lm[lm_order[i],]
      nb2 <- lm[lm_order[j],]
      
      if(!if_bb_only){
        if(method==1){
          if(i<j){
            d <- dist(rbind(nb1,nb2),dist_method)
          }else{
            d <- -dist(rbind(nb1,nb2),dist_method)
          }
          
          cell <- rownames(nb[nb[,"neighbor1"]==row.names(nb1) & nb[,"neighbor2"]==row.names(nb2),]) 
          cell_x <- apply(log2rpkm[,cell],2,function(x) transition_cor(nb1,nb2,x,d)) 
          
          res[[k]] = list(nb1=row.names(nb1),nb2=row.names(nb2),cell_x=sort(cell_x),order=names(cell_x[order(cell_x)]))
        }else{
          if(i!=j){
            cell_close2nb1 <- row.names(nb[nb[,"neighbor1"]==row.names(nb1) & nb[,"neighbor2"]==row.names(nb2),])
            log2rpkm_close2nb1 <- log2rpkm[,cell_close2nb1]
            dist_close2nb1 <- apply(log2rpkm_close2nb1,2,function(x) dist(rbind(nb1,x)))
            if(i<j){
              res[[k]] = list(nb1=row.names(nb1),nb2=row.names(nb2),cell_x=sort(dist_close2nb1),order=names(sort(dist_close2nb1)))                
            }else{
              res[[k]] = list(nb1=row.names(nb1),nb2=row.names(nb2),cell_x=sort(dist_close2nb1,decreasing=TRUE),order=names(sort(dist_close2nb1,decreasing=TRUE)))                
            }        
          }
        }
        k <- k+1
      }else{
        if(j==i+1 | i==j+1){
          if(method==1){
            if(i<j){
              d <- dist(rbind(nb1,nb2),dist_method)
            }else{
              d <- -dist(rbind(nb1,nb2),dist_method)
            }
            
            cell <- rownames(nb[nb[,"neighbor1"]==row.names(nb1) & nb[,"neighbor2"]==row.names(nb2),]) 
            cell_x <- apply(log2rpkm[,cell],2,function(x) transition_cor(nb1,nb2,x,d)) 
            
            res[[k]] = list(nb1=row.names(nb1),nb2=row.names(nb2),cell_x=sort(cell_x),order=names(cell_x[order(cell_x)]))
          }else{
            cell_close2nb1 <- row.names(nb[nb[,"neighbor1"]==row.names(nb1) & nb[,"neighbor2"]==row.names(nb2),])
            log2rpkm_close2nb1 <- log2rpkm[,cell_close2nb1]
            dist_close2nb1 <- apply(log2rpkm_close2nb1,2,function(x) dist(rbind(nb1,x)))
            if(i<j){
              res[[k]] = list(nb1=row.names(nb1),nb2=row.names(nb2),cell_x=sort(dist_close2nb1),order=names(sort(dist_close2nb1)))                
            }else{
              res[[k]] = list(nb1=row.names(nb1),nb2=row.names(nb2),cell_x=sort(dist_close2nb1,decreasing=TRUE),order=names(sort(dist_close2nb1,decreasing=TRUE)))                
            }
          }
          k <- k+1
        }
      }
    }
  }
  return(res)
}


transition_cor <- function(nb1,nb2,cell,nb2_x,dist_method="euclidean"){
  #d <- dist(rbind(lm1,lm2))
  
  d1 <- dist(rbind(nb1,cell),dist_method)
  d2 <- dist(rbind(nb2,cell),dist_method)
  
  x <- (d1*d1-d2*d2+nb2_x*nb2_x)/(2*nb2_x)
  return(x)
}

#' color_code_node_2 plot state transition network in which nodes i.e. landmarks are color-coded by average expression of the given gene
#' @param networkFile a tab delimited file containing a matrix of trimmed state transition network, column and row names are landmarks, the values are 0 or 1 indicating whether the two landmarks are connected.
#' @param rpkmFile a tab delimited txt file of expression data, containing cells in columns and genes in rows 
#' @param lmFile a tab delimited file of landmark cluster assignment of individual cells. The first column indicates cell ID, the second column indicates the landmark cluster which the cell was assigned to.  
#' @param geneName gene name or a vector of gene names
#' @param baseName prefix name of resulting files
#' @param seed the seed
#' @import plotrix
#' @importFrom stats aggregate
#' @export
#' @examples
#' \dontrun{
#' rpkmFile="TPM_GSE60783_noOutlier.txt";
#' lmFile="TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG_landmark_cluster.txt";
#' network="TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG_state_transition_mst.txt";
#' # network can be the return value of trim_net 
#' color_code_node_2(networkFile=network,
#'                   rpkmFile=rpkmFile,
#'                   lmFile=lmFile,
#'                   geneName=c("Irf8","Id2","Batf3"),
#'                   baseName="cDC1_marker",
#'                   seed=NULL)
#' }
color_code_node_2 <- function(networkFile, rpkmFile,lmFile, geneName,baseName = NULL, seed=NULL){
  if(is.character(networkFile)){
    if(file.exists(networkFile)){
       network_matrix <- read.table(networkFile,sep="\t",header=T,row.names=1,check.names=F)
    }
  }else{
    network_matrix <- networkFile
  }
  
  lm <- read.table(lmFile,sep="\t",header=T,row.names=1)  
  rpkm <- read.table(rpkmFile,sep="\t",header=T,row.names=1)
  
  if(length(geneName)<2){
    rpkm_g <- rpkm[geneName,]
    rpkm_g_order <- rpkm_g[row.names(lm)]
    
    rpkm_g_median <- aggregate(x=t(rpkm_g_order),by=list(lm[,1]),FUN="median")
    if(max(rpkm_g_median[,2])==0){
      rpkm_g_median <- aggregate(x=t(rpkm_g_order),by=list(lm[,1]),FUN="mean")
    }
    row.names(rpkm_g_median) <- rpkm_g_median[,1]
    rpkm_g_median <- rpkm_g_median[colnames(network_matrix),geneName]
    
    #pdf(paste(baseName,"_",geneName,".pdf",sep=""))
    if(is.null(seed)){
      set.seed(3959)
    }
    igraph_annot(adjMatrix=network_matrix,attrVec=rpkm_g_median,main=geneName)
    ##dev.off()
  }else{     
    for(i in 1:length(geneName)){
      rpkm_g <- rpkm[geneName[i],]
      rpkm_g_order <- rpkm_g[row.names(lm)]
      
      rpkm_g_median <- aggregate(x=t(rpkm_g_order),by=list(lm[,1]),FUN="median")
      if(max(rpkm_g_median[,2])==0){
        rpkm_g_median <- aggregate(x=t(rpkm_g_order),by=list(lm[,1]),FUN="mean")
      }
      row.names(rpkm_g_median) <- rpkm_g_median[,1]
      rpkm_g_median <- rpkm_g_median[colnames(network_matrix),geneName[i]]
      
      #pdf(paste(baseName,"_",geneName[i],".pdf",sep=""))
      if(is.null(seed)){
        set.seed(3959)
      }
      igraph_annot(adjMatrix=network_matrix,attrVec=rpkm_g_median,main=geneName[i])
      #dev.off()
    }       
  }
}

#' @importFrom igraph graph.adjacency plot.igraph layout.kamada.kawai
#' @importFrom grDevices colorRampPalette palette
#' @importFrom stats quantile
#' @importFrom graphics title image
igraph_annot <- function(adjMatrix, attrVec, main = "", mode = "directed", xlab="Expression Value",
                         palette = "bluered", layout = layout.kamada.kawai, pctile_color=c(0.02,0.98) ){  
  graph <- graph.adjacency(as.matrix(adjMatrix), mode=mode, weighted = TRUE, diag=FALSE)
  graph_l <- layout
  
  if (palette == "jet")
    palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  else if (palette == "bluered")
    palette <- colorRampPalette(c( "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
  else
    stop("Please use a supported color palette.  Options are 'bluered' or 'jet'") 
  if (!is.vector(pctile_color) || length(pctile_color) != 2) 
    stop("pctile_color must be a two element vector with values in [0,1]")
  colorscale <- palette(100)
  boundary <- quantile(attrVec, probs=pctile_color, na.rm=TRUE)
  boundary <- round(boundary, 2)
  grad <- seq(boundary[1], boundary[2], length.out=length(colorscale))
  color <- colorscale[findInterval(attrVec, grad, all.inside=TRUE)]
  graph$color <- color
  
  #pdf(paste(out_dir,basename(f),".",name,".pdf",sep=""))
  plot.igraph(graph, layout = graph_l,
       vertex.color= graph$color,
       vertex.size = 25,
       vertex.shape="circle",
       vertex.frame.color= color,
       vertex.label.cex = 0.8,
       vertex.label.color = "black",
       vertex.label.family = "sans",
       #edge.width=graph$weight, 
       #edge.color="grey",
       #edge.label=round(E(graph)$weight),
       #edge.label.cex = 0.6,
       edge.arrow.size = 0.5,
       edge.arrow.width=1,
       main=main
  )
  
  title(main=main)
  subplot(
    image(
      grad, c(1), matrix(1:length(colorscale),ncol=1), col=colorscale,
      xlab=xlab,
      ylab="", yaxt="n", xaxp=c(boundary,1)
    ),
    x="right,bottom",size=c(1,.20)
  )
  ##dev.off()
}



# ## 2D interactive plot
# tkplot(mst_graph, vertex.size=3,vertex.label.cex = 0.6)
# 
# ## 3D plot
# coords <- layout.fruchterman.reingold(mst_graph, dim=3)
# rglplot(mst_graph, layout=coords,vertex.size=5, vertex.label=NA)
# #rglplot(mst_graph, layout=coords,vertex.size=2, vertex.label=NA, vertex.color=V(mst_graph)$exprs)

# The following functions are copied from the TeachingDemos package, 
# authored by Greg Snow <greg.snow at imail.org>, and licensed under
# the Artistic-2.0 license.

#' @importFrom graphics par
#' @importFrom grDevices xy.coords
cnvrt.coords <- function(x,y=NULL,input=c('usr','plt','fig','dev','tdev')) {
  
  input <- match.arg(input)
  xy <- xy.coords(x,y, recycle=TRUE)
  
  cusr <- par('usr')
  cplt <- par('plt')
  cfig <- par('fig')
  cdin <- par('din')
  comi <- par('omi')
  cdev <- c(comi[2]/cdin[1],(cdin[1]-comi[4])/cdin[1],
            comi[1]/cdin[2],(cdin[2]-comi[3])/cdin[2])
  
  if(input=='usr'){
    usr <- xy
    
    plt <- list()
    plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
    plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])
    
    fig <- list()
    fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
    fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
    
    dev <- list()
    dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
    dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]
    
    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]
    
    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }
  
  if(input=='plt') {
    
    plt <- xy
    
    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]
    
    fig <- list()
    fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
    fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
    
    dev <- list()
    dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
    dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]
    
    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]
    
    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }
  
  if(input=='fig') {
    
    fig <- xy
    
    plt <- list()
    plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
    plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])
    
    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]
    
    dev <- list()
    dev$x <- fig$x*(cfig[2]-cfig[1])+cfig[1]
    dev$y <- fig$y*(cfig[4]-cfig[3])+cfig[3]
    
    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]
    
    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }
  
  if(input=='dev'){
    dev <- xy
    
    fig <- list()
    fig$x <- (dev$x-cfig[1])/(cfig[2]-cfig[1])
    fig$y <- (dev$y-cfig[3])/(cfig[4]-cfig[3])
    
    plt <- list()
    plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
    plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])
    
    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]
    
    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]
    
    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }
  
  if(input=='tdev'){
    tdev <- xy
    
    dev <- list()
    dev$x <- (tdev$x-cdev[1])/(cdev[2]-cdev[1])
    dev$y <- (tdev$y-cdev[3])/(cdev[4]-cdev[3])
    
    fig <- list()
    fig$x <- (dev$x-cfig[1])/(cfig[2]-cfig[1])
    fig$y <- (dev$y-cfig[3])/(cfig[4]-cfig[3])
    
    plt <- list()
    plt$x <- (fig$x-cplt[1])/(cplt[2]-cplt[1])
    plt$y <- (fig$y-cplt[3])/(cplt[4]-cplt[3])
    
    usr <- list()
    usr$x <- plt$x*(cusr[2]-cusr[1])+cusr[1]
    usr$y <- plt$y*(cusr[4]-cusr[3])+cusr[3]
    
    tdev <- list()
    tdev$x <- dev$x*(cdev[2]-cdev[1])+cdev[1]
    tdev$y <- dev$y*(cdev[4]-cdev[3])+cdev[3]
    
    return( list( usr=usr, plt=plt, fig=fig, dev=dev, tdev=tdev ) )
  }
  
}

#' @importFrom graphics par
#' @importFrom grDevices xy.coords
subplot <- function(fun, x, y=NULL, size=c(1,1), vadj=0.5, hadj=0.5,
                    inset=c(0,0), type=c('plt','fig'), pars=NULL){
  
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  
  type <- match.arg(type)
  
  if(missing(x)) x <- locator(2)
  
  if(is.character(x)) {
    if(length(inset) == 1) inset <- rep(inset,2)
    x.char <- x
    tmp <- par('usr')
    x <- (tmp[1]+tmp[2])/2
    y <- (tmp[3]+tmp[4])/2
    
    if( length(grep('left',x.char, ignore.case=TRUE))) {
      x <- tmp[1] + inset[1]*(tmp[2]-tmp[1])
      if(missing(hadj)) hadj <- 0
    }
    if( length(grep('right',x.char, ignore.case=TRUE))) {
      x <- tmp[2] - inset[1]*(tmp[2]-tmp[1])
      if(missing(hadj)) hadj <- 1
    }
    if( length(grep('top',x.char, ignore.case=TRUE))) {
      y <- tmp[4] - inset[2]*(tmp[4]-tmp[3])
      if(missing(vadj)) vadj <- 1
    }
    if( length(grep('bottom',x.char, ignore.case=TRUE))) {
      y <- tmp[3] + inset[2]*(tmp[4]-tmp[3])
      if(missing(vadj)) vadj <- 0
    }
  }
  
  xy <- xy.coords(x,y)
  
  if(length(xy$x) != 2){
    pin <- par('pin')
    tmp <- cnvrt.coords(xy$x[1],xy$y[1],'usr')$plt
    
    x <- c( tmp$x - hadj*size[1]/pin[1],
            tmp$x + (1-hadj)*size[1]/pin[1] )
    y <- c( tmp$y - vadj*size[2]/pin[2],
            tmp$y + (1-vadj)*size[2]/pin[2] )
    
    xy <- cnvrt.coords(x,y,'plt')$fig
  } else {
    xy <- cnvrt.coords(xy,'usr')$fig
  }
  
  par(pars)
  if(type=='fig'){
    par(fig=c(xy$x,xy$y), new=TRUE)
  } else {
    par(plt=c(xy$x,xy$y), new=TRUE)
  }
  fun
  tmp.par <- par(no.readonly=TRUE)
  
  return(invisible(tmp.par))
}