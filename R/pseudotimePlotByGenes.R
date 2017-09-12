#' pseudotimeplotByGenes
#' 
#' pseudotimeplotByGenes
#' 
#' @param exprs a data frame or matrix of log transformed expression data, with row of cells and column of genes
#' @param if_log2 whether exprs is log2 transformed
#' @param cell_annotation a two column data frame of matrix annotating cells(cell ID or name) with cell types
#' @param cell_order a vector stroing the order of cells with cell ID or name, same as appeared in row names of \code{exprs}
#' @param plot_genes a vector storing the genes selected for plot, same as appeared in the column names of \code{exprs}, if NULL, all genes in exprs will be selected
#' @param reverse_order reverse the order of the pseudotime
#' @param min_expr the threshold for cutting of the cell expressions in regression values, values lower than this will be forced to \code{min_exprs}
#' @param cell_size the size of cells in the plot
#' @param plot_cols plot colours
#' @param trend_formula the formula for regression analysis
#' @return a object of ggplot
#' @importFrom reshape2 melt
#' @import ggplot2
#' @importFrom VGAM vgam
#' @importFrom plyr ddply
#' @export
#' @examples
#' \dontrun{
#' pseudotimePlotByGenes(exprs = "FULL.log2TPM.txt", 
#'                       cell_annotation = "splAnnotation_outlierRemoved.txt", 
#'                       cell_order = "uspin.PCA.0.03.seed1.txt", 
#'                       plot_genes = "genes.txt", 
#'                       reverse_order = TRUE, plot_cols = 3)
#' }
pseudotimePlotByGenes <-function(exprs, if_log2=TRUE, cell_annotation, cell_order, plot_genes=NULL, 
                                 reverse_order = FALSE, min_expr=-3, cell_size=2, plot_cols = NULL, 
                                 trend_formula="expression ~ sm.ns(Pseudotime, df=3)"){
    
    ## read data if file names provided
    if(is.character(exprs)){
        if(file.exists(exprs))
            exprs <- read.table(exprs, header = TRUE, row.names = 1, check.names = FALSE)
    }
    if(if_log2){
        exprs <- apply(exprs,c(1,2),function(x) if(x>1) log2(x) else 0)
      
    }
        
    if(is.character(cell_annotation)){
        if(file.exists(cell_annotation))
            cell_annotation <- read.table(cell_annotation, header = TRUE)
    }
    
    if(is.character(cell_order)){
        if(file.exists(cell_order))
            cell_order <- as.character(read.table(cell_order)[,1])
    }
    
    if(is.character(plot_genes)){
        if(file.exists(plot_genes))
            plot_genes <- as.character(read.table(plot_genes)[,1])
    }
    
    ## checking the data type
    if(all(cell_order %in% colnames(exprs))){
        exprs <- data.frame(t(exprs),check.names = FALSE)
    }else if (any(!(cell_order %in% row.names(exprs)))){
        stop("cell IDs in cell_order doesn's matches in exprs!")
    }
    
    if(is.null(plot_genes)){
        plot_genes <- colnames(exprs)
    }else if(any(!(plot_genes %in% colnames(exprs)))){
        stop("plot_genes doesn't match in exprs!")
    }
        
    if(any(cell_annotation[ ,1] %in% row.names(exprs))){
        cell_anno_cellid <- cell_annotation[,1]
        cell_anno_celltype <- cell_annotation[,2]
    }else if(any(cell_annotation[ ,2] %in% row.names(exprs))){
        cell_anno_cellid <- cell_annotation[,2]
        cell_anno_celltype <- cell_annotation[,1]
    }else{
        stop("cell ID in cell_annotation doesn't match in exprs")
    }
    
    if(reverse_order)
        cell_order <- rev(cell_order)
    plotExprs <- exprs[row.names(exprs) %in% cell_order, plot_genes]
    plotExprs$Pseudotime <- match(row.names(plotExprs), cell_order)
    plotExprs$celltype <- cell_anno_celltype[match(row.names(plotExprs), cell_anno_cellid)]
    plotExprs[,which(apply(plotExprs,2,max)==0)] <- 0.000000000000000000001     
    cds_exprs <- melt(plotExprs, id.vars = c("Pseudotime", "celltype"), 
                      variable.name = "genes", value.name = "expression")
    cds_exprs$genes <- factor(cds_exprs$genes)
    
    merged_df_with_vgam <- ddply(cds_exprs, .(genes), function(x) { 
        fit_res <- tryCatch({
            vg <- suppressWarnings(vgam(formula = as.formula(trend_formula), 
                                        family = VGAM::tobit(Lower = min_expr, lmu = "identitylink"), 
                                        data = x, maxit=30, checkwz=FALSE))
            res <- predict(vg, type="response")
            res[res < min_expr] <- min_expr
            res
        }
        ,error = function(e) {
            print("Error!")
            print(e)
            res <- rep(NA, nrow(x))
            res
        }
        )
        expectation = fit_res
        data.frame(Pseudotime=x$Pseudotime, expectation=expectation)
    })
    
    ## gene profile according to pseudotime
    monocle_theme_opts <- function()
    {
        theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
            theme(panel.border = element_blank(), axis.line = element_line()) +
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
            theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
            theme(panel.background = element_rect(fill='white')) +
            theme(legend.position = "right") +
            theme(axis.title = element_text(size = 15))
    }
    color_by="celltype"
    if(is.null(plot_cols)){
        plot_cols <- round(sqrt(length(plot_genes)))
    }
    
    q <- ggplot(aes(Pseudotime, expression), data=cds_exprs) 
    q <- q + geom_point(aes_string(color=color_by), size=I(cell_size))
    q <- q + geom_line(aes(Pseudotime, expectation), data=merged_df_with_vgam)
    q <- q + facet_wrap(~genes, ncol=plot_cols, scales="free_y")
    q <- q + ylab("Expression (log2)") + xlab("Cell order") + theme_bw()
    q <- q + guides(colour = guide_legend(title = "Cell Type", override.aes = list(size = cell_size*3)))
    q <- q + monocle_theme_opts() 
    
    return(q)
    
}

