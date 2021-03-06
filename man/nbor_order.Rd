% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mpath.R
\name{nbor_order}
\alias{nbor_order}
\title{nbor_order sorts individual cells according to their various stages during transition to resemble the landmark-to-landmark continuum}
\usage{
nbor_order(exprs, ccFile, lm_order = c("CD115+CDP_1", "PreDC_9", "PreDC_10",
  "PreDC_11"), if_bb_only = FALSE, method = 1, writeRes = TRUE)
}
\arguments{
\item{exprs}{a data frame or matrix of expression data(ie. rpkm, TPM, fpkm) or a tab delimited txt file of expression data, containing cells in columns and genes in rows}

\item{ccFile}{a data frame or matrix of two columns or a tab delimited file of landmark cluster assignment of individual cells. The first column indicates cell ID, the second column indicates the landmark cluster which the cell was assigned to.}

\item{lm_order}{a vector of landmark IDs indicating along which path the cells are to be sorted}

\item{if_bb_only}{a boolean to indicate if only cells on backbone will be sorted. Default is FALSE}

\item{method}{1 or 2 to indicate which method to be used for sorting. Default is 1}

\item{writeRes}{a boolean to indicate whether to save result files}
}
\value{
a vector of  re-orderd cell IDs
}
\description{
nbor_order sorts individual cells according to their various stages during transition to resemble the landmark-to-landmark continuum
}
\examples{
\dontrun{
exprs = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG.txt";
ccFile = "TPM_GSE60783_noOutlier_geneQC0.05anyGroup_CD4vsCD8DEG_landmark_cluster.txt";
order <- nbor_order(exprs = exprs,
                    ccFile = ccFile, 
                    lm_order = c("CDP_2","CDP_1","preDC_9","preDC_3"),
                    if_bb_only=FALSE,
                    method=1)
}
}
