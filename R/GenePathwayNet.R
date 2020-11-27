#  Generation gene-pathway bipartite matrix

#' Generate the bipartite matrix from gene-pathway interaction networks
#' @description The function \code{genepathwayAdj} is to generate the adjacency
#'   matrix from gene-pathway interaction network
#' @param drugAdj drug adjacency matrix (sparse matrix) generated from drug-drug
#'   association network via \code{AdjMatrix_drugs}
#' @param pathwayadj pathway-pathway association adjacency matrix (sparse
#'   matrix)
#' @param geneadj gene-gene association network adjacency matrix (sparse matrix)
#' @param gpweight the weight of gene-pathway edges
#' @param load_dir path to load or save modeling data files
#' @return \code{genepathwayAdj} gene-pathway adjacency matrix
#' @export
#'
#' @examples
#'

genepathwayAdj <- function(drugAdj,
                           pathwayadj,
                           geneadj,
                           gpweight,
                           load_dir){

  # load(system.file("data", "KEGG_pw_gene.RData", package = "DComboNet"))

  genepathwayAdj = matrix(0,nrow(pathwayadj),nrow(geneadj))
  rownames(genepathwayAdj) = rownames(pathwayadj)
  colnames(genepathwayAdj) = colnames(geneadj)
  pathways <- rownames(pathwayadj)
  N = nrow(pathwayadj)

  # genepathway = read.table(paste0(load_dir,"pathway/KEGG_ID_gene.txt"),sep="\t",header = T,stringsAsFactors = F)
  #  names(genepathway) = c("KEGG_ID","pathway_name","Gene")
  names(genepathway) = c("Gene","KEGG_ID")
  for(i in 1:N){
    pathway <- pathways[i]
    gene <- genepathway[which(genepathway$KEGG_ID == pathway),]$Gene

    if(length(gene) > 0){
      for(j in 1:length(gene)){
        if(!is.na(gene[j])){
          # identify the drug position on the matrix.
          index_gene <- which(colnames(genepathwayAdj) %in% gene[j])
          if(length(index_gene) == 1){
            genepathwayAdj[i,index_gene] <- gpweight
          }
        }
      }
    }
  }
  return(genepathwayAdj)
}
