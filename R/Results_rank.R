#' Drug result rank
#' @description The function \code{rank_drugs} is to rank drugs according to rwr
#'   result
#' @param Num.Gene Number of nodes in gene-gene association network
#' @param Num.Drug Number of nodes in drug-drug association network
#' @param RWR.result data.frame, global proximity between drug seed and other
#'   drug in drug-gene/pathway interaction network computed via \code{rwr}
#'   function
#' @param Drug_seeds drug that take as seed in \code{seedscore} function
#' @return sorted global proximity between drug seed and other drug in
#'   drug-gene/pathway interaction network
#' @export



rank_drugs <- function(Num.Gene,
                       Num.Drug,
                       RWR.result,
                       Drug_seeds){

  drugs_rank <- data.frame(DrugID = character(length = Num.Drug), Score = 0)
  drugs_rank$DrugID <- row.names(RWR.result)[1:(Num.Drug)]
  drugs_rank$Score <- RWR.result[(1):(Num.Drug), 1]
  drugs_rank_sort <- drugs_rank[with(drugs_rank, order(-Score, DrugID)), ]
  drugs_rank_sort <- drugs_rank_sort[which(!drugs_rank_sort$DrugID %in% Drug_seeds),]

  return(drugs_rank_sort)
}

#' Gene result rank
#' @description The function \code{rank_genes} is to rank genes according to rwr result
#' @param Num.Drug Number of nodes in drug-drug association network
#' @param Num.Gene Number of nodes in gene-gene association network
#' @param RWR.result data.frame, global proximity between drug seed and other
#'   drug in drug-gene/pathway interaction network computed via \code{rwr}
#'   function
#' @param Gene_seeds cancer related gene seeds (optional)
#' @return global proximity between drug seed and other
#'   genes in drug-gene/pathway interaction network
#' @export
#'
#'
rank_genes <- function(Num.Drug,
                       Num.Gene,
                       RWR.result,
                       Gene_seeds = NULL){

  genes_rank <- data.frame(GeneNames = character(length = Num.Gene), Score = 0)
  genes_rank$GeneNames <- row.names(RWR.result)[(1+Num.Drug):(Num.Drug+Num.Gene)]
  genes_rank$Score <- RWR.result[(1+Num.Drug):(Num.Drug+Num.Gene)]
  genes_rank_sort <- genes_rank[with(genes_rank, order(-Score, GeneNames)), ]

  if(is.null(Gene_seeds)!=0){
    return(genes_rank_sort)
  }else{
    genes_rank_sort <- genes_rank_sort[which(!genes_rank_sort$GeneNames %in% Gene_seeds),]
    return(genes_rank_sort)

  }
}


#' Pathway result rank
#' @description The function \code{rank_pathways} is to rank pathways according
#'   to rwr result
#' @param Num.Gene Number of nodes in gene-gene association network
#' @param Num.Drug Number of nodes in drug-drug association network
#' @param RWR.result data.frame, global proximity between drug seed and other
#'   drug in drug-gene/pathway interaction network computed via \code{rwr}
#'   function
#' @return global proximity between drug seed and other pathways in
#'   drug-gene/pathway association network
#' @export
#'

rank_pathways <- function(Num.Gene,
                          Num.Drug,
                          RWR.result){

  pathways_rank <- data.frame(PathwayID = character(length = nrow(RWR.result) - (Num.Gene+Num.Drug)), Score = 0)
  pathways_rank$PathwayID <- row.names(RWR.result)[(Num.Gene+Num.Drug+1):nrow(RWR.result)]
  pathways_rank$Score <- RWR.result[(Num.Gene+Num.Drug+1):nrow(RWR.result), 1]
  pathways_rank_sort <- pathways_rank[with(pathways_rank, order(-Score, PathwayID)), ]

  return(pathways_rank_sort)
}


