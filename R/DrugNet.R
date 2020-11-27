#' Generate drug adjacency matrix
#'
#' @description The function \code{AdjMatrix_drugs} is to generate adjacency
#'   matrix for drug-drug interaction network weighted by the integrated
#'   pharmacological similarity values obtained from \code{DrugNetFeature}.
#' @param x drug-drug similarity feature file generated from
#'   \code{DrugNetFeature}
#' @param weighted \code{TRUE/FALSE}, choose if the edge of drug network should
#'   be weighted by feature values. The abbreviationS "\code{T/F}" are also
#'   acceptable
#' @param  weight.type if parameter "\code{weighted = TRUE/T}", choose the
#'   feature type to weight the edge of drug-drug network. The options are:
#'   "ATCsim","STsim","SEsim","integrated_score","integrated_score2"
#'
#' @import igraph
#' @return \code{AdjMatrix_drugs} Drug adjacency matrix generated from drug-drug
#'   association network
#' @export
#'
#' @examples
#'
#' # drugfeature = read.table(paste0(load_dir, "drug_net/features.csv", sep =
#' # ",", header = TRUE, stringsAsFactors = FALSE))
#' \dontrun{
#' AdjMatrix_drugs(x = drugfeature, weighted = TRUE, weight.type = "integrated_score")
#' }

AdjMatrix_drugs <- function(x,
                            weighted = c("TRUE","T","FALSE","F"),
                            weight.type = c("ATCsim","STsim","SEsim","integrated_score","integrated_score2")){

  if (requireNamespace("igraph")){

    if(weighted == FALSE & weighted == F){
      x$basic = 0
      x[x$TAG=="P",]$basic = 1
      Drug_Network <- igraph::graph.data.frame(x[c(1:2)],directed=FALSE)
      igraph::E(Drug_Network)$weight <- x$basic
      Drug_Network <- igraph::simplify(Drug_Network, remove.multiple = T, remove.loops = T)
      AdjMatrix_drugs <- igraph::as_adjacency_matrix(Drug_Network,sparse = TRUE)

    }else if(weighted == TRUE & weighted == T){

      Drug_Network <- igraph::graph.data.frame(x[c(1:2)],directed=FALSE)
      if(weight.type == "ATCsim"){
        igraph::E(Drug_Network)$weight <- x$ATCsim
        Drug_Network <- igraph::simplify(Drug_Network, remove.multiple = T, remove.loops = T)
        AdjMatrix_drugs = igraph::get.adjacency(Drug_Network, attr="weight")
      }else if(weight.type == "STsim"){
        igraph::E(Drug_Network)$weight <- x$STsim
        Drug_Network <- igraph::simplify(Drug_Network, remove.multiple = T, remove.loops = T)
        AdjMatrix_drugs = igraph::get.adjacency(Drug_Network, attr="weight")
      }else if(weight.type == "SEsim"){
        igraph::E(Drug_Network)$weight <- x$SEsim
        Drug_Network <- igraph::simplify(Drug_Network, remove.multiple = T, remove.loops = T)
        AdjMatrix_drugs = igraph::get.adjacency(Drug_Network, attr="weight")
      }else if(weight.type == "integrated_score"){
        igraph::E(Drug_Network)$weight <- x$integrated_score
        Drug_Network <- igraph::simplify(Drug_Network, remove.multiple = T, remove.loops = T)
        AdjMatrix_drugs = igraph::get.adjacency(Drug_Network, attr="weight")
      }else if(weight.type == "integrated_score2"){
        igraph::E(Drug_Network)$weight <- x$integrated_score2
        Drug_Network <- igraph::simplify(Drug_Network, remove.multiple = T, remove.loops = T)
        AdjMatrix_drugs = igraph::get.adjacency(Drug_Network, attr="weight")
      }
    }

    #rm(Drug_Network);gc()
    return(AdjMatrix_drugs)
  }
}
