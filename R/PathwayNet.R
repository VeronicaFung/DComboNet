#' Generate the adjacency matrix from pathway-pathway association network
#' @description The function \code{pathwayNet.Adj} is to generate adjacency
#'   matrix from pathway-pathway interaction network
#' @param PathwayNetwork \code{data.frame} of pathway-pathway associations
#' @import Matrix
#' @import igraph
#' @return \code{pathwayadj} pathway-pathway association adjacency matrix
#'   (sparse matrix)
#' @export
#'
#' @examples
#' # require(igraph)
#' p <- data.frame(V1 = c("p1", "p2", "p3", "p4"), V2 = c("p5", "p2", "p6",
#' "p1"))
#'
#' pathwayNet.Adj(p)
#'

pathwayNet.Adj <- function(PathwayNetwork){
  if(requireNamespace(c("Matrix","igraph"))){

    WWI.net <- igraph::graph.data.frame(PathwayNetwork)
    N = length(igraph::V(WWI.net)$name)
    pathwayadj <- Matrix::Matrix(0,ncol=N,nrow=N,sparse = TRUE)
    pathwayadj <-  igraph::as_adjacency_matrix(WWI.net,sparse = TRUE)
    pathwayadj <- pathwayadj[order(rownames(pathwayadj)),order(colnames(pathwayadj))]
    return(pathwayadj)

  }
}
