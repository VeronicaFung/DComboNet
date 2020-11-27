#' Generate transition matrix from drug-gene/pathway interaction networks
#' @description The function \code{TransitionMatrix} is to generate the
#'   adjacency matrix from drug-drug association network, drug-gene association
#'   network, gene-gene association network, drug-pathway association network,
#'   pathway-pathway association network and gene-pathway association network.
#' @param drugAdj drug adjacency matrix (sparse matrix) generated from drug-drug
#'   association network via \code{AdjMatrix_drugs}
#' @param geneAdj adjacency matrix of gene-gene association network (sparse
#'   matrix) via \code{geneNet.Adj}
#' @param pathwayAdj adjacency matrix of pathway-pathway association network
#'   (sparse matrix) via \code{pathwayNet.Adj}
#' @param druggeneAdj adjacency matrix of drug-gene association network (sparse
#'   matrix) via \code{drugGeneNet.L1} or \code{drugGeneNet.L2} depends on which
#'   model is called
#' @param drugpathwayAdj adjacency matrix of drug-pathway association network
#'   (sparse matrix) via \code{DrugPathwayMatrix.L1} or
#'   \code{DrugPathwayMatrix.L2} depends on which model is called
#' @param genepathwayAdj adjacency matrix of gene-pathway association network
#'   (sparse matrix) via \code{genepathwayAdj}
#' @param x numeric, the probability that the inter-network crossing event
#'   happens between drug and gene nodes
#' @param y numeric, the probability that the inter-network crossing event
#'   happens between drug and pathway nodes
#' @param z numeric, the probability that the inter-network crossing event
#'   happens between gene and pathway nodes
#' @param A numeric, the probability of jumping within drugnet, when the
#'   inter-network crossing events between drug and gene nodes and between drug
#'   and pathway nodes are neither existing
#' @param B numeric, the probability of jumping within genenet, when the
#'   inter-network crossing events between drug and gene nodes and between gene
#'   and pathway nodes are neither existing
#' @import Matrix
#' @return Transition matrix generated from drug-gene/pathway interaction
#'   networks
#' @export
#'
#' @seealso Functions to generate adjacency matrix for different subnetworks:
#'   \code{\link{AdjMatrix_drugs}} for drug-drug association network,
#'   \code{\link{geneNet.Adj}} for gene-gene association network,
#'   \code{\link{pathwayNet.Adj}} for pathway-pathway association network,
#'   \code{\link{drugGeneNet.L1}} or \code{\link{drugGeneNet.L2}} for drug-gene
#'   association network, \code{\link{DrugPathwayMatrix.L1}} or
#'   \code{\link{DrugPathwayMatrix.L2}} for drug-pathway association network,
#'   \code{\link{genepathwayAdj}} for gene-pathway association network.
#'
#' @examples
#' \dontrun{
#' # This function replies on adjacency matrices from different subnetworks
#'   drugAdj = AdjMatrix_drugs(...)
#'   geneAdj = geneNet.Adj(...)
#'   pathwayAdj = pathwayNet.Adj(...)
#'   # For level one model
#'   druggeneAdj = drugGeneNet.L1(...)
#'   drugpathwayAdj = DrugPathwayMatrix.L1(...)
#'   # For level two model
#'   druggeneAdj = drugGeneNet.L2(...)
#'   drugpathwayAdj = DrugPathwayMatrix.L2(...)
#'   genepathwayAdj = genepathwayAdj(...)
#'   x = 0.5
#'   y = 0.5
#'   z = 0
#'   A = 0.5
#'   B = 0.5
#'   TransitionMatrix(drugAdj = drugAdj, geneAdj = geneAdj, pathwayAdj =
#'   pathwayAdj, druggeneAdj = druggeneAdj, drugpathwayAdj = drugpathwayAdj,
#'   genepathwayAdj = genepathwayAdj, x = 0.5, y = 0.5, z = 0, A = 0.5, B = 0.5)
#' }

# Transition matrix
TransitionMatrix <- function(drugAdj,
                             geneAdj,
                             pathwayAdj,
                             druggeneAdj,
                             drugpathwayAdj,
                             genepathwayAdj,
                             x = 0.5,
                             y = 0.5,
                             z = 0,
                             A = 0.5,
                             B = 0.5){

  if (requireNamespace("Matrix")){

    # drugAdj = drugnet_adj
    # geneAdj = geneadj
    # pathwayAdj = pathwayadj
    # druggeneAdj = dgAdj
    # genepathwayAdj = gpAdj
    # drugpathwayAdj = drugpathwayMatrix
    #x = 1/3
    #y = 1/3
    #z = 1/3

    M = nrow(drugAdj)
    N = nrow(geneAdj)
    Transition_Gene_Drugs <- Matrix::Matrix(0,nrow=N,ncol=M,sparse = TRUE)
    colnames(Transition_Gene_Drugs) <- colnames(druggeneAdj)
    row.names(Transition_Gene_Drugs) <- row.names(druggeneAdj)
    Col_Sum_Bipartite <- Matrix::colSums(druggeneAdj, na.rm = FALSE, dims = 1)#, sparseResult = FALSE)
    for(j in 1:M){
      if(Col_Sum_Bipartite[j] != 0){
        Transition_Gene_Drugs[,j] <- (x * druggeneAdj[,j] /Col_Sum_Bipartite[j])
      }else if(Col_Sum_Bipartite[j] == 0){
        Transition_Gene_Drugs[,j] <- 0
      }
    }

    #write.csv(as.matrix(Transition_Gene_Drugs,"Output_Files/test.file/Transition_Gene_Drugs.csv")

    Transition_Drugs_Gene <- Matrix::Matrix(0,nrow=M,ncol=N,sparse = TRUE)
    colnames(Transition_Drugs_Gene) <- rownames(druggeneAdj)
    rownames(Transition_Drugs_Gene) <- colnames(druggeneAdj)
    Row_Sum_Bipartite <- Matrix::rowSums(druggeneAdj, na.rm = FALSE, dims = 1)#, sparseResult = FALSE)
    for(i in 1:N){
      if(Row_Sum_Bipartite[i] != 0){
        Transition_Drugs_Gene[,i] <- (x * druggeneAdj[i,] /Row_Sum_Bipartite[i])
      }else if(Row_Sum_Bipartite[i] == 0){
        Transition_Drugs_Gene[,i] <- 0
      }
    }
    #write.csv(as.matrix(Transition_Drugs_Gene),"Output_Files/test.file/Transition_Drugs_Gene.csv")

    #drugAdj = drugnet_adj
    #geneAdj = pathwayadj
    M = nrow(drugAdj)
    L = nrow(pathwayAdj)
    Transition_Pathway_Drugs <- Matrix::Matrix(0,nrow=L,ncol=M,sparse = TRUE)
    colnames(Transition_Pathway_Drugs) <- colnames(drugpathwayAdj)
    row.names(Transition_Pathway_Drugs) <- row.names(drugpathwayAdj)
    Col_Sum_Bipartite <- Matrix::colSums(drugpathwayAdj, na.rm = FALSE, dims = 1)
    for(j in 1:M){
      if(Col_Sum_Bipartite[j] != 0){
        Transition_Pathway_Drugs[,j] <- (y*drugpathwayAdj[,j] /Col_Sum_Bipartite[j])
      }else if(Col_Sum_Bipartite[j] == 0){
        Transition_Pathway_Drugs[,j] <- 0
      }
    }

    #write.csv(as.matrix(Transition_Gene_Drugs,"Output_Files/test.file/Transition_Gene_Drugs.csv")

    Transition_Drugs_Pathway <- Matrix::Matrix(0,nrow=M,ncol=L,sparse = TRUE)
    colnames(Transition_Drugs_Pathway) <- row.names(drugpathwayAdj)
    row.names(Transition_Drugs_Pathway) <- colnames(drugpathwayAdj)

    Row_Sum_Bipartite <- Matrix::rowSums(drugpathwayAdj, na.rm = FALSE, dims = 1)#, sparseResult = FALSE)
    for(i in 1:L){
      if(Row_Sum_Bipartite[i] != 0){
        Transition_Drugs_Pathway[,i] <- (y*drugpathwayAdj[i,] /Row_Sum_Bipartite[i])
      }else if(Row_Sum_Bipartite[i] == 0){
        Transition_Drugs_Pathway[,i] <- 0
      }
    }

    #geneAdj = geneadj
    #pathwayAdj = pathwayadj
    N = nrow(geneAdj)
    L = nrow(pathwayAdj)
    Transition_Pathway_Genes <- Matrix::Matrix(0,nrow=L,ncol=N,sparse = TRUE)
    colnames(Transition_Pathway_Genes) <- colnames(genepathwayAdj)
    row.names(Transition_Pathway_Genes) <- row.names(genepathwayAdj)
    Col_Sum_Bipartite <- Matrix::colSums(genepathwayAdj)
    for(j in 1:N){
      if(Col_Sum_Bipartite[j] != 0){
        Transition_Pathway_Genes[,j] <- (z*genepathwayAdj[,j] /Col_Sum_Bipartite[j])
      }else if(Col_Sum_Bipartite[j] == 0){
        Transition_Pathway_Genes[,j] <- 0
      }
    }

    #write.csv(as.matrix(Transition_Gene_Drugs,"Output_Files/test.file/Transition_Gene_Drugs.csv")

    Transition_Genes_Pathway <- Matrix::Matrix(0,nrow=N,ncol=L,sparse = TRUE)
    colnames(Transition_Genes_Pathway) <- row.names(genepathwayAdj)
    row.names(Transition_Genes_Pathway) <- colnames(genepathwayAdj)
    Row_Sum_Bipartite <- Matrix::rowSums(genepathwayAdj)
    for(i in 1:L){
      if(Row_Sum_Bipartite[i] != 0){
        Transition_Genes_Pathway[,i] <- (z*genepathwayAdj[i,]) /Row_Sum_Bipartite[i]
      }else if(Row_Sum_Bipartite[i] == 0){
        Transition_Genes_Pathway[,i] <- 0
      }
    }




    Transition_Pathway_Network <- Matrix::Matrix(0,nrow=L,ncol=L,sparse = TRUE)
    rownames(Transition_Pathway_Network) <- rownames(pathwayAdj)
    colnames(Transition_Pathway_Network) <- colnames(pathwayAdj)

    Col_Sum_Multiplex <- Matrix::colSums(pathwayAdj,na.rm=FALSE,dims=1, sparseResult=FALSE)
    Row_Sum_Bipartite <- Matrix::rowSums(drugpathwayAdj, na.rm = FALSE, dims = 1)#, sparseResult = FALSE)

    for (j in 1:L){
      if(Row_Sum_Bipartite[j] != 0){
        Transition_Pathway_Network[,j] <- ((1-y)*pathwayAdj[,j] /Col_Sum_Multiplex[j])
      } else {
        Transition_Pathway_Network[,j] <- pathwayAdj[,j] /Col_Sum_Multiplex[j]
      }
    }


    Transition_Gene_Network <- Matrix::Matrix(0,nrow=N,ncol=N,sparse = TRUE)
    rownames(Transition_Gene_Network) <- rownames(geneAdj)
    colnames(Transition_Gene_Network) <- colnames(geneAdj)

    Col_Sum_Multiplex <- Matrix::colSums(geneAdj,na.rm=FALSE,dims=1, sparseResult=FALSE)
    Row_Sum_Bipartite1 <- Matrix::rowSums(druggeneAdj, na.rm = FALSE, dims = 1) #, sparseResult = FALSE)
    Row_Sum_Bipartite2 <- Matrix::colSums(genepathwayAdj)

    for (j in 1:N){
      if(Row_Sum_Bipartite1[j] != 0 & Row_Sum_Bipartite2[j] != 0){

        Transition_Gene_Network[,j] <- (B*geneAdj[,j]/Col_Sum_Multiplex[j])

      } else if(Row_Sum_Bipartite1[j] != 0 & Row_Sum_Bipartite2[j] == 0){

        Transition_Gene_Network[,j] <- ((1-x)*geneAdj[,j] /Col_Sum_Multiplex[j])

      } else if(Row_Sum_Bipartite1[j] == 0 & Row_Sum_Bipartite2[j] != 0){

        Transition_Gene_Network[,j] <- ((1-z)*geneAdj[,j] /Col_Sum_Multiplex[j])

      } else{
        Transition_Gene_Network[,j] <- geneAdj[,j] /Col_Sum_Multiplex[j]
      }
    }

    Transition_Drugs_Network <- Matrix::Matrix(0,nrow=M,ncol=M,sparse = TRUE)

    rownames(Transition_Drugs_Network) <- rownames(drugAdj)
    colnames(Transition_Drugs_Network) <- colnames(drugAdj)

    Col_Sum_Drug <- Matrix::colSums(drugAdj,na.rm=FALSE,dims=1, sparseResult=FALSE)
    Col_Sum_Bipartite1 <- Matrix::colSums(druggeneAdj, na.rm = FALSE, dims = 1)#, sparseResult = FALSE)
    Col_Sum_Bipartite2 <- Matrix::colSums(drugpathwayAdj, na.rm = FALSE, dims = 1)#, sparseResult = FALSE)

    for (j in 1:M){
      if(Col_Sum_Bipartite1[j] != 0 & Col_Sum_Bipartite2[j] !=0){

        Transition_Drugs_Network[,j] <- (A*drugAdj[,j] /Col_Sum_Drug[j])

      }else if(Col_Sum_Bipartite1[j] != 0 & Col_Sum_Bipartite2[j] ==0 ){

        Transition_Drugs_Network[,j] <- ((1-x)*drugAdj[,j]/Col_Sum_Drug[j])

      } else if(Col_Sum_Bipartite1[j] == 0  &  Col_Sum_Bipartite2[j] != 0 ){

        Transition_Drugs_Network[,j] <- ((1-y)*drugAdj[,j]/Col_Sum_Drug[j])

      }else{

        Transition_Drugs_Network[,j] <- drugAdj[,j] /Col_Sum_Drug[j]

      }
    }



    Transition_Gene_Network[( is.na(Transition_Gene_Network))] = 0
    Transition_Gene_Drugs[( is.na(Transition_Gene_Drugs))] = 0
    Transition_Drugs_Gene[( is.na(Transition_Drugs_Gene))] = 0
    Transition_Drugs_Network[( is.na(Transition_Drugs_Network))] = 0
    Transition_Pathway_Network[( is.na(Transition_Pathway_Network))] = 0
    Transition_Genes_Pathway[( is.na(Transition_Genes_Pathway))] = 0
    Transition_Pathway_Genes[( is.na(Transition_Pathway_Genes))] = 0
    Transition_Drugs_Pathway[( is.na(Transition_Drugs_Pathway))] = 0
    Transition_Pathway_Drugs[( is.na(Transition_Pathway_Drugs))] = 0


    Transition_Matrix_1 <- cbind(Transition_Drugs_Network,Transition_Drugs_Gene, Transition_Drugs_Pathway)
    Transition_Matrix_2 <- cbind(Transition_Gene_Drugs,Transition_Gene_Network, Transition_Genes_Pathway)
    Transition_Matrix_3 <- cbind(Transition_Pathway_Drugs, Transition_Pathway_Genes,Transition_Pathway_Network)

    Transition_Matrix <- rbind(Transition_Matrix_1,Transition_Matrix_2,Transition_Matrix_3)
    return(Transition_Matrix)
  }
}



