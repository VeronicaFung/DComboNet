#' Drug similarity based on Drug target distancce
#' @description The function \code{Targets.Dis} is to generate the distance
#'   matrix between drug targets
#' @param targets target genes list for drug candidates, target gene names
#'   should convert to official gene name (i.e. HGNC gene symbol)
#' @param PPI protein-protein interaction network. User can use either the
#'   default network from inBioMap PPI network or offer your own PPT network
#' @param load_dir path to load or save modelling data files
#' @return Distance matrix between target genes in protein-protein interaction
#'   network.
#' @export
#'


Targets.Dis <- function(targets,
                        PPI,
                        load_dir){

  #source("https://bioconductor.org/biocLite.R")
  #biocLite(c('supraHex','Rgraphviz','RBGL','graph','MASS', 'mgcv', 'nlme', 'rpart'))

  if(requireNamespace('igraph')){

  #Distance
  #Backgroud PPI: InBio

  PPI <- igraph::read.graph(paste0(load_dir, '/data/network/background.PPI/INBIO_PPI.net'), format='pajek')
  nodes <- read.csv(paste0(load_dir, '/data/network/background.PPI/INBIO_PPI_node.txt'), sep=' ',header=F, stringsAsFactors = F)
  targets <- unique(as.vector(sort(targets)))

  targets.index <- rep(0,length(targets))
  for(i in 1:length(targets)){
    if(targets[i] %in% nodes[,2])
      targets.index[i] <- nodes[nodes[,2]==targets[i],1]
  }
  Shortest.path <- igraph::shortest.paths(PPI,v=targets.index[targets.index!=0],to=targets.index[targets.index!=0])
  rownames(Shortest.path) <- targets[targets.index!=0]
  colnames(Shortest.path) <- targets[targets.index!=0]
  SP <- matrix(0,nrow=length(targets),ncol=length(targets))
  rownames(SP) <- targets
  colnames(SP) <- targets
  for(i in rownames(SP)){
    for(j in rownames(SP)){
      if(i %in% rownames(Shortest.path) && j %in% rownames(Shortest.path))
        SP[which(rownames(SP)==i),which(rownames(SP)==j)] <- Shortest.path[which(rownames(Shortest.path)==i),which(rownames(Shortest.path)==j)]
      else
        SP[which(rownames(SP)==i),which(rownames(SP)==j)] <- 0
    }
  }
  return(SP)
  }
}

#' Drug similarity based on Drug target distance
#' @description The function \code{DT.Dis} is to calculate the drug-drug
#'   similarity based on drug target distance in PPI network.
#' @param drug_target drug-target gene interaction file in \code{data.frame}
#'   format (otherwise \code{as.data.frame} can be used to transform the format
#'   of file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol)
#' @param PPI protein-protein interaction network. User can use either the
#'   default network from inBioMap PPI network or offer your own PPT network
#' @return Drug average target gene distance matrix
#' @export
#'

DT.Dis <- function(drug_target,
                   PPI){

  DT.M <- as.matrix(table(drug_target))
  targets <- as.vector(sort(drug_target[,2]))
  ITD.M <- Targets.Dis(targets, PPI)
  N.drug <- nrow(DT.M) #number of drug
  Dis.matrix <- matrix(0,nrow=N.drug,ncol=N.drug) #initiat the distance matrix
  rownames(Dis.matrix) <- rownames(DT.M)
  colnames(Dis.matrix) <- rownames(DT.M)

  for(i in 1:N.drug){
    for(j in i:N.drug){
      Dis.matrix[i,j] <- mean(ITD.M[which(DT.M[i,]>0),which(DT.M[j,]>0)])
      Dis.matrix[j,i] <- Dis.matrix[i,j]
    }
  }
  return(Dis.matrix)
}


#' Drug similarity based on pathway similarity
#' @description The function \code{WWI.sim} is to generate the distance matrix
#'   between drug targets
#' @param drugpair drug pairs generated from \code{druglist} you provided (or
#'   the subset)
#' @param drugtarget drug-target gene interaction file in \code{data.frame}
#'   format (otherwise \code{as.data.frame} can be used to transform the format
#'   of file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol)
#' @param targetPW gene and their related pathways table
#' @param PPI protein-protein interaction network. You can either use the
#'   default network from INBIO or offer your own PPT network.
#' @importFrom dplyr arrange
#' @import igraph
#' @return drug.matrix Drug-drug similarity matrix based on the pathway-pathway
#'   interaction
#' @export
#'

WWI.sim <- function(drugpair,
                    drugtarget,
                    targetPW,
                    PPI){

  if(requireNamespace(c('igraph','dplyr'))){

  #drugpair <- read.csv('new_data_integration/new_data/drug_drugpair/drugpair.new.csv', header = T, stringsAsFactors = F)
  #names(drugpair) = c('Drug1','Drug2')
  #drugtarget <- read.csv('new_data_integration/new_data/drugtarget/drugtarget.new.csv', header = T, stringsAsFactors = F)
  #targetPW <- unique(read.csv('new_data_integration/new_data/pathway_KEGG_DAVID/gene_pathway.csv', header = T, stringsAsFactors = F)[c(1,2)])


  #names(drugtarget) <- c('Drug1','Target1')
  #drugpair_target <- merge(drugpair,drugtarget,by = 'Drug1')
  #names(drugtarget) <- c('Drug2','Target2')
  #drugpair_target <- merge(drugpair_target,drugtarget,by = 'Drug2')
  #names(targetPW) <- c('Target1','Pathway1')
  #drugpair_target_wwi <- merge(drugpair_target,targetPW,by = 'Target1')
  #names(targetPW) <- c('Target2','Pathway2')
  #drugpair_target_wwi <- merge(drugpair_target_wwi,targetPW,by = 'Target2')
  #drugpair_target_wwi <- drugpair_target_wwi[c(4,2,5,3,1,6)]
  #write.csv(drugpair_target_wwi,'testset_Oneil/data/KEGG/drugpair_target_wwi.csv', row.names = F)

  drugpair_target_wwi <- read.csv(paste0(load_dir, '/data/KEGG/drugpair_target_wwi.csv'), header =T, stringsAsFactors = F)
  WWI1 <- read.csv(paste0(load_dir,'/pathway/WWI.txt'), sep = '\t', header =T, stringsAsFactors = F)
  WWI <- igraph::graph.data.frame(WWI1)
  drug.wwi1 <- unique(c(drugpair_target_wwi$Pathway1, drugpair_target_wwi$Pathway2))
  wwilist <- unique(c(WWI1$Pathway1, WWI1$Pathway2))
  drug.wwi <- intersect(drug.wwi1, wwilist)
  PPI <- read.csv(paste0(load_dir,'/data/network/background.PPI/INBIO_PPI.txt'), sep = '\t', header =F, stringsAsFactors = F)
  ppi_node <- unique((c(PPI$V1,PPI$V2)))
  PPI <- igraph::graph.data.frame(PPI)
  targets1 <- sort(unique(as.vector(drugtarget[,2])))
  targets <- intersect(targets1, ppi_node)

  drugpair_wwi <- unique(drugpair_target_wwi[c(1,3,4,6)]) #drugpair and the pathway-pathway interaction between them
  drugpair_tti <- unique(drugpair_target_wwi[c(1,2,4,5)]) #drugpair and the target-target interaction between them
  drugpair_wwi = dplyr::arrange(drugpair_wwi,Drug1,Drug2,Pathway1,Pathway2)
  drugpair_tti = dplyr::arrange(drugpair_tti,Drug1,Drug2,Target1,Target2)
  drugpair <- unique(drugpair_target_wwi[c(1,4)])
  drugpair<- dplyr::arrange(drugpair, Drug1, Drug2)
  druglist <- unique(c(drugpair$Drug1,drugpair$Drug2))
  drug.matrix = matrix(0,nrow = length(druglist), ncol = length(druglist))
  row.names(drug.matrix) = druglist
  colnames(drug.matrix) = druglist
  n = length(druglist)

  for(i in 1:n){
    for(j in i:n){
      #i = 1
      #j = 14
      if(druglist[i] == druglist[j]){
        drug.matrix[i,j] = 0
      }else{
        DTW = dplyr::filter(drugpair_target_wwi, Drug1 == druglist[i], Drug2 == druglist[j])
        if(nrow(DTW) != 0){
          DTW1 = setdiff(unique(DTW$Pathway1),unique(DTW$Pathway2))
          DTW2 = setdiff(unique(DTW$Pathway2),unique(DTW$Pathway1))
          DW1 <- intersect(wwilist,DTW1)
          DW2 <- intersect(wwilist,DTW2)
          if(length(DW1) != 0 & length(DW2)!=0){
            DW_dis = igraph::shortest.paths(WWI, DW1, DW2)
          }
          DTW_2 = unique(DTW[DTW$Pathway1 == DTW$Pathway2,])[c(1:6)]
          if(nrow(DTW.2) != 0){
            DTW2.W1 = unique(DTW_2$Pathway1)
            DTW2.W2 = unique(DTW_2$Pathway2)
            T1 = length(unique(DTW_2$Target1))
            T2 = length(unique(DTW_2$Target2))
            W = sort(unique(c(DTW_2$Pathway1,DTW_2$Pathway2)))
            w.matrix <- matrix(0,nrow =length(W), ncol = length(W))
            row.names(w.matrix) = W
            colnames(w.matrix) = W
            for(k in 1:length(W)){
              DT = dplyr::filter(DTW_2,Pathway1 == W[k])[c(2,5)]
              DT1 <- intersect(ppi_node,unique(DT[,1]))
              DT2 <- intersect(ppi_node,unique(DT[,2]))
              if(nrow(DT) != 0){
                DT_dis = igraph::shortest.paths(PPI, DT1, DT2)
                w.matrix[k,k] = exp(-sum(DT_dis)/(T1 * T2))
              }
            }
            nDW1 <- length(union(DW1,DTW2.W1))
            nDW2 <- length(union(DW2,DTW2.W2))
            drug.matrix[i,j] = exp(-(sum(DW_dis) + sum(w.matrix))/(nDW1*nDW2))
            drug.matrix[j,i] = drug.matrix[i,j]
          }
          else{
            nDW1 <- length(DW1)
            nDW2 <- length(DW2)
            drug.matrix[i,j] = exp(-sum(DW_dis) /(nDW1*nDW2))
            drug.matrix[j,i] = drug.matrix[i,j]
          }
        }
      }
    }
  }
  return(drug.matrix)
  }
}

