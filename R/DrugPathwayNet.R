#  generation drug-pathway bipartite matrix

#' Generate adjacency matrix from drug-pathway interaction networks for level
#' one model
#' @description The function "DrugPathwayMatrix.L1" is to generate the adjacency
#'   matrix from drug-pathway interaction network for level one model
#' @param dt drug-target gene interaction file in \strong{data.frame} format
#'   (otherwise \code{as.data.frame} can be used to transform the format of
#'   file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol).
#' @param dg drug-gene interaction file (Optional). drug and genes that either
#'   induced by drug or may influence drug effect can be provided for construct
#'   gene network for level one model. Genes can be collected from publications
#'   or from IPA tool. Gene name should be converted to official gene names
#'   (i.e. HGNC gene symbols).
#' @param drugAdj drug adjacency matrix (sparse matrix) generated from drug-drug
#'   association network via \code{AdjMatrix_drugs}
#' @param PathwayNetwork \code{data.frame} of pathway pathway association network
#'   network, sparse matrix generated via \code{pathwayNet.Adj}
#' @param dgpweight numeric, weight of drug-pathway associations through drug
#'   related genes
#' @param load_dir path to load or save modeling data files
#' @import Matrix
#' @import igraph
#' @import Hmisc
#' @return \code{drugpathwayMatrix} drug-pathway adjacency matrix (sparse
#'   matrix)
#'
#' @export
#'
#' @examples
#' #
#' # loading files
#' \dontrun{
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' drugnet_feature = read.table(paste0(load_dir,"drug_net/features.csv"),
#' sep=",",header=TRUE, stringsAsFactors = FALSE)
#' drugtarget = read.table(paste0(load_dir,"data/drugtarget.csv"), sep =
#' ",",header = TRUE, stringsAsFactors = FALSE)
#' druggene = data.frame(drug = "Sorafenib", target  = "RAF1")
#' # or druggene = NULL
#' drugnet_adj = AdjMatrix_drugs(x = drugnet_feature, weighted = TRUE,
#' weight.type = "integrated_score")
#' PathwayNetwork = read.table(paste0(load_dir,"pathway/WWI.txt"), sep = "\t",
#' header =TRUE, stringsAsFactors = FALSE)
#' dgpweight = 1
#' #
#' # construct drug-pathway adjacency matrix for level one model#'
#' DrugPathwayMatrix.L1(dt = drugtarget, dg = druggene, drugAdj = drugnet_adj,
#' dgpweight = dgpweight, PathwayNetwork = PathwayNetwork, load_dir = load_dir)
#' }
#'


DrugPathwayMatrix.L1 <- function(dt = NULL,
                                 dg = NULL,
                                 drugAdj,
                                 PathwayNetwork,
                                 dgpweight,
                                 load_dir){

  if (requireNamespace(c("Matrix","Hmisc","igraph"))){

    # load(system.file("data", "drugtarget_lib.Rdata", package = "DComboNet"))
    # load(system.file("data", "druggene_lib.Rdata", package = "DComboNet"))
    # load(system.file("data", "KEGG_pw_gene.Rdata", package = "DComboNet"))
    # dt_lib = read.csv(paste0(load_dir,"data/drugtarget.csv"),stringsAsFactors = F)

    names(dt_lib) = c("Drug","Target")
    if(is.null(dt)){
      dt = dt_lib
    }else{
      names(dt) = c("Drug","Target")
      dt = unique(rbind(dt_lib,dt))
    }

    if(is.null(dg)){
      print("       No new drug related genes received       ")
      print("       using only drug related genes In Our Library       ")

      # dg_lib = unique(read.csv(paste0(load_dir,"data/drug_gene/drug_ipa.inPPI.csv"), header = T, stringsAsFactors = F)[c(1,3)])
      names(dg_lib) = c("Drug","Target")
      dg = unique(dg_lib)

    }else{
      print("       New drug related genes have been found       ")
      print("       Merge new input with default drug related genes       ")
      # dg_lib = unique(read.csv(paste0(load_dir,"data/drug_gene/drug_ipa.inPPI.csv"), header = T, stringsAsFactors = F)[c(1,3)])
      names(dg_lib) = c("Drug","Target")
      names(dg) = c("Drug","Target")
      dg = unique(rbind(dg_lib,dg))
    }

    # genepathway = read.table(paste0(load_dir,"pathway/KEGG_ID_gene.txt"),sep="\t",header = T,stringsAsFactors = F)

    genepathway = genepathway[c(3,1)]
    names(genepathway) = c("Target","pathway")
    dtw = merge(dt, genepathway,by="Target")
    dw = unique(dtw[c(2,3)])

    dg_pw = merge(dg, genepathway, by = "Target")
    dg_pw =  dg_pw[c(2,3)]
    names(dg_pw) = names(dw)
    dg_pw$Drug = Hmisc::capitalize(dg_pw$Drug)
    names(dg_pw) = c("Drug","Pathway")

    WWI.net <- igraph::graph.data.frame(PathwayNetwork)
    pathways <- sort(igraph::V(WWI.net)$name)
    N = length(pathways)
    M = nrow(drugAdj)
    names(dw) = c("Drug","Pathway")
    drugpathwayMatrix = Matrix::Matrix(data = 0, nrow = N, ncol = M)
    rownames(drugpathwayMatrix) = pathways
    colnames(drugpathwayMatrix) = colnames(drugAdj)

    for(i in 1:N){
      pathway <- pathways[i]
      drug <- dg_pw[which(dg_pw$Pathway == pathway),]$Drug
      if(length(drug) > 0){
        for(j in 1:length(drug)){
          if(!is.na(drug[j])){
            # identify the drug position on the matrix.
            index_drug <- which(colnames(drugpathwayMatrix) %in% drug[j])
            # check if that index is present in the matrix.
            if(length(index_drug) == 1){
              drugpathwayMatrix[i,index_drug] <- dgpweight #
            }
          }
        }
      }
    }


    for(i in 1:N){
      pathway <- pathways[i]
      drug <- dw[which(dw$Pathway == pathway),]$Drug

      if(length(drug) > 0){
        for(j in 1:length(drug)){
          if(!is.na(drug[j])){
            # identify the drug position on the matrix.
            index_drug <- which(colnames(drugpathwayMatrix) %in% drug[j])
            # check if that index is present in the matrix.
            if(length(index_drug) == 1){
              drugpathwayMatrix[i,index_drug] <- length(dtw[dtw$Drug == drug & dtw$pathway == pathway,]$Target)
            }
          }
        }
      }
    }

    return(drugpathwayMatrix)
  }
}



#  generation drug-pathway bipartite matrix

#' Generate the bipartite matrix from drug-pathway interaction networks
#' @description The function \code{DrugPathwayMatrix.L2} is to generate
#'   adjacency matrix from drug-pathway interaction network for level two model.
#' @param dt drug-target gene interaction file in \strong{data.frame} format
#'   (otherwise \code{as.data.frame} can be used to transform the format of
#'   file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol).
#' @param dDEG if \code{cellline} is not included in provided cancer cell line
#'   list, in \code{as.data.frame} format must be provided, otherwise
#'   \emph{ERROR} information will be returned
#' @param drugDEP if \code{cellline} is not included in provided cancer cell
#'   line list, in \code{as.data.frame} format must be provided, otherwise
#'   \emph{ERROR} information will be returned
#' @param drugAdj drug adjacency matrix (sparse matrix) generated from drug-drug
#'   association network via \code{AdjMatrix_drugs}
#' @param PathwayNetwork \code{data.frame} of pathway pathway association
#'   network network, sparse matrix generated via \code{pathwayNet.Adj}
#' @param dataset GEO ID of LINCS datasets, two datasets are provided here:
#'   "70138" and "92742"
#' @param cellline cancer cell lines name. If level two model is called
#'   (\code{model = "L2"} in \code{DComboNet} or \code{L2.LOCOCV}). If requested
#'   cell line is included in the default cancer cell line list, \code{dDEG} do
#'   not need to provided unless particular drug induced differentially
#'   expressed gene list you want to include in the gene network, then provided
#'   \code{dDEG} in \code{as.data.frame} format. If the cancer cell line is not
#'   in the current cell line list, you should provide \code{cellline} and
#'   \code{dDEG} for construct specific gene network, otherwise \emph{ERROR}
#'   information will be returned
#' @param treatment_time drug treatment time provided in LINCS database. For
#'   example, MCF7 cell line (\code{cellline = "MCF7"}) provided two time points
#'   in LINCS GSE70138 dataset - 6 hour and 12 hour - \code{6} or \code{12}
#'   should be taken as input here, as numeric or character
#' @param dDEpweight numeric, weight of drug-drug-differentially expressed
#'   pathway edges
#' @param load_dir path to load or save modeling data files
#' @import Matrix
#' @import igraph
#' @import Hmisc
#' @examples
#' #
#' # loading files
#' \dontrun{
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' drugnet_feature = read.table(paste0(load_dir,"drug_net/features.csv"),
#' sep=",",header=TRUE, stringsAsFactors = FALSE)
#' drugnet_feature = drugnet_feature[drugnet_feature$integrated_score2 >=0.2,]
#' drugtarget = read.table(paste0(load_dir,"data/drugtarget.csv"), sep =
#' ",",header = TRUE, stringsAsFactors = FALSE)
#' drugDEG = NULL
#' # or drugDEG = data.frame(drug = "Sorafenib", target = c("FOSL1", "TUBB6"))
#' dataset = "92742"
#' cellline = "HEPG2"
#' t="6"
#' drugnet_adj = AdjMatrix_drugs(x = drugnet_feature, weighted = TRUE,
#' weight.type = "integrated_score")
#' PathwayNetwork = read.table(paste0(load_dir,"pathway/WWI.txt"), sep = "\t",
#' header =TRUE, stringsAsFactors = FALSE)
#' dDEpweight = 1
#' #
#' # construct drug-pathway adjacency matrix for level two model
#' DrugPathwayMatrix.L2(dt = drugtarget, dDEG = drugDEG, dataset=dataset,
#' cellline = cellline, treatment_time=t, drugAdj = drugnet_adj, dDEpweight =
#' dDEpweight, PathwayNetwork = PathwayNetwork, load_dir = load_dir)
#' }
#'
#' @return \code{drugpathwayMatrix} drug-pathway adjacency matrix (sparse
#'   matrix)
#' @export
#'


DrugPathwayMatrix.L2 <- function(dt = NULL,
                                 dg = NULL,
                                 dDEG = NULL,
                                 dataset=NULL,
                                 cellline,
                                 treatment_time = NULL,
                                 drugDEP = NULL,
                                 drugAdj,
                                 PathwayNetwork,
                                 dDEpweight = 2,
                                 load_dir){

  if (requireNamespace(c("Matrix","Hmisc","igraph"))){

    # dt_lib = read.csv(paste0(load_dir,"data/drugtarget.csv"),stringsAsFactors = F)
    # load(system.file("data", "drugtarget_lib.Rdata", package = "DComboNet"))
    # load(system.file("data", "druggene_lib.Rdata", package = "DComboNet"))
    # load(system.file("data", "KEGG_pw_gene.Rdata", package = "DComboNet"))
    # load(system.file("data", "cellline_name.Rdata", package = "DComboNet"))
    # genepathway = read.table(paste0(load_dir,"pathway/KEGG_ID_gene.txt"),sep="\t",header = T,stringsAsFactors = F)
    genepathway = genepathway[c(3,1)]
    names(genepathway) = c("Target","pathway")

    names(dt_lib) = c("Drug","Target")

    if(is.null(dt)){
      dt = dt_lib
    }else{
      names(dt) = c("Drug","Target")
      dt = unique(rbind(dt_lib,dt))
    }

    dtw = merge(dt, genepathway,by="Target")
    # dw = unique(dtw[c(2,3)])
    # names(dw) = c("Drug","Pathway")
    dw_t=data.frame(table(dtw[c(2,3)]))
    dw = dw_t[dw_t$Freq > 0,]
    names(dw) = c('Drug','Pathway','Target')

    if(is.null(dg)){
      print("       No new drug related genes received       ")
      print("       using only drug related genes In Our Library       ")

      # dg_lib = unique(read.csv(paste0(load_dir,"data/drug_gene/drug_ipa.inPPI.csv"), header = T, stringsAsFactors = F)[c(1,3)])
      names(dg_lib) = c("Drug","Target")
      dg = unique(dg_lib)

    }else{
      print("       New drug related genes have been found       ")
      print("       Merge new input with default drug related genes       ")
      # dg_lib = unique(read.csv(paste0(load_dir,"data/drug_gene/drug_ipa.inPPI.csv"), header = T, stringsAsFactors = F)[c(1,3)])
      names(dg_lib) = c("Drug","Target")
      names(dg) = c("Drug","Target")
      dg = unique(rbind(dg_lib,dg))
    }

    # celllineDEP_lib = read.csv(paste0(load_dir,"data/celllineDEG_lib.csv"),stringsAsFactors = F)
    dg_pw = merge(dg, genepathway, by = "Target")
    # dg_pw =  dg_pw[c(2,3)]
    # names(dg_pw) = c("Drug","Pathway")
    # dg_pw$Drug = Hmisc::capitalize(dg_pw$Drug)
    dg_pw_t=data.frame(table(dg_pw[c(2,3)]))
    dg_pw = dg_pw_t[dg_pw_t$Freq > 2,]
    names(dg_pw) = c('Drug','Pathway','Target')

    if(cellline %in% celllineDEP_lib[,1]){

      if(is.null(drugDEP)){

        dp_lib.path =  system.file("extdata", paste0("LINCS_data/GSE",dataset,"/",cellline,"_DEgene_",treatment_time,"h.csv"), package = "DComboNet")
        if(file.exists(dp_lib.path)){
          d_pw <- unique(read.csv(dp_lib.path, header = T, stringsAsFactors = F))
        }

        # d_pw = unique(read.csv(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEgene_",treatment_time,"h.csv"))[c(2,4,3)])
        d_pw = unique(d_pw[c(2,4,3)])
        names(d_pw) = c("Drug", "Pathway", "logFC")
        d_pw$Drug = Hmisc::capitalize(d_pw$Drug)
        druglist = data.frame(Drug = colnames(drugAdj))
        d_pw = merge(druglist, d_pw, by = "Drug")

      }else{

        dp_lib.path =  system.file("extdata", paste0("LINCS_data/GSE",dataset,"/",cellline,"_DEgene_",treatment_time,"h.csv"), package = "DComboNet")
        # d_pw = unique(read.csv(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEgene_",treatment_time,"h.csv"))[c(2,4,3)])
        if(file.exists(dp_lib.path)){
          d_pw <- unique(read.csv(dp_lib.path, header = T, stringsAsFactors = F))
        }
        d_pw = unique(d_pw[c(2,4,3)])

        names(d_pw) = c("Drug","Pathway","logFC")
        d_pw$Drug = Hmisc::capitalize(d_pw$Drug)

        druglist = data.frame(Drug = colnames(drugAdj))
        d_pw = merge(druglist, d_pw,by="Drug")

        if(("logFC" %in% names(drugDEP)) != TRUE){
          drugDEP$logFC = log2(dDEpweight)
          names(drugDEP) = c("Drug","Target","logFC")

        }

        d_pw = unique(rbind(d_pw,drugDEP))
        druglist = data.frame(Drug = colnames(drugAdj))
        d_pw = merge(druglist, d_pw,by="Drug")
      }
      dDEG_pw = d_pw
      names(dDEG_pw) = c('Drug','Pathway','Target')
      dDEG_pw$Target = 1
      dDEG_pw = unique(dDEG_pw)
    }else{

      if(is.null(drugDEP)){
        if(is.null(dDEG)){
          print("       ERROR: Calling Level 2 Model, Drug-DEG should be provided if Drug-DEP File are not available       ")
        }else{
          #
          # if(("logFC" %in% names(drugDEG)) != TRUE){
          #   drugDEP$logFC = log2(dDEpweight)
          # }
          # names(drugDEP) = c("Drug","Target","logFC")
          # d_pw = drugDEP
          druglist = data.frame(Drug = unique(dt$Drug))
          # druglist = data.frame(Drug=colnames(drugAdj))# unique(dt[1])
          names(dDEG) = names(dt)
          # dDEG$Drug = capitalize(dDEG$Drug)
          dDEG = merge(dDEG,druglist, by='Drug')
          names(dDEG) = c("Drug","Target","logFC")
          dDEG_pw = unique(merge(dDEG, genepathway, by = 'Target'))
          tmp=data.frame(table(dDEG_pw[c(2,4)]))
          dDEG_pw = tmp[tmp$Freq > 2,]
          #dDEG_pw =  unique(dDEG_pw[c(2,3,1)])
          #names(dDEG_pw) = names(dw)
          names(dDEG_pw) = c('Drug','Pathway','Target')

          #dDEG = dDEG[c('Drug','Target','logFC')]
        }
      }
    }


    druglist = data.frame(Drug = unique(dt$Drug))
    # druglist = data.frame(Drug=colnames(drugAdj))# unique(dt[1])
    dDEG_lib.path =  system.file("extdata", paste0("LINCS_data/GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), package = "DComboNet")

    if(file.exists(dDEG_lib.path)){
      dDEG_lib <- unique(read.csv(dDEG_lib.path, header = T, stringsAsFactors = F))
    }
    # dDEG_lib <- unique(read.csv(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), header = T, stringsAsFactors = F))

    names(dDEG_lib) = c("Drug","Target","logFC")
    dDEG_lib$Drug = Hmisc::capitalize(dDEG_lib$Drug)
    druglist = unique(dt[1])
    dDEG = merge(dDEG_lib,druglist, by="Drug")
    dDEG = dDEG[c("Drug","Target","logFC")]

    dDEG_pw = unique(merge(dDEG, genepathway, by = 'Target'))
    tmp=data.frame(table(dDEG_pw[c(2,4)]))
    dDEG_pw = tmp[tmp$Freq > 2,]
    #dDEG_pw =  unique(dDEG_pw[c(2,3,1)])
    #names(dDEG_pw) = names(dw)
    names(dDEG_pw) = c('Drug','Pathway','Target')

    WWI.net <- igraph::graph.data.frame(PathwayNetwork)
    pathways <- sort(igraph::V(WWI.net)$name)
    N = length(pathways)
    M = nrow(drugAdj)
    drugpathwayMatrix = Matrix::Matrix(data = 0, nrow = N, ncol = M)
    rownames(drugpathwayMatrix) = pathways
    colnames(drugpathwayMatrix) = colnames(drugAdj)
    #dtw <- dtw[which(dtw$pathway %in% genes), ]

    for(i in 1:N){
      pathway <- pathways[i]
      drug <- dg_pw[which(dg_pw$Pathway == pathway),]$Drug
      if(length(drug) > 0){
        for(j in 1:length(drug)){
          if(!is.na(drug[j])){
            # identify the drug position on the matrix.
            index_drug <- which(colnames(drugpathwayMatrix) %in% drug[j])
            # check if that index is present in the matrix.
            if(length(index_drug) == 1){
              drugpathwayMatrix[i,index_drug] <- dg_pw[dg_pw$Drug == drug[j] & dg_pw$Pathway == pathway,]$Target # / max(dg_pw$Target) *1
            }
          }
        }
      }
    }

    for(i in 1:N){
      pathway <- pathways[i]
      drug <- dDEG_pw[which(dDEG_pw$Pathway == pathway),]$Drug
      if(length(drug) > 0){
        for(j in 1:length(drug)){
          if(!is.na(drug[j])){
            # identify the drug position on the matrix.
            index_drug <- which(colnames(drugpathwayMatrix) %in% drug[j])
            # check if that index is present in the matrix.
            if(length(index_drug) == 1){
              drugpathwayMatrix[i,index_drug] <- dDEG_pw[dDEG_pw$Drug == drug[j] & dDEG_pw$Pathway == pathway,]$Target # / max(dDEG_pw$Target) *2# dDEGpWeight
            }
          }
        }
      }
    }

    #dtw <- dtw[which(dtw$pathway %in% genes), ]
    for(i in 1:N){
      pathway <- pathways[i]
      drug <- dw[which(dw$Pathway == pathway),]$Drug

      if(length(drug) > 0){
        for(j in 1:length(drug)){
          if(!is.na(drug[j])){
            # identify the drug position on the matrix.
            index_drug <- which(colnames(drugpathwayMatrix) %in% drug[j])
            # check if that index is present in the matrix.
            if(length(index_drug) == 1){
              drugpathwayMatrix[i,index_drug] <- dw[dw$Drug == drug[j] & dw$Pathway == pathway,]$Target # / max(dw$Target) *3
            }
          }
        }
      }
    }
    return(drugpathwayMatrix)
  }
}


