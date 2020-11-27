#' Generate adjacency matrix from drug-gene interaction networks
#' for level one model
#' @description The function \code{drugGeneAdj} is to generate the adjacency
#'   matrix from gene-gene interaction network
#' @param dt drug-target gene interaction file in \strong{data.frame} format
#'   (otherwise \code{as.data.frame} can be used to transform the format of
#'   file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol).
#' @param dg drug-gene interaction file (Optional). drug and genes that either
#'   induced by drug or may influence drug effect can be provided for construct
#'   gene network for level one model. Genes can be collected from publications
#'   or from IPA tool. Gene name should be converted to official gene names
#'   (i.e. HGNC gene symbols).
#' @param drugAdj drug adjacency matrix generated from drug-drug
#'   association network via \code{AdjMatrix_drugs}
#' @param geneNet gene-gene association network adjacency matrix (\bold{sparse}
#'   matrix format)
#' @param dtweight the weight of drug-target gene edges
#' @param dgweight the weight of drug-related gene edge
#' @param load_dir path to load or save modeling data files
#' @import Matrix
#' @import igraph
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
#' druggene = data.frame(drug = "Sorafenib", target  = "RAF1")
#' drugnet_adj = AdjMatrix_drugs(x = drugnet_feature, weighted = TRUE,
#' weight.type = "integrated_score")
#' gene.net = geneNet(dt = drugtarget, dg = druggene, load_dir = load_dir)
#' dtweight = 2
#' dgweight = 1
#' #
#' # construct drug-gene adjacency matrix for level one model
#' drugGeneNet.L1(dt=drugtarget, dg=druggene, drugAdj=drugnet_adj, geneNet =
#' gene.net, dtweight = dtweight, dgweight = dgweight, load_dir = load_dir)
#' }
#'
#' @return \code{druggeneMatrix} drug-gene adjacency matrix
#' @export
#'


drugGeneNet.L1 <- function(dt = NULL,
                           dg=NULL,
                           drugAdj,
                           geneNet,
                           dtweight,
                           dgweight,
                           load_dir){

  # x: Gene_Drug_relation.target
  # y: Gene_Drug_relation.ipa
  # remove Gene_Drug relations that are not in our background PPI network .
  if(requireNamespace(c("igraph","Matrix"))){

    # load(system.file("data", "drugtarget_lib.RData", package = "DComboNet"))
    # load(system.file("data", "druggene_lib.RData", package = "DComboNet"))

    genes <- sort(igraph::V(geneNet)$name)
    N = length(genes)
    M = nrow(drugAdj)

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

    dg <- dg[which(dg$Target %in% genes), ]
    druggeneMatrix = Matrix::Matrix(data = 0, nrow = N, ncol = M)
    rownames(druggeneMatrix) = genes
    colnames(druggeneMatrix) = colnames(drugAdj)

    for (i in 1:N){
      gene <- genes[i]
      drug <- dg$Drug[which(dg$Target == gene)]
      if (length(drug) > 0){
        for(j in 1:length(drug)){
          if(!is.na(drug[j])){
            # identify the drug position on the matrix.
            index_drug <- which(colnames(druggeneMatrix) %in% drug[j])
            # check if that index is present in the matrix.
            if(length(index_drug) == 1){
              druggeneMatrix[i,index_drug] <- dgweight
            }
          }
        }
      }
    }

    dt <- dt[which(dt$Target %in% genes), ]
    for(i in 1:N){
      gene <- genes[i]
      drug <- dt$Drug[which(dt$Target == gene)]

      if(length(drug) > 0){
        for(j in 1:length(drug)){
          if(!is.na(drug[j])){
            # identify the drug position on the matrix.
            index_drug <- which(colnames(druggeneMatrix) %in% drug[j])
            # check if that index is present in the matrix.
            if(length(index_drug) == 1){
              druggeneMatrix[i,index_drug] <- dtweight
            }
          }
        }
      }
    }
    return(druggeneMatrix)
  }
}



#  generation drug-gene bipartite matrix

#' Generate the bipartite matrix from drug-gene interaction
#' networks for level two model
#' @description The function \code{drugGeneNet.L2} is to generate the adjacency
#'   matrix from gene-gene interaction network for level two model
#' @param dt drug-target gene interaction file in \strong{data.frame} format
#'   (otherwise \code{as.data.frame} can be used to transform the format of
#'   file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol).
#' @param dg drug-gene interaction file (Optional). drug and genes that either
#'   induced by drug or may influence drug effect can be provided for construct
#'   gene network for level one model. Genes can be collected from publications
#'   or from IPA tool. Gene name should be converted to official gene names
#'   (i.e. HGNC gene symbols).
#' @param dDEG if \code{cellline} is not included in provided cancer cell line
#'   list, in \code{as.data.frame} format must be provided, otherwise
#'   \emph{ERROR} information will be returned
#' @param drugAdj drug adjacency matrix generated from drug-drug association
#'   network via \code{AdjMatrix_drugs}
#' @param geneNet gene-gene association network adjacency matrix (\bold{sparse}
#'   matrix format)
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
#' @param dtweight the weight of drug-target gene edges
#' @param dgweight the weight of drug-related gene edge
#' @param dDEGweight the weight of drug-differentially expressed genes(DEGs)
#'   edge
#' @param load_dir path to load or save modeling data files
#' @import Matrix
#' @import igraph
#' @import Hmisc
#' @examples
#' #
#' # loading files
#'
#' \dontrun{
#'   load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' drugnet_feature = read.table(paste0(load_dir,"drug_net/features.csv"),
#' sep=",",header=TRUE, stringsAsFactors = FALSE)
#' drugnet_feature = drugnet_feature[drugnet_feature$integrated_score2 >=0.2,]
#' drugtarget = read.table(paste0(load_dir,"data/drugtarget.csv"), sep =
#' ",",header = TRUE, stringsAsFactors = FALSE)
#' druggene = data.frame(drug = "Sorafenib", target  = "RAF1")
#' drugnet_adj = AdjMatrix_drugs(x = drugnet_feature, weighted = TRUE,
#' weight.type = "integrated_score")
#' drugDEG = NULL
#' # or drugDEG = data.frame(drug = "Sorafenib", target = c("FOSL1", "TUBB6"))
#' dataset = "92742"
#' cellline = "HEPG2"
#' t="6"
#' cancergene = NULL
#' # or cancergene = c("AKT3","CDH2")
#' gene.net = geneNet(dt = drugtarget, dg = druggene, dDEG = drugDEG, dataset=
#' dataset, cellline = cellline, treatment_time = t, cancer_gene = cancergene,
#' load_dir = load_dir)
#' dtweight = 2
#' dgweight = 1
#' dDEGweight = 1
#' #
#' # construct drug-gene adjacency matrix for level two model
#' drugGeneNet.L2(dt=drugtarget, dg=druggene, dDEG = drugDEG,
#' drugAdj=drugnet_adj, geneNet = gene.net, dataset = dataset, cellline =
#' cellline, treatment_time = t, dtweight = dtweight, dgweight = dgweight,
#' dDEGweight = dDEGweight, load_dir = load_dir)
#' }
#'
#' @return \code{druggeneMatrix} drug-gene adjacency matrix
#' @export
#'


drugGeneNet.L2 <- function(dt = NULL,
                           dg = NULL,
                           dDEG = NULL,
                           drugAdj,
                           geneNet,
                           dataset=NULL,
                           cellline,
                           treatment_time = NULL,
                           dtweight,
                           dgweight,
                           dDEGweight,
                           load_dir){

  if (requireNamespace(c("igraph","Matrix","Hmisc"))){

    # load(system.file("data", "drugtarget_lib.RData", package = "DComboNet"))
    # load(system.file("data", "druggene_lib.RData", package = "DComboNet"))
    # load(system.file("data", "cellline_name.RData", package = "DComboNet"))

    genes <- sort(igraph::V(geneNet)$name)
    N = length(genes)
    M = nrow(drugAdj)

    # dt_lib = read.csv(paste0(load_dir,"data/drugtarget.csv"),stringsAsFactors = F)
    names(dt_lib) = c("Drug","Target")
    if(is.null(dt)){
      dt = dt_lib
    }else{
      names(dt) = c("Drug","Target")
      dt = unique(rbind(dt_lib,dt))
    }

    # celllineDEG_lib = read.csv(paste0(load_dir,"data/celllineDEG_lib.csv"),stringsAsFactors = F)
    #names(dDEG) = names(dt)

    if(cellline %in% celllineDEG_lib[,1]){

      if(is.null(dDEG)){

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

      }else{

        dDEG_lib.path =  system.file("extdata", paste0("LINCS_data/GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), package = "DComboNet")
        if(file.exists(dDEG_lib.path)){
          dDEG_lib <- unique(read.csv(dDEG_lib.path, header = T, stringsAsFactors = F))
        }
        names(dDEG_lib) = c("Drug","Target","logFC")
        # dDEG_lib <- unique(read.csv(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), header = T, stringsAsFactors = F))
        #
        if(("logFC" %in% names(drugDEG)) != TRUE){
          dDEG$logFC = log(dDEGweight)
          names(dDEG) = c("Drug","Target","logFC")

        }
        dDEG = unique(rbind(dDEG_lib,dDEG))

        dDEG_lib$Drug = Hmisc::capitalize(dDEG_lib$Drug)
        druglist = unique(dt[1])
        dDEG = merge(dDEG,druglist, by="Drug")
        dDEG = dDEG[c("Drug","Target","logFC")]
      }
    }else{
      if(is.null(dDEG)){
        print("       ERROR: Calling Level 2 Model, Drug-DEG File Must Be Provided       ")
      }else{


        druglist = data.frame(Drug = unique(dt$Drug))
        # druglist = data.frame(Drug = colnames(drugAdj))
        if(("logFC" %in% names(dDEG)) != TRUE){
          dDEG$logFC = log(dDEGweight)
        }
        names(dDEG) = c("Drug","Target","logFC")
        dDEG$Drug =  Hmisc::capitalize(dDEG$Drug)
        dDEG = merge(dDEG,druglist, by="Drug")
        dDEG = dDEG[c("Drug","Target","logFC")]
      }
    }

    druggeneMatrix = Matrix::Matrix(data = 0, nrow = N, ncol = M)
    rownames(druggeneMatrix) = genes
    colnames(druggeneMatrix) = colnames(drugAdj)
    #
    dDEG <- dDEG[which(dDEG$Target %in% genes), ]

    for(i in 1:N){
      gene <- genes[i]
      drug <- dDEG$Drug[which(dDEG$Target == gene)]
      if (length(drug) > 0){
        for(j in 1:length(drug)){
          if(!is.na(drug[j])){
            index_drug <- which(colnames(druggeneMatrix) %in% drug[j])
            if(length(index_drug) == 1){
              druggeneMatrix[i,index_drug] <- 2^dDEG[dDEG$Drug == drug[j] & dDEG$Target == gene,]$logFC
            }
          }
        }
      }
    }

    dt <- dt[which(dt$Target %in% genes), ]
    for (i in 1:N){
      gene <- genes[i]
      drug <- dt$Drug[which(dt$Target == gene)]
      if (length(drug) > 0){
        for(j in 1:length(drug)){
          if(!is.na(drug[j])){
            index_drug <- which(colnames(druggeneMatrix) %in% drug[j])
            if(length(index_drug) == 1){
              druggeneMatrix[i,index_drug] <-  2^max(dDEG$logFC) +1
            }
          }
        }
      }
    }

    return(druggeneMatrix)
  }
}

