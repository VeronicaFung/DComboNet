#' Generate cancer gene list
#' @description The function \code{cancerGene} is to generate cancer sample
#'   specific expressed gene list from CCLE database. Cancer sample specific
#'   expressed genes are obtained by comparing gene expression of certain cancer
#'   sample with the rest, fold-change above 1.5 times is considered as specific
#'   expressed genes.
#' @param cellline cancer cell lines name, if requested cell line name is
#'   included in the default cell line list, cancer sample specific expressed
#'   gene list will be computed and saved in given path; otherwise
#'   \emph{ERROR} report and \emph{SUGGESTION} will be returned
#' @param load_dir defined the path to load and save model data files
#' @return \code{cell_genelist} cancer cell line specific expressed gene list,
#'   which will be returned and saved in provided \code{load_dir} under the
#'   folder named "\emph{gene_network}" in \emph{.txt} format, file will be
#'   named as " "\code{cellline}"_genelist.txt" name you provide
#' @export
#'
#' @examples
#' \dontrun{
#' cellline = "HEPG2"
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#'
#' cancerGene(cellline = cellline, load_dir = load_dir)
#' }

cancerGene = function(cellline,
                      load_dir){
  # ccleexpr = read.table(paste0(load_dir,"data/celllineExpr/CCLEexprMatrix.txt"),
  # sep = "\t",header=T, stringsAsFactors = F)
  load(paste0(load_dir, "celllineExpr/CCLEexprMatrix.RData"))

  if(cellline %in% colnames(ccleexpr)){
    cell <- which(colnames(ccleexpr) == cellline)
    cellExp <- ccleexpr[c(2,cell)]
    otherExp <-  ccleexpr[,-c(1,cell)]
    ExpMean <- data.frame(ExpMean = apply(otherExp[2:ncol(otherExp)], 1, mean))
    ExpMean <- cbind(otherExp[1],ExpMean)
    cellFC <- abs(cellExp[2] - ExpMean[2])
    cell <- cbind(otherExp[1],cellFC)
    cell_genelist <- cell[cell[2] > 1.5,][1]
    write.table(cell_genelist, paste0(load_dir,"cellline_genes/",cellline,"_genelist.txt", sep = ""), sep = "\t", col.names = T, row.names = F)
  }else{
    print("ERROR: The cell line name you put in does not match the name in CCLE or it did not contain in CCLE.")
    print("SUGGESTION: Please check the name with CCLE or upload your own cancer gene list. ")
  }
  return(cell_genelist)
}



#' Generate gene-gene interaction networks
#' @description \code{geneNet} is to generate the gene-gene interaction
#'   networks
#' @param dt drug-target gene interaction file in \code{data.frame} format
#'   (otherwise \code{as.data.frame} can be used to transform the format of
#'   file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol)
#' @param dg drug-gene interaction file (Optional). drug and genes that either
#'   induced by drug or may influence drug effect can be provided for construct
#'   gene network for level one model. Genes can be collected from publications
#'   or from IPA tool. Gene name should be converted to official gene names
#'   (i.e. HGNC gene symbols)
#' @param dDEG if \code{cellline} is not included in provided cancer cell line
#'   list, in \code{as.data.frame} format must be provided, otherwise
#'   \emph{ERROR} information will be returned
#' @param dataset GEO ID of LINCS datasets, two datasets are provided here:
#'   "70138" and "92742"
#' @param cellline cancer cell lines name. If level two model is called
#'   (\code{model = "L2"} in \code{\link{DComboNet}} or
#'   \code{\link{L2.LOCOCV}}). If requested cell line is included in the default
#'   cancer cell line list, \code{dDEG} do not need to provided unless
#'   particular drug induced differentially expressed gene list you want to
#'   include in the gene network, then provided \code{dDEG} in
#'   \code{as.data.frame} format. If the cancer cell line is not in the current
#'   cell line list, you should provide \code{cellline} and \code{dDEG} for
#'   construct specific gene network, otherwise \emph{ERROR} information will be
#'   returned
#' @param treatment_time drug treatment time provided in LINCS database. For
#'   example, MCF7 cell line (\code{cellline = "MCF7"}) provided two time points
#'   in LINCS GSE\code{70138} dataset - 6 hour and 12 hour - \code{6} or
#'   \code{12} should be taken as input here, as numeric or character
#' @param cancer_gene a list of cancer related genes. If Level one model is
#'   called (\code{model = "L1"} in \code{\link{DComboNet}} or
#'   \code{\link{L1.LOCOCV}}), \code{cancer_gene} list is provided within the
#'   function, it can be \code{NULL}. In level two model, if cell line name is
#'   included in our provided cancer specific gene lists, \code{NULL} can take
#'   as input, otherwise \code{cancer_gene} should be prepared in
#'   \code{data.frame} format
#' @param load_dir path to load and save model data file(s)
#' @import igraph
#' @import dnet
#' @return \code{gene_network} \code{igraph} network file
#' @export
#' @examples
#' # 1) gene network built according to level one model
#' \dontrun{
#' drugtarget = data.frame(drug = "Sorafenib", target = "BRAF")
#' druggene = data.frame(drug = "Sorafenib", target = "RAF1")
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' geneNet(dt = drugtarget, dg = druggene, load_dir = load_dir)
#' }
#'
#' # 2) gene network built according to level two model
#' \dontrun{
#' drugtarget = data.frame(drug = "Sorafenib", target = "BRAF")
#' druggene = data.frame(drug = "Sorafenib", target = "RAF1")
#' drugDEG = NULL
#' # or drugDEG = data.frame(drug = "Sorafenib", target = c("FOSL1", "TUBB6", "RRP8"))
#' dataset = "92742"
#' cellline = "HEPG2"
#' t="6"
#' cancergene = NULL
#' # or cancergene = c("AKT3","CDH2")
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' geneNet(dt = drugtarget, dg = druggene, dDEG = drugDEG, dataset= dataset,
#' cellline = cellline, treatment_time = t, cancer_gene = cancergene, load_dir =
#' load_dir)
#' }
#'

geneNet <- function(dt = NULL,
                    dg = NULL,
                    dDEG = NULL,
                    dataset = NULL,
                    cellline = NULL,
                    treatment_time = NULL,
                    cancer_gene = NULL,
                    load_dir){

  if (requireNamespace(c("igraph","dnet"))){

    # load(system.file("data", "drugtarget_lib.RData", package = "DComboNet"))
    # load(system.file("data", "druggene_lib.RData", package = "DComboNet"))
    # load(system.file("data", "cancer_gene_L1.RData", package = "DComboNet"))
    # load(system.file("data", "inBio_PPI_net.RData", package = "DComboNet"))
    # load(system.file("data", "cellline_name.RData", package = "DComboNet"))

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

    # PPI_table <- data.table::fread("G:/lab/Projects/p1_DComboNet/DComboNet_improve/Rpackage/DComboNet-master/inst/extdata/PPI_inbio_score.csv",sep = ',', header = T, stringsAsFactors = F)
    # PPI_tablePPI_table <- data.frame(PPI_table)


    PPI_table = PPI_table[PPI_table$confidence.score >= 0.15,]
    PPI_Network1 <- igraph::graph.data.frame(PPI_table[c(1,2)],directed=FALSE)
    PPI_Network1 <- igraph::simplify(PPI_Network1, remove.multiple = TRUE, remove.loops = TRUE)

    if(is.null(cellline)){

      print("       No Cancer cell lines information received       ")
      print("       Creating Gene Network for Level 1 Model       ")

      # cancer_gene = read.csv(paste0(load_dir,"data/network/protein-coding cancer genes name.csv"),header = F, stringsAsFactors = F)
      cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,dg$Target)))

    }else{

      print("       Cancer cell lines information has been found       ")
      print("       Creating Gene Network for Level 2 Model       ")

      # cellline_lib = read.csv(paste0(load_dir,"data/celllines_lib.csv"),stringsAsFactors = F)
      # celllineDEG_lib = read.csv(paste0(load_dir,"data/celllineDEG_lib.csv"),stringsAsFactors = F)

      if(!(cellline %in% celllineDEG_lib[,1]) && is.null(dDEG)){

        print("       The Differentially Expressed Genes After Drugs Treatment Is Neither Inputed Nor In Our Library       ")

        if((cellline %in% cellline_lib[,1]) && is.null(cancer_gene)){

          print("       You Did Not Provide Cancer Specific Genes List, But It has been found In Our Library        ")
          print("       Creating Basic Cancer Specific Gene Network (Without Differentially Expressed Genes After Drugs Treatment)       ")

          file.path =  system.file("extdata", paste0("cellline_genes/",cellline,"_genelist.csv"), package = "DComboNet")
          if( file.exists(file.path)){
            cancer_gene = read.csv(file.path,header = F, stringsAsFactors = F)
          }else{
            print("       ERROR: cancer cell line specific expressed gene list should be provided for level-two model!       ")
          }

          cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,dg$Target)))

        }else if(nrow(cancer_gene)!=0){
          print("       Cancer Specific Genes Has Been Privided        ")
          print("       Creating Customized Basic Cancer Specific Gene Network (Without Differentially Expressed Genes After Drugs Treatment)       ")

          cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,dg$Target)))

        }

      }else if(cellline %in% celllineDEG_lib[,1]){

        print("       The Differentially Expressed Genes After Drugs Treatment Can Be Found In Our Library       ")


        if(is.null(dDEG)){
          print("       The Differentially Expressed Genes After Drugs Treatment Has Not Been Provided       ")

          if((cellline %in% cellline_lib[,1]) && is.null(cancer_gene)){
            print("       You Did Not Provide Cancer Specific Genes List, But It has been found In Our Library        ")
            print("       Creating Cancer Specific Gene Network       ")


            file.path =  system.file("extdata", paste0("LINCS_data/GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), package = "DComboNet")
            if( file.exists(file.path)){
              dDEG_lib <- unique(read.csv(file.path, header = T, stringsAsFactors = F))[c(1,2)]
            }else{
              print("       ERROR: drug induced differentially expressed gene list can not be found for level-two model!       ")
            }

            # dDEG_lib <- unique(read.csv(paste0(load_dir,"data/drug_DEG/",cellline,"_L5DEG3.29_2.txt"), sep = "\t", header = T, stringsAsFactors = F)[c(1,2)])
            names(dDEG_lib) = c("Drug","Target")
            dDEG_lib$Drug = Hmisc::capitalize(dDEG_lib$Drug)

            druglist = unique(dt[1])
            dDEG = merge(dDEG_lib,druglist, by="Drug")
            dDEG = dDEG[c("Drug","Target")]

            # cancer_gene = read.csv(paste0(load_dir,"cellline_genes/",cellline,"_genelist.csv", sep = ""),header = F, stringsAsFactors = F)
            file.path =  system.file("extdata", paste0("cellline_genes/",cellline,"_genelist.csv"), package = "DComboNet")
            if( file.exists(file.path)){
              cancer_gene = read.csv(file.path,header = F, stringsAsFactors = F)
            }else{
              print("       ERROR: cancer cell line specific expressed gene list should be provided for level-two model!       ")
            }
            # cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,dg$Target)))
            cancer_gene = union(cancer_gene[,1], unique(union(dt$Target,  union(dg$Target,dDEG$Target))))
            # cancer_gene = union(cancer_gene[,1], unique(union(dt$Target,  dDEG$Target)))

          }else if(!(cellline %in% cellline_lib[,1]) && nrow(cancer_gene)!=0){
            print("       Cancer Specific Genes Has Been Privided        ")
            print("       Creating Customized Cancer Specific Gene Network       ")

            file.path =  system.file("extdata", paste0("LINCS_data/GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), package = "DComboNet")
            if( file.exists(file.path)){
              dDEG_lib <- unique(read.csv(file.path, header = T, stringsAsFactors = F))[c(1,2)]
            }else{
              print("       ERROR: drug induced differentially expressed gene list can not be found for level-two model!       ")
            }
            # dDEG_lib <- unique(read.csv(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), header = T, stringsAsFactors = F))
            names(dDEG_lib) = c("Drug","Target")
            dDEG_lib$Drug = Hmisc::capitalize(dDEG_lib$Drug)
            druglist = unique(dt[1])
            dDEG = merge(dDEG_lib,druglist, by="Drug")
            dDEG = dDEG[c("Drug","Target")]
            #gene_network
            cancer_gene = union(cancer_gene[,1],  unique(union(dt$Target, union(dg$Target,dDEG$Target))))
            # cancer_gene = union(cancer_gene[,1],  unique(union(dt$Target, dDEG$Target)))
          }
        }else if(nrow(dDEG)!=0){

          print("       The Differentially Expressed Genes After Drugs Treatment Has Been Provided       ")
          print("       The Differentially Expressed Genes After Drugs Treatment Can Also Be Found In Our Library       ")
          print("       Creating Customized Cancer Specific Gene Network       ")

          if((cellline %in% cellline_lib[,1]) && is.null(cancer_gene)){
            print("       You Did Not Provide Cancer Specific Genes List, But It has been found In Our Library        ")
            print("       Creating Cancer Specific Gene Network       ")

            # dDEG_lib <- unique(read.csv(paste0(load_dir,"data/drug_DEG/",cellline,"_L5DEG3.29_2.txt"), sep = "\t", header = T, stringsAsFactors = F)[c(1,2)])
            # dDEG_lib <- unique(read.csv(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), header = T, stringsAsFactors = F)[c(1,2)])
            file.path =  system.file("extdata", paste0("LINCS_data/GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), package = "DComboNet")
            # dDEG_lib <- unique(read.csv(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), header = T, stringsAsFactors = F))
            if( file.exists(file.path)){
              dDEG_lib <- unique(read.csv(file.path, header = T, stringsAsFactors = F))[c(1,2)]
            }else{
              print("       ERROR: drug induced differentially expressed gene list can not be found for level-two model!       ")
            }

            names(dDEG_lib) = c("Drug","Target")
            dDEG_lib$Drug = Hmisc::capitalize(dDEG_lib$Drug)
            dDEG = unique(rbind(dDEG_lib,dDEG))
            druglist = unique(dt[1])
            dDEG = merge(dDEG,druglist, by="Drug")
            dDEG = dDEG[c("Drug","Target")]

            file.path =  system.file("extdata", paste0("cellline_genes/",cellline,"_genelist.csv"), package = "DComboNet")
            if( file.exists(file.path)){
              cancer_gene = read.csv(file.path,header = F, stringsAsFactors = F)
            }else{
              print("       ERROR: cancer cell line specific expressed gene list should be provided for level-two model!       ")
            }
            # cancer_gene = read.csv(paste0(load_dir,"cellline_genes/",cellline,"_genelist.csv", sep = ""),header = F, stringsAsFactors = F)
            cancer_gene = union(cancer_gene[,1],unique(union(dt$Target, union(dg$Target,dDEG$Target))))
            # cancer_gene = union(cancer_gene[,1],unique(union(dt$Target, dDEG$Target)))

          }else if(!(cellline %in% cellline_lib[,1]) && nrow(cancer_gene)!=0){
            print("       Cancer Specific Genes Has Been Privided        ")
            print("       Creating Customized Cancer Specific Gene Network       ")

            # dDEG_lib <- unique(read.csv(paste0(load_dir,"data/drug_DEG/",cellline,"_L5DEG3.29_2.txt"), sep = "\t", header = T, stringsAsFactors = F)[c(1,2)])
            file.path =  system.file("extdata", paste0("LINCS_data/GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), package = "DComboNet")
            if( file.exists(file.path)){
              dDEG_lib <- unique(read.csv(file.path, header = T, stringsAsFactors = F))[c(1,2)]
            }else{
              print("       ERROR: drug induced differentially expressed gene list can not be found for level-two model!       ")
            }
            # dDEG_lib <- unique(read.csv(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"), header = T, stringsAsFactors = F)[c(1,2)])
            names(dDEG_lib) = c("Drug","Target")
            dDEG_lib$Drug = Hmisc::capitalize(dDEG_lib$Drug)
            #names(dDEG) = names(dt)
            dDEG = unique(rbind(dDEG_lib,dDEG))
            druglist = unique(dt[1])
            dDEG = merge(dDEG,druglist, by="Drug")
            dDEG = dDEG[c("Drug","Target")]
            #gene_network
            cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,union(dg$Target,dDEG$Target))))
            # cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,dDEG$Target)))
          }
        }
      }else if(!(cellline %in% celllineDEG_lib[,1])){
        if(is.null(dDEG)){
          print("       ERROR: Drug-DEG File Must Be Provided, If Cancer Information Has Been Provided       ")
        }else{
          if((cellline %in% cellline_lib[,1]) && is.null(cancer_gene)){
            print("       You Did Not Provide Cancer Specific Genes List, But It has been found In Our Library        ")
            print("       Creating Cancer Specific Gene Network       ")

            names(dDEG) = c("Drug","Target")
            druglist = unique(dt[1])
            dDEG = merge(dDEG,druglist, by="Drug")
            dDEG = dDEG[c("Drug","Target")]

            file.path =  system.file("extdata", paste0("cellline_genes/",cellline,"_genelist.csv"), package = "DComboNet")
            if( file.exists(file.path)){
              cancer_gene = read.csv(file.path,header = F, stringsAsFactors = F)
            }else{
              print("       ERROR: cancer cell line specific expressed gene list should be provided for level-two model!       ")
            }
            # cancer_gene = read.csv(paste0(load_dir,"cellline_genes/",cellline,"_genelist.csv", sep = ""),header = F, stringsAsFactors = F)
            cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,union(dg$Target,dDEG$Target))))
            # cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,dDEG$Target)))

          }else if(!(cellline %in% cellline_lib[,1]) && nrow(cancer_gene)!=0){
            print("       Cancer Specific Genes Has Been Privided        ")
            print("       Creating Customized Cancer Specific Gene Network       ")

            names(dDEG) = c("Drug","Target")
            druglist = unique(dt[1])
            dDEG = merge(dDEG,druglist, by="Drug")
            dDEG = dDEG[c("Drug","Target")]
            #gene_network
            cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,union(dg$Target,dDEG$Target))))
            # cancer_gene = union(cancer_gene[,1],unique(union(dt$Target,dDEG$Target)))
          }
        }
      }
    }

    gene_network = dnet::dNetInduce(PPI_Network1, nodes_query=cancer_gene, knn=0, remove.loops=T, largest.comp=T)
    gene_network <- igraph::simplify(gene_network, remove.multiple = TRUE, remove.loops = TRUE)
    #gene_network = gene_network %u% PPI_Network
    igraph::E(gene_network)$weight = 1
    names(gene_network) <- "gene_network"
    #rm(PPI_table,PPI_Network1,dt,dg,cancer_gene);gc()
    return(gene_network)
  }
}





#' Generate the adjacency matrix from gene-gene interaction networks
#' @description The function \code{geneNet.Adj} is used to generate adjacency
#'   matrix from gene-gene interaction network built from \code{\link{geneNet}}
#' @param GeneNetwork \code{igraph} network file generated from
#'   \code{\link{geneNet}}
#' @import Matrix
#' @import igraph
#' @return \code{adjMatrix} gene-gene association network adjacency matrix
#'   (sparse matrix format)
#' @export
#' @seealso \code{\link{geneNet}}
#' @examples
#' # require(igraph)
#' g <- data.frame(V1 = c("g1", "g2", "g3", "g4"), V2 = c("g5", "g2", "g6",
#' "g1"))
#' g.net <- igraph::graph.data.frame(g)
#' geneNet.Adj(g.net)

geneNet.Adj <- function(GeneNetwork){
  if (requireNamespace(c("Matrix","igraph"))){
    N = length(igraph::V(GeneNetwork)$name)
    print(paste0("       The Gene Network has ",N," Nodes and ", length(igraph::E((GeneNetwork)))," Edges      "))
    adjMatrix <- Matrix::Matrix(0,ncol=N,nrow=N,sparse = TRUE)
    adjMatrix <- igraph::as_adjacency_matrix(GeneNetwork,sparse = TRUE)
    adjMatrix <- adjMatrix[order(rownames(adjMatrix)),order(colnames(adjMatrix))]
    return(adjMatrix)
  }
}


