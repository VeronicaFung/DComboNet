
#' L1.LOCOCV

#' @description the function "L1.LOCOCV" is to predict possible drug to pair
#'   with the drug candidate without expression data
#'
#' @param load_dir path to load or save modelin data files
#' @param resultdir path to save result files
#' @param model choose model type, two models here to choose: "L1", "L2"
#' @param drugcandidate drug and its corresponding "drugbank ID" in \code{data
#'   frame}
#' @param CDK_FP Drug fingerprint matrix generated from PEDAL
#' @param pubchemFP Drug Pubchem fingerprint matrix generated from PEDAL
#' @param MACCS_FP Drug MACCS fingerprint matrix generated from PEDAL
#' @param drugnetWeight \code{TRUE/FALSE}, choose if the edge of drug network
#'   should be weighted by feature values. The abbreviationS "\code{T/F}" are
#'   also acceptable
#' @param featuretype if parameter "\code{weighted = TRUE/T}", choose the
#'   feature type to weight the edge of drug-drug network. The options are:
#'   "ATCsim","STsim","SEsim","integrated_score","integrated_score2"
#' @param drugtarget drug-target gene interaction file in \code{data.frame}
#'   format (otherwise \code{as.data.frame} can be used to transform the format
#'   of file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol)
#' @param druggene drug-gene interaction file (Optional). drug and genes that
#'   either induced by drug or may influence drug effect can be provided for
#'   construct gene network for level one model. Genes can be collected from
#'   publications or from IPA tool. Gene name should be converted to official
#'   gene names (i.e. HGNC gene symbols)
#' @param cancergene cancer_gene a list of cancer related genes. For level one
#'   model here, \code{cancer_gene} list is provided within the function, it can
#'   be \code{NULL}.
#' @param dtweight the weight of drug-target gene edges
#' @param dgweight the weight of drug-related gene edge
#' @param dgpweight numeric, weight of drug-pathway associations through drug
#'   related genes
#' @param gpweight the weight of gene-pathway edges
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
#' @param eta numeric parameter to controls the probability of restarting in the
#'   corresponding network
#' @param r numeric, global restart parameter
#' @import igraph
#' @import compiler
#' @import Matrix
#' @import reshape2
#' @return Ranked result tables of global proximity between drug seed and other
#'   drug nodes computed from drug-gene/pathway interaction network. File will
#'   be saved under the path \code{resultdir} provided by user.
#' @export
#'

L1.LOCOCV <- function(load_dir,
                      resultdir,
                      model = "L1",
                      drugcandidate = NULL,
                      CDK_FP = NULL,
                      pubchemFP = NULL,
                      MACCS_FP = NULL,
                      drugnetWeight = c("TRUE","T","FALSE","F"),
                      featuretype = c("ATCsim","STsim","SEsim","integrated_score"),
                      drugtarget,
                      druggene = NULL,
                      cancergene = NULL,
                      dtweight = 2,
                      dgweight = 1,
                      dgpweight = 1,
                      gpweight = 1,
                      x = 0.5,
                      y = 0.5,
                      z = 0,
                      A = 0.5,
                      B = 0.5,
                      eta = 1,
                      r = 0.7){

  if (requireNamespace(c("igraph","data.table","compiler","Matrix","Hmisc","reshape2","devtools"))){

    dir.create(resultdir)
    dir.create(paste0(resultdir,model,"_result/"))
    dir.create(paste0(resultdir,model,"_result/drugrank"))
    dir.create(paste0(resultdir,model,"_result/generank"))
    dir.create(paste0(resultdir,model,"_result/pathwayrank"))
    dir.create(paste0(resultdir,model,"_result/potential_net"))
    #drugseeds = readline("Which drug are you interested in? ")

    c.DrugNetFeature = compiler::cmpfun(DrugNetFeature)
    c.geneNet = compiler::cmpfun(geneNet)
    c.geneNet.Adj= compiler::cmpfun(geneNet.Adj)
    c.drugGeneNet.L1= compiler::cmpfun(drugGeneNet.L1)
    c.drugGeneNet.L2= compiler::cmpfun(drugGeneNet.L2)
    c.transitionMatrix = compiler::cmpfun(TransitionMatrix)

    #loading internal data
    # load(system.file("data", "druglistlib.RData", package = "DComboNet"))
    # load(system.file("data", "pathway_pathway_interaction.RData", package = "DComboNet"))


    #timestart<-Sys.time()
    print("----- Drug Net Features Calculating -----")

    # druglistlib = read.csv(paste0(load_dir,"data/druglist.csv"),stringsAsFactors = F)

    if(is.null(drugcandidate) |length(setdiff(drugcandidate[,1], druglistlib[,1])) == 0 ){

      if(model == "L1" | is.null(cellline)){

        feature.path = system.file("extdata", "drug_net/features.csv", package = "DComboNet")

        if(file.exists(feature.path)){

          print("       The drug network features has been found in Pre-calculated L1 model       ")
          drugnet_feature = read.table(feature.path, sep = ",", header = T, stringsAsFactors = F)

        }else{

          drugnet_feature = c.DrugNetFeature(druglist = drugcandidate,
                                             model = model,
                                             cellline = cellline,
                                             CDK_FP = CDK_FP,
                                             pubchemFP = pubchemFP,
                                             MACCS_FP = MACCS_FP,
                                             load_dir = load_dir)
        }

      }else if(model == "L2" | (is.null(cellline)==FALSE)){

        feature.path =  system.file("extdata", paste0("drug_net/features_",cellline,".csv"), package = "DComboNet")

        if( file.exists(feature.path)){

          print(paste0("       The drug network for ",cellline," features has been found in Pre-calculated L2 model       "))
          drugnet_feature = read.table(feature.path, sep = ",", header = T, stringsAsFactors = F)

        }else{

          drugnet_feature = c.DrugNetFeature(druglist = drugcandidate,
                                             model = model,
                                             cellline = cellline,
                                             CDK_FP = CDK_FP,
                                             pubchemFP = pubchemFP,
                                             MACCS_FP = MACCS_FP,
                                             load_dir = load_dir)
        }
      }

    }else{

      drugnet_feature = c.DrugNetFeature(druglist = drugcandidate,
                                         model = model,
                                         cellline = cellline,
                                         CDK_FP = CDK_FP,
                                         pubchemFP = pubchemFP,
                                         MACCS_FP = MACCS_FP,
                                         load_dir = load_dir)
    }

    drugnet_feature = drugnet_feature[drugnet_feature$integrated_score2 >=0.2,]

    #timeend<-Sys.time()
    #runningtime<-timeend-timestart
    #print(runningtime)
    # cost about 6.968375 mins
    print("----- Drug Net Features Calculation Finished -----")
    print(" ")
    print("----- Drug Net Adjacency Matrix Generating -----")
    drugnet_adj = AdjMatrix_drugs(x = drugnet_feature,
                                  weighted = drugnetWeight,
                                  weight.type = featuretype)


    if(is.null(drugDEG) & (model == "L2" | (is.null(cellline)==FALSE))){

      drugDEG_preparation(cellline = cellline,
                          dataset = dataset,
                          treatment_time =  treatment_time,
                          foldchange = foldchange_DEG,
                          pvalue = pvalue_DEG,
                          load_dir = load_dir)

      drugDEP_preparation(cellline = cellline,
                          dataset = dataset,
                          treatment_time =  treatment_time,
                          foldchange = foldchange_DEP,
                          pvalue = pvalue_DEP,
                          load_dir = load_dir)
    }

    print(" ")
    print("----- Gene Net Generating -----")
    print(" ")

    if(model == "L1" | is.null(cellline)){

      gene.net <- c.geneNet(dt= drugtarget,
                            dg = druggene,
                            dDEG = drugDEG,
                            cellline = cellline,
                            cancer_gene = cancergene,
                            load_dir = load_dir)

    }else if(model == "L2" | (is.null(cellline)==FALSE)){

      gene.net <- c.geneNet(dt= drugtarget,
                            dg = druggene,
                            dDEG = drugDEG,
                            dataset = dataset,
                            cellline = cellline,
                            treatment_time = treatment_time,
                            cancer_gene = cancergene,
                            load_dir = load_dir)

    }

    print("-----  Gene Net Adjacency Matrix Generating -----")
    print(" ")
    geneadj <- c.geneNet.Adj(GeneNetwork = gene.net)

    print("----- Drug Gene Adjacency Matrix Generating -----")


    if(model == "L1" | is.null(cellline)){

      print("------ Creating Drug Gene Adjacency Matrix for Level 1 Model ------")
      dgAdj = drugGeneNet.L1(dt = drugtarget,
                             dg = druggene,
                             drugAdj = drugnet_adj,
                             geneNet = gene.net,
                             dtweight = dtweight,
                             dgweight = dgweight,
                             load_dir = load_dir)

    }else if(model == "L2" | (is.null(cellline)==FALSE)){

      print("------ Creating Drug Gene Adjacency Matrix for Level 2 Model ------")

      dgAdj = drugGeneNet.L2(dt = drugtarget,
                             dg = druggene,
                             dDEG = drugDEG,
                             drugAdj = drugnet_adj,
                             geneNet = gene.net,
                             dataset = dataset,
                             cellline = cellline,
                             treatment_time = treatment_time,
                             dtweight = dtweight,
                             dgweight = dgweight,
                             dDEGweight = dDEGweight,
                             load_dir = load_dir)
    }


    # Pathway_pathway net -> Adj Matrix

    # WWI_all <- read.csv(paste0(load_dir,"pathway/WWI.txt"), sep = "\t", header =T, stringsAsFactors = F)

    print("----- Pathway Adjacency Matrix Generating -----")

    pathwayadj = pathwayNet.Adj(PathwayNetwork = WWI_all[c(1,2)])

    # Drug_Pathway_Matrix

    print("----- Drug Pathway Adjacency Matrix Generating -----")

    if(model == "L1" | is.null(cellline)){

      print("------ Creating Drug Pathway Adjacency Matrix for Level 1 Model ------")
      drugpathwayMatrix <- DrugPathwayMatrix.L1(dt = drugtarget,
                                                dg = druggene,
                                                drugAdj = drugnet_adj,
                                                dgpweight = dgpweight,
                                                PathwayNetwork = WWI_all,
                                                load_dir = load_dir)

    }else if(model == "L2" | (is.null(cellline)==FALSE)){

      print("------ Creating Drug Pathway Adjacency Matrix for Level 2 Model ------")

      drugpathwayMatrix <- DrugPathwayMatrix.L2(dt = drugtarget,
                                                dDEG = drugDEG,
                                                dataset=dataset,
                                                cellline = cellline,
                                                treatment_time=treatment_time,
                                                drugDEP = drugDEP,
                                                drugAdj = drugnet_adj,
                                                dDEpweight = dDEpweight,
                                                PathwayNetwork = WWI_all,
                                                load_dir = load_dir)
    }


    # genepathwayMatrix

    gpAdj = genepathwayAdj(drugAdj = drugnet_adj,
                           pathwayadj = pathwayadj,
                           geneadj = geneadj,
                           gpweight = gpweight,
                           load_dir = load_dir)



    N.gene = nrow(geneadj)
    N.drug = nrow(drugnet_adj)

    dir.create(resultdir)
    dir.create(paste0(resultdir,model,"_result/"))
    dir.create(paste0(resultdir,model,"_result/drugrank"))
    dir.create(paste0(resultdir,model,"_result/generank"))
    dir.create(paste0(resultdir,model,"_result/pathwayrank"))
    #drugseeds = readline("Which drug are you interested in? ")


    print("----- Preparaing Drug Seed -----")

    if(is.null(cellline)){

      drugpair = drugnet_feature[drugnet_feature$TAG=="P",][c(1,2)]

    }else if(!is.null(cellline)){
      drugpair = drugnet_feature[drugnet_feature$TAG=="0",][c(1,2)]
    }
    names(drugpair) = c("A","B")
    drugpair$A = capitalize(drugpair$A)
    drugpair$B = capitalize(drugpair$B)
    drugpair1 = drugpair[c(2,1)]
    names(drugpair1) = c("A","B")
    drugpair = rbind(drugpair, drugpair1)




    print("----- Preparaing Drug Seed -----")


    for(i in 1:nrow(drugpair)){

      drugseeds <- drugpair[i,1]

      if(drugseeds %in% colnames(drugnet_adj)){


        if(model == "L1" | is.null(cellline)){

          print("------ Level 1 Model ------")

          seed_score = seedscore(seeds = drugseeds,
                                 eta = eta)
          D.N1 <- which(colnames(drugnet_adj) %chin% drugseeds)
          D.N2 <- which(colnames(drugnet_adj) %chin% drugpair[i,2])

          drugnet_adj1 = drugnet_adj

          if(nrow( drugnet_feature[drugnet_feature$A == drugseeds & drugnet_feature$B == drugpair[i,2],])!=0){
            drugnet_adj1[D.N1, D.N2] = drugnet_feature[drugnet_feature$A == drugseeds & drugnet_feature$B == drugpair[i,2],]$integrated_score2
            drugnet_adj1[D.N2, D.N1] = drugnet_adj1[D.N1, D.N2]
          }else if(nrow( drugnet_feature[drugnet_feature$A == drugpair[i,2] & drugnet_feature$B == drugseeds,])!=0){
            drugnet_adj1[D.N1, D.N2] = drugnet_feature[drugnet_feature$A == drugpair[i,2] & drugnet_feature$B == drugseeds,]$integrated_score2
            drugnet_adj1[D.N2, D.N1] =drugnet_adj1[D.N1, D.N2]

          }


          dgTranMatrix <- TransitionMatrix(drugAdj = drugnet_adj1,
                                           geneAdj = geneadj,
                                           pathwayAdj = pathwayadj,
                                           druggeneAdj = dgAdj,
                                           genepathwayAdj=gpAdj,
                                           drugpathwayAdj = drugpathwayMatrix,
                                           x=x,
                                           y=y,
                                           z=z,
                                           A=A,
                                           B=B)



          N.gene = nrow(geneadj)
          N.drug = nrow(drugnet_adj)


          rwr_result = rwr(tm = dgTranMatrix,
                           r = r,
                           seeds_score = seed_score)

          drugs_rank = rank_drugs(Num.Gene = N.gene,
                                  Num.Drug = N.drug,
                                  RWR.result = rwr_result,
                                  Drug_seeds = drugseeds)

          drugs_rank$rank = 1:nrow(drugs_rank)


        }

        write.csv(drugs_rank,paste0(resultdir,model,"_result/drugrank/",drugseeds,"_",drugpair[i,2],"_rank.csv"), row.names = F, quote = F)
        print(paste0("The prediction result of drug name ",drugseeds," has been saved!"))
        #return(drugs_rank)

      }else{
        print(paste0("The drug name ",drugseeds," can not be found in the drug network! "))
      }
    }
  }
}




#' L2.LOCOCV

#' @description the function "L2.LOCOCV" is to predict possible drug to pair
#'   with the drug candidate without expression data
#'
#' @param load_dir path to load or save modelin data files
#' @param resultdir path to save result files
#' @param model choose model type, two models here to choose: "L1", "L2"
#' @param drugcandidate drug and its corresponding "drugbank ID" in \code{data
#'   frame}
#' @param CDK_FP Drug fingerprint matrix generated from PEDAL
#' @param pubchemFP Drug Pubchem fingerprint matrix generated from PEDAL
#' @param MACCS_FP Drug MACCS fingerprint matrix generated from PEDAL
#' @param drugnetWeight \code{TRUE/FALSE}, choose if the edge of drug network
#'   should be weighted by feature values. The abbreviationS "\code{T/F}" are
#'   also acceptable
#' @param featuretype if parameter "\code{weighted = TRUE/T}", choose the
#'   feature type to weight the edge of drug-drug network. The options are:
#'   "ATCsim","STsim","SEsim","integrated_score","integrated_score2"
#' @param drugtarget drug-target gene interaction file in \code{data.frame}
#'   format (otherwise \code{as.data.frame} can be used to transform the format
#'   of file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol)
#' @param druggene drug-gene interaction file (Optional). drug and genes that
#'   either induced by drug or may influence drug effect can be provided for
#'   construct gene network for level one model. Genes can be collected from
#'   publications or from IPA tool. Gene name should be converted to official
#'   gene names (i.e. HGNC gene symbols)
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
#' @param foldchange numeric, fold change of gene expression between
#'   drug-treated cancer cell line and DMSO treated cell line
#' @param pvalue numeric, p-values corresponding to the t-statistics
#' @param drugDEG if \code{cellline} is not included in provided cancer cell
#'   line list, in \code{as.data.frame} format must be provided, otherwise
#'   \emph{ERROR} information will be returned
#' @param cancergene cancer_gene a list of cancer related genes. In level two
#'   model, if cell line name is included in our provided cancer specific gene
#'   lists, \code{NULL} can take as input, otherwise \code{cancer_gene} should
#'   be prepared in \code{data.frame} format
#' @param dtweight the weight of drug-target gene edges
#' @param dgweight the weight of drug-related gene edge
#' @param dDEGweight the weight of drug-differentially expressed genes(DEGs)
#'   edge
#' @param dgpweight numeric, weight of drug-pathway associations through drug
#'   related genes
#' @param dDEpweight numeric, weight of drug-drug-differentially expressed
#'   pathway edges
#' @param gpweight the weight of gene-pathway edges
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
#' @param eta numeric parameter to controls the probability of restarting in the
#'   corresponding network
#' @param r numeric, global restart parameter
#' @import igraph
#' @import compiler
#' @import Matrix
#' @import reshape2
#' @return Ranked result tables of global proximity between drug seed and other
#'   drug nodes computed from drug-gene/pathway interaction network. File will
#'   be saved under the path \code{resultdir} provided by user.
#' @export


L2.LOCOCV <- function(load_dir,
                      resultdir,
                      model = "L2",
                      drugcandidate = NULL,
                      CDK_FP = NULL,
                      pubchemFP = NULL,
                      MACCS_FP = NULL,
                      drugnetWeight = c("TRUE","T","FALSE","F"),
                      featuretype = c("ATCsim","STsim","SEsim","integrated_score"),
                      drugtarget,
                      druggene = NULL,
                      dataset = NULL,
                      cellline = NULL,
                      treatment_time = NULL,
                      foldchange = NULL,
                      pvalue = NULL,
                      drugDEG = NULL,
                      cancergene = NULL,
                      dtweight = 2,
                      dgweight = 1,
                      dDEGweight = NULL,
                      dgpweight = 1,
                      dDEpweight = NULL,
                      gpweight=1,
                      x = 0.5,
                      y = 0.5,
                      z = 0,
                      A = 0.5,
                      B = 0.5,
                      eta = 1,
                      r = 0.7){

  if (requireNamespace(c("igraph","data.table","compiler","Matrix","Hmisc","reshape2","devtools"))){

    dir.create(resultdir)
    dir.create(paste0(resultdir,model,"_result/"))
    dir.create(paste0(resultdir,model,"_result/drugrank"))
    dir.create(paste0(resultdir,model,"_result/generank"))
    dir.create(paste0(resultdir,model,"_result/pathwayrank"))
    dir.create(paste0(resultdir,model,"_result/potential_net"))
    #drugseeds = readline("Which drug are you interested in? ")

    #Loading functions
    devtools::load_all()
    # source(paste0(load_dir,"/DComboNet/R/DEG_DEpathway_preparation.R"))
    # source(paste0(load_dir,"/DComboNet/R/DrugNetFeatures.R"))
    # source(paste0(load_dir,"/DComboNet/R/DrugNet.R"))
    # source(paste0(load_dir,"/DComboNet/R/DrugGeneNet.R"))
    # source(paste0(load_dir,"/DComboNet/R/GeneNet.R"))
    # source(paste0(load_dir,"/DComboNet/R/DrugPathwayNet.R"))
    # source(paste0(load_dir,"/DComboNet/R/PathwayNet.R"))
    # source(paste0(load_dir,"/DComboNet/R/GenePathwayNet.R"))
    # source(paste0(load_dir,"/DComboNet/R/DrugGeneHeteroNet.R"))
    # source(paste0(load_dir,"/DComboNet/R/RWR_fun.R"))
    # source(paste0(load_dir,"/DComboNet/R/SeedsPreparation.R"))
    # source(paste0(load_dir,"/DComboNet/R/Results_rank.R"))


    c.DrugNetFeature = compiler::cmpfun(DrugNetFeature)
    c.geneNet = compiler::cmpfun(geneNet)
    c.geneNet.Adj= compiler::cmpfun(geneNet.Adj)
    c.drugGeneNet.L1= compiler::cmpfun(drugGeneNet.L1)
    c.drugGeneNet.L2= compiler::cmpfun(drugGeneNet.L2)
    c.transitionMatrix = compiler::cmpfun(TransitionMatrix)

    #loading internal data
    # load(system.file("data", "druglistlib.RData", package = "DComboNet"))
    # load(system.file("data", "pathway_pathway_interaction.RData", package = "DComboNet"))


    #timestart<-Sys.time()
    print("----- Drug Net Features Calculating -----")

    # druglistlib = read.csv(paste0(load_dir,"data/druglist.csv"),stringsAsFactors = F)

    if(is.null(drugcandidate) |length(setdiff(drugcandidate[,1], druglistlib[,1])) == 0 ){

      if(model == "L1" | is.null(cellline)){

        feature.path = system.file("extdata", "drug_net/features.csv", package = "DComboNet")

        if(file.exists(feature.path)){

          print("       The drug network features has been found in Pre-calculated L1 model       ")
          drugnet_feature = read.table(feature.path, sep = ",", header = T, stringsAsFactors = F)

        }else{

          drugnet_feature = c.DrugNetFeature(druglist = drugcandidate,
                                             model = model,
                                             cellline = cellline,
                                             CDK_FP = CDK_FP,
                                             pubchemFP = pubchemFP,
                                             MACCS_FP = MACCS_FP,
                                             load_dir = load_dir)
        }

      }else if(model == "L2" | (is.null(cellline)==FALSE)){

        feature.path =  system.file("extdata", paste0("drug_net/features_",cellline,".csv"), package = "DComboNet")

        if( file.exists(feature.path)){

          print(paste0("       The drug network for ",cellline," features has been found in Pre-calculated L2 model       "))
          drugnet_feature = read.table(feature.path, sep = ",", header = T, stringsAsFactors = F)

        }else{

          drugnet_feature = c.DrugNetFeature(druglist = drugcandidate,
                                             model = model,
                                             cellline = cellline,
                                             CDK_FP = CDK_FP,
                                             pubchemFP = pubchemFP,
                                             MACCS_FP = MACCS_FP,
                                             load_dir = load_dir)
        }
      }

    }else{

      drugnet_feature = c.DrugNetFeature(druglist = drugcandidate,
                                         model = model,
                                         cellline = cellline,
                                         CDK_FP = CDK_FP,
                                         pubchemFP = pubchemFP,
                                         MACCS_FP = MACCS_FP,
                                         load_dir = load_dir)
    }

    drugnet_feature = drugnet_feature[drugnet_feature$integrated_score2 >=0.2,]

    #timeend<-Sys.time()
    #runningtime<-timeend-timestart
    #print(runningtime)
    # cost about 6.968375 mins
    print("----- Drug Net Features Calculation Finished -----")
    print(" ")
    print("----- Drug Net Adjacency Matrix Generating -----")
    drugnet_adj = AdjMatrix_drugs(x = drugnet_feature,
                                  weighted = drugnetWeight,
                                  weight.type = featuretype)


    if(is.null(drugDEG) & (model == "L2" | (is.null(cellline)==FALSE))){

      drugDEG_preparation(cellline = cellline,
                          dataset = dataset,
                          treatment_time =  treatment_time,
                          foldchange = foldchange_DEG,
                          pvalue = pvalue_DEG,
                          load_dir = load_dir)

      drugDEP_preparation(cellline = cellline,
                          dataset = dataset,
                          treatment_time =  treatment_time,
                          foldchange = foldchange_DEP,
                          pvalue = pvalue_DEP,
                          load_dir = load_dir)
    }

    print(" ")
    print("----- Gene Net Generating -----")
    print(" ")

    if(model == "L1" | is.null(cellline)){

      gene.net <- c.geneNet(dt= drugtarget,
                            dg = druggene,
                            dDEG = drugDEG,
                            cellline = cellline,
                            cancer_gene = cancergene,
                            load_dir = load_dir)

    }else if(model == "L2" | (is.null(cellline)==FALSE)){

      gene.net <- c.geneNet(dt= drugtarget,
                            dg = druggene,
                            dDEG = drugDEG,
                            dataset = dataset,
                            cellline = cellline,
                            treatment_time = treatment_time,
                            cancer_gene = cancergene,
                            load_dir = load_dir)

    }

    print("-----  Gene Net Adjacency Matrix Generating -----")
    print(" ")
    geneadj <- c.geneNet.Adj(GeneNetwork = gene.net)

    print("----- Drug Gene Adjacency Matrix Generating -----")


    if(model == "L1" | is.null(cellline)){

      print("------ Creating Drug Gene Adjacency Matrix for Level 1 Model ------")
      dgAdj = drugGeneNet.L1(dt = drugtarget,
                             dg = druggene,
                             drugAdj = drugnet_adj,
                             geneNet = gene.net,
                             dtweight = dtweight,
                             dgweight = dgweight,
                             load_dir = load_dir)

    }else if(model == "L2" | (is.null(cellline)==FALSE)){

      print("------ Creating Drug Gene Adjacency Matrix for Level 2 Model ------")

      dgAdj = drugGeneNet.L2(dt = drugtarget,
                             dg = druggene,
                             dDEG = drugDEG,
                             drugAdj = drugnet_adj,
                             geneNet = gene.net,
                             dataset = dataset,
                             cellline = cellline,
                             treatment_time = treatment_time,
                             dtweight = dtweight,
                             dgweight = dgweight,
                             dDEGweight = dDEGweight,
                             load_dir = load_dir)
    }


    # Pathway_pathway net -> Adj Matrix

    # WWI_all <- read.csv(paste0(load_dir,"pathway/WWI.txt"), sep = "\t", header =T, stringsAsFactors = F)

    print("----- Pathway Adjacency Matrix Generating -----")

    pathwayadj = pathwayNet.Adj(PathwayNetwork = WWI_all[c(1,2)])

    # Drug_Pathway_Matrix

    print("----- Drug Pathway Adjacency Matrix Generating -----")

    if(model == "L1" | is.null(cellline)){

      print("------ Creating Drug Pathway Adjacency Matrix for Level 1 Model ------")
      drugpathwayMatrix <- DrugPathwayMatrix.L1(dt = drugtarget,
                                                dg = druggene,
                                                drugAdj = drugnet_adj,
                                                dgpweight = dgpweight,
                                                PathwayNetwork = WWI_all,
                                                load_dir = load_dir)

    }else if(model == "L2" | (is.null(cellline)==FALSE)){

      print("------ Creating Drug Pathway Adjacency Matrix for Level 2 Model ------")

      drugpathwayMatrix <- DrugPathwayMatrix.L2(dt = drugtarget,
                                                dDEG = drugDEG,
                                                dataset=dataset,
                                                cellline = cellline,
                                                treatment_time=treatment_time,
                                                drugDEP = drugDEP,
                                                drugAdj = drugnet_adj,
                                                dDEpweight = dDEpweight,
                                                PathwayNetwork = WWI_all,
                                                load_dir = load_dir)
    }


    # genepathwayMatrix

    gpAdj = genepathwayAdj(drugAdj = drugnet_adj,
                           pathwayadj = pathwayadj,
                           geneadj = geneadj,
                           gpweight = gpweight,
                           load_dir = load_dir)



    N.gene = nrow(geneadj)
    N.drug = nrow(drugnet_adj)

    dir.create(resultdir)
    dir.create(paste0(resultdir,model,"_result/"))
    dir.create(paste0(resultdir,model,"_result/drugrank"))
    dir.create(paste0(resultdir,model,"_result/generank"))
    dir.create(paste0(resultdir,model,"_result/pathwayrank"))
    #drugseeds = readline("Which drug are you interested in? ")


    if(!is.null(cellline)){
      drugpair = drugnet_feature[drugnet_feature$TAG=="0",][c(1,2)]
    }
    names(drugpair) = c("A","B")
    drugpair$A = capitalize(drugpair$A)
    drugpair$B = capitalize(drugpair$B)
    drugpair1 = drugpair[c(2,1)]
    names(drugpair1) = c("A","B")
    drugpair = rbind(drugpair, drugpair1)



    print("----- Preparaing Drug Seed -----")


    for(i in 1:nrow(drugpair)){

      drugseeds <- drugpair[i,1]

      if(drugseeds %in% colnames(drugnet_adj)){

        if(model == "L2" | (is.null(cellline)==FALSE)){

          print("------ Level 2 Model ------")

          seed_score = seedscore(seeds = drugseeds,
                                 eta = eta)
          D.N1 <- which(colnames(drugnet_adj) %chin% drugseeds)
          D.N2 <- which(colnames(drugnet_adj) %chin% drugpair[i,2])

          drugnet_adj1 = drugnet_adj

          if(nrow( drugnet_feature[drugnet_feature$A == drugseeds & drugnet_feature$B == drugpair[i,2],])!=0){
            drugnet_adj1[D.N1, D.N2] = drugnet_feature[drugnet_feature$A == drugseeds & drugnet_feature$B == drugpair[i,2],]$integrated_score2
            drugnet_adj1[D.N2, D.N1] = drugnet_adj1[D.N1, D.N2]
          }else if(nrow( drugnet_feature[drugnet_feature$A == drugpair[i,2] & drugnet_feature$B == drugseeds,])!=0){
            drugnet_adj1[D.N1, D.N2] = drugnet_feature[drugnet_feature$A == drugpair[i,2] & drugnet_feature$B == drugseeds,]$integrated_score2
            drugnet_adj1[D.N2, D.N1] = drugnet_adj1[D.N1, D.N2]

          }


          dgTranMatrix <- TransitionMatrix(drugAdj = drugnet_adj1,
                                           geneAdj = geneadj,
                                           pathwayAdj = pathwayadj,
                                           druggeneAdj = dgAdj,
                                           genepathwayAdj=gpAdj,
                                           drugpathwayAdj = drugpathwayMatrix,
                                           x=x,
                                           y=y,
                                           z=z,
                                           A=A,
                                           B=B)

          N.gene = nrow(geneadj)
          N.drug = nrow(drugnet_adj)

          rwr_result = rwr(tm = dgTranMatrix,
                           r = r,
                           seeds_score = seed_score)

          drugs_rank = rank_drugs(Num.Gene = N.gene,
                                  Num.Drug = N.drug,
                                  RWR.result = rwr_result,
                                  Drug_seeds = drugseeds)

          drugs_rank$rank = 1:nrow(drugs_rank)


        }


        write.csv(drugs_rank,paste0(resultdir,model,"_result/drugrank/",drugseeds,"_",drugpair[i,2],"_rank.csv"), row.names = F, quote = F)

        print(paste0("The prediction result of drug name ",drugseeds," has been saved!"))

      }else{
        print(paste0("The drug name ",drugseeds," can not be found in the drug network! "))
      }
    }
  }
}
