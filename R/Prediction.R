#' Prediction drug combinations
#'
#' @description the function \code{DComboNet} is to predict possible drug to
#'   pair with the drug candidate without expression data
#' @param load_dir path to load or save modelin data files
#' @param resultdir path to save result files
#' @param model choose model type, two models here to choose: "L1", "L2"
#' @param drugcandidate drug and its corresponding "drugbank ID" in
#'   \code{data.frame}
#' @param manual_input logical (\code{TRUE(T)/FALSE(F)}), if user wants to
#'   manually input drug(s)
#' @param CDK_FP Drug fingerprint matrix generated from PEDAL
#' @param pubchemFP Drug Pubchem fingerprint matrix generated from PEDAL
#' @param MACCS_FP Drug MACCS fingerprint matrix generated from PEDAL
#' @param drugnetWeight logical (\code{TRUE(T)/FALSE(F)}), choose if the edge of
#'   drug network should be weighted by feature values. The abbreviationS
#'   "\code{T/F}" are also acceptable
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
#' @param foldchange_DEG numeric, fold change of gene expression between
#'   drug-treated cancer cell line and DMSO treated cell line
#' @param pvalue_DEG numeric, p-values corresponding to the t-statistics
#' @param foldchange_DEP numeric, fold change of enriched pathways via GSVA
#'   between drug-treated cancer cell line and DMSO treated cell line
#' @param pvalue_DEP numeric, p-values corresponding to the t-statistics
#' @param drugDEG if \code{cellline} is not included in provided cancer cell
#'   line list, in \code{as.data.frame} format must be provided, otherwise
#'   \emph{ERROR} information will be returned
#' @param drugDEP if \code{cellline} is not included in provided cancer cell
#'   line list, in \code{as.data.frame} format must be provided, otherwise
#'   \emph{ERROR} information will be returned
#' @param cancergene cancer_gene a list of cancer related genes. If Level one
#'   model is called (\code{model = "L1"} in \code{\link{DComboNet}} or
#'   \code{\link{L1.LOCOCV}}), \code{cancer_gene} list is provided within the
#'   function, it can be \code{NULL}. In level two model, if cell line name is
#'   included in our provided cancer specific gene lists, \code{NULL} can take
#'   as input, otherwise \code{cancer_gene} should be prepared in
#'   \code{data.frame} format
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
#' @importFrom utils read.table read.csv write.csv write.table
#' @return Ranked result tables. Each drug seed will create three files that
#'   contains global proximity between drug seed and other drug/gene/pathway
#'   nodes computed from drug-gene/pathway interaction network. File will be
#'   saved under the path \code{resultdir} provided by user.
#' @export
#'
#' @seealso \code{\link{L1.LOCOCV}}, \code{\link{L2.LOCOCV}}
#'
#' @examples
#'
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' resultdir = "G:/lab/DCcomboNet/Rpackage/tryout_result/"
#' drugcandidate = data.frame(drug = "Sorafenib", drug = "DB00398")
#' manual_input = FALSE
#'
#' # 1) runing level one model (L1)
#' \dontrun{
#' DComboNet(load_dir = load_dir,
#'           resultdir = resultdir,
#'           model = "L1",
#'           drugcandidate = drugcandidate,
#'           drugnetWeight = TRUE,
#'           featuretype = "integrated_score")
#' }
#'
#' # Runing level one model with extra data
#' dt_table = read.table(paste0(load_dir,"data/drugtarget.csv"), sep =
#' ",",header = TRUE, stringsAsFactors = FALSE)
#' dg_table = data.frame(drug = "Sorafenib", target  = "RAF1")
#'
#' \dontrun{
#' DComboNet(load_dir = load_dir,
#'           resultdir = resultdir,
#'           model = "L1",
#'           manual_input = FALSE,
#'           drugcandidate = drugcandidate,
#'           drugnetWeight = TRUE,
#'           featuretype = "integrated_score",
#'           drugtarget = dt_table,
#'           druggene = dg_table)
#'
#'
#' # 2) runing level two model (L2)
#' # Calling provided data from LINCS database
#'
#'\dontrun{
#' DComboNet(load_dir = load_dir,
#'           resultdir = resultdir,
#'           model = "L2",
#'           manual_input = FALSE,
#'           drugcandidate = drugcandidate,
#'           drugnetWeight = TRUE,
#'           featuretype = "integrated_score",
#'           dataset = "92742",
#'           cellline = "HEPG2",
#'           treatment_time = 6,
#'           foldchange_DEG = 0.5,
#'           pvalue_DEG = 0.05,
#'           foldchange_DEP = 0.5,
#'           pvalue_DEP = 0.05)
#' }
#'
#'
#'# Providing drug-DEG table, drug-DEP table and cancer sample/cell line specfic
#'# expressed gene list.
#'
#'\dontrun{
#'
#' drugDEG = read.table(paste0(load_dir, "/test/Sorafenib_DEG_test.txt"), sep =
#' "\t", header = T, stringsAsFactors = F)
#' drugDEP = read.table(paste0(load_dir, "/test/Sorafenib_DEP_test.txt"), sep =
#' "\t", header = T, stringsAsFactors = F)
#' cancergene = read.table(paste0(load_dir, "/test/HEPG2_genelist.csv"), sep =
#' ",", header = T, stringsAsFactors = F)
#' DComboNet(load_dir = load_dir,
#'           resultdir = resultdir,
#'           model = "L2",
#'           manual_input = FALSE,
#'           drugcandidate = drugcandidate,
#'           drugnetWeight = TRUE,
#'           featuretype = "integrated_score",
#'           cellline = "HEPG2",
#'           drugDEG = drugDEG,
#'           drugDEP = drugDEP,
#'           cancergene = cancergene)
#' }
#'
#'


DComboNet <- function(load_dir,
                      resultdir,
                      model = c("L1","L2"),
                      drugcandidate = NULL,
                      manual_input = c("TRUE","T","FALSE","F"),
                      CDK_FP = NULL,
                      pubchemFP = NULL,
                      MACCS_FP = NULL,
                      drugnetWeight = c("TRUE","T","FALSE","F"),
                      featuretype = c("ATCsim","STsim","SEsim","integrated_score"),
                      drugtarget = NULL,
                      druggene = NULL,
                      dataset = NULL,
                      cellline = NULL,
                      treatment_time = NULL,
                      foldchange_DEG = NULL,
                      pvalue_DEG = NULL,
                      foldchange_DEP = NULL,
                      pvalue_DEP = NULL,
                      drugDEG = NULL,
                      drugDEP = NULL,
                      cancergene = NULL,
                      dtweight = 2,
                      dgweight = 1,
                      dDEGweight = 1,
                      dgpweight = 1,
                      dDEpweight = 1,
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
    # devtools::load_all()


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
    c.transitionMatrix = compiler::cmpfun(TransitionMatrix)

    #loading internal data
    # load(system.file("data", "druglistlib.RData", package = "DComboNet"))
    # load(system.file("data", "pathway_pathway_interaction.RData", package = "DComboNet"))


    #timestart<-Sys.time()
    print("----- Drug Net Features Calculating -----")

    # druglistlib = read.csv(paste0(load_dir,"data/druglist.csv"),stringsAsFactors = F)

    # if(is.null(drugcandidate) |length(setdiff(drugcandidate[,1], druglistlib[,1])) == 0 ){

      if(model == "L1"){

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

#     feature.path =  "G:/lab/Projects/p1_DComboNet/DComboNet_improve/Rpackage/DComboNet-master/inst/extdata/drug_net/features_OCILY3.csv"
#     drugnet_feature = read.table(feature.path, sep = ",", header = T, stringsAsFactors = F)

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

    # if(is.null(drugDEG) & (model == "L2" | (is.null(cellline)==FALSE))){
    #
    #
    #   drugDEG_preparation(cellline = cellline,
    #                       dataset = dataset,
    #                       treatment_time =  treatment_time,
    #                       foldchange = foldchange_DEG,
    #                       pvalue = pvalue_DEG,
    #                       load_dir = load_dir)
    #
    #   drugDEP_preparation(cellline = cellline,
    #                       dataset = dataset,
    #                       treatment_time =  treatment_time,
    #                       foldchange = foldchange_DEP,
    #                       pvalue = pvalue_DEP,
    #                       load_dir = load_dir)
    # }

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
                                                dDEpweight = 1,
                                                PathwayNetwork = WWI_all,
                                                load_dir = load_dir)
    }


    # genepathwayMatrix

    gpAdj = genepathwayAdj(drugAdj = drugnet_adj,
                           pathwayadj = pathwayadj,
                           geneadj = geneadj,
                           gpweight = 1,
                           load_dir = load_dir)



    dgTranMatrix <- TransitionMatrix(drugAdj = drugnet_adj,
                                     geneAdj = geneadj,
                                     pathwayAdj = pathwayadj,
                                     druggeneAdj = dgAdj,
                                     genepathwayAdj = gpAdj,
                                     drugpathwayAdj = drugpathwayMatrix,
                                     x=x,
                                     y=y,
                                     z=z,
                                     A=A,
                                     B=B)

    multiplex_m =  reshape2::melt(as.matrix(dgTranMatrix))
    multiplex_df=multiplex_m[multiplex_m$value!=0,]
    genenet = reshape2::melt(as.matrix(geneadj))
    genenet=genenet[genenet$value!=0,]
    genepathwaynet = reshape2::melt(as.matrix(gpAdj))
    genepathwaynet=genepathwaynet[genepathwaynet$value!=0,]

    if(model =="L1"){

      save(multiplex_df, genenet, genepathwaynet, file = paste0(resultdir,model,"_result/potential_net/net_data.RData"))

    }else if(model == "L2"){

      save(multiplex_df, genenet, genepathwaynet, file = paste0(resultdir,model,"_result/potential_net/",cellline,"_net_data.RData"))
    }

    N.gene = nrow(geneadj)
    N.drug = nrow(drugnet_adj)



    print("----- Preparaing Drug Seed -----")
    if(is.null(drugcandidate)){

      if(manual_input %in% c("T","TRUE")  ){
        print("Please type in the drug you are interested: ")
        drugcandidate = data.frame(readline())
        drugcandidate = Hmisc::capitalize(tolower(drugcandidate))
      }else{
        drugcandidate = druglistlib # read.csv(paste0(load_dir,"data/druglist.csv"), header = T, stringsAsFactors = F)

      }
    }

    for(drugseeds in drugcandidate[,1]){

      if(drugseeds %in% colnames(drugnet_adj)){

        if(model == "L1" | is.null(cellline)){

          print("------ Level 1 Model ------")

          seed_score = seedscore(seeds = drugseeds, eta = eta)

          rwr_result = rwr(tm = dgTranMatrix,
                           r = r,
                           seeds_score = seed_score)

          drugs_rank = rank_drugs(Num.Gene = N.gene,
                                  Num.Drug = N.drug,
                                  RWR.result = rwr_result,
                                  Drug_seeds = drugseeds)

          genes_rank = rank_genes(Num.Gene = N.gene,
                                  Num.Drug = N.drug,
                                  RWR.result = rwr_result)

          pathways_rank <- rank_pathways(Num.Gene=N.gene,
                                         Num.Drug=N.drug,
                                         RWR.result=rwr_result)

        }else if(model == "L2" | (is.null(cellline)==FALSE)){

          print("------ Level 2 Model ------")

          seed_score = seedscore(seeds = drugseeds,
                                 eta = eta)

          rwr_result = rwr(tm = dgTranMatrix,
                           r = r,
                           seeds_score = seed_score)
          drugs_rank = rank_drugs(Num.Gene = N.gene,
                                  Num.Drug = N.drug,
                                  RWR.result = rwr_result,
                                  Drug_seeds = drugseeds)

          genes_rank = rank_genes(Num.Drug = N.drug,
                                  Num.Gene = N.gene,
                                  RWR.result = rwr_result)

          pathways_rank <- rank_pathways(Num.Gene=N.gene,
                                         Num.Drug=N.drug,
                                         RWR.result=rwr_result)


        }
        drugs_rank$rank = 1:nrow(drugs_rank)
        genes_rank$rank = 1:nrow(genes_rank)
        pathways_rank$rank = 1:nrow(pathways_rank)

        #colnames(rwr_result) = drugseeds
        write.csv(drugs_rank,paste0(resultdir,model,"_result/drugrank/",drugseeds,"_rank.csv"), row.names = F, quote = F)
        write.csv(genes_rank,paste0(resultdir,model,"_result/generank/",drugseeds,"_rank.csv"), row.names = F, quote = F)
        write.csv(pathways_rank,paste0(resultdir,model,"_result/pathwayrank/",drugseeds,"_rank.csv"), row.names = F, quote = F)
        print(paste0("The prediction result of drug name ",drugseeds," has been saved!"))
        #return(drugs_rank)

      }else{
        print(paste0("The drug name ",drugseeds," can not be found in the drug network! "))
      }
    }
  }
}

#' Ranking system of DComboNet
#' @description the function \code{DComboNet} is to predict possible drug to
#'   pair with the drug candidate without expression data
#' @param model choose model type, two models here to choose: "L1", "L2"
#' @param drugcandidate drug and its corresponding "drugbank ID" in
#'   \code{data.frame}
#' @param resultdir path to where result files saved
#' @param filename
#' @param outputdir
#' @param methods choose method used to set the final rank. Two methods are
#'   are available to choose: "sum" and "two_threshold". Default method is
#'   "sum".
#' @param threshold_syn if \code{methods = "two_threshold"}, the threshold used
#'   to determine the combinable pairs.
#' @param threshold_ant if \code{methods = "two_threshold"}, the threshold used
#'   to determine the uncombinable pairs.
#' @return Ranked result tables. Each drug seed will create three files that
#'   contains global proximity between drug seed and other drug/gene/pathway
#'   nodes computed from drug-gene/pathway interaction network. File will be
#'   saved under the path \code{resultdir} provided by user.
#' @export
#'
#' @examples
#' \dontrun{
#' result_dir = 'G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_tryout_results/result_overlap/'
#' outputdir = 'G:/lab/Projects/p1_DComboNet/DComboNet_improve/OCI_LY3_tryout_results/'
#' Ranking(
#'  drugcandidate = druglist,
#'  model = "L2",
#'  resultdir = result_dir,
#'  filename = 'test',
#'  outputdir = outputdir,
#'  methods = "sum",
#'  threshold_syn = 0.1,
#'  threshold_ant = 0.5
#' )}
#'

Ranking <- function(
  drugcandidate,
  model,
  resultdir,
  filename = NULL,
  outputdir = NULL,
  methods = c("sum", "two_threshold"),
  threshold_syn = 0.1,
  threshold_ant = 0.5
) {

  if(methods == "two_threshold" & (is.null(threshold_syn) | is.null(threshold_ant) )){
    print("ERROR: if ranking method 'two_threshold' has been chosen, 'threshold_syn' and 'threshold_ant' must be provided.")
  }


  a = list.files(paste0(resultdir,model,'_result/drugrank/'))
  dir = paste0(resultdir,model,'_result/drugrank/',a)
  n = length(a)
  b = gsub('_rank.csv','',a)
  drugrank = read.csv(dir[1],stringsAsFactors = F)
  drugrank$drugseed = b[1]
  for(i in 2:n){
    drugrank1 = read.csv(dir[i],stringsAsFactors = F)
    drugrank1$drugseed = b[i]
    drugrank = rbind(drugrank,drugrank1)
  }
  drugrank = unique(drugrank[c(4,1,2,3)])
  names(drugrank) = c('drugseed','DrugID','Score','Rank')
  drugrank$drugseed = toupper(drugrank$drugseed)
  drugrank$DrugID = toupper(drugrank$DrugID)
  names(drugrank) = c('A','B','Score','Rank')

  drugpair = matrix(0,length(drugcandidate),length(drugcandidate));
  colnames(drugpair) = drugcandidate
  rownames(drugpair) = drugcandidate
  drugpair[lower.tri(drugpair)] = 1
  #diag(drugpair) <- 0
  drugpair = reshape2::melt(as.matrix(drugpair))
  drugrank1 = drugpair[drugpair$value == 1,][c(1,2)]
  names(drugrank1) = c('A','B')
  drugrank1$A = toupper(drugrank1$A)
  drugrank1$B = toupper(drugrank1$B)

  drugrank1 = merge(drugrank[c(1,2,3)], drugrank1, by =c('A','B'))
  drugrank1$A_B = as.numeric(0)
  drugrank1$B_A = as.numeric(0)
  drugrank1$A_B_score = as.numeric(0)
  drugrank1$B_A_score = as.numeric(0)
  drugrank1$pred_TAG = 0
  drugrank1$Score_two_threshold = 0
  drugrank1$Rank_two_threshold = 0
  num_drug = max(drugrank$Rank)

  for(i in 1:nrow(drugrank1)){

    seedA = drugrank[drugrank$A == drugrank1[i,1] & drugrank$B == drugrank1[i,2],]
    seedB = drugrank[drugrank$B == drugrank1[i,1] & drugrank$A == drugrank1[i,2],]

    if(nrow(seedA)>0 & nrow(seedB)>0){

      drugrank1[i,]$A_B = round(mean(seedA$Rank),0)
      drugrank1[i,]$B_A = round(mean(seedB$Rank),0)

      drugrank1[i,]$A_B_score = mean(seedA$Score)
      drugrank1[i,]$B_A_score = mean(seedB$Score)

      if(drugrank1[i,]$A_B <= num_drug*threshold_syn & drugrank1[i,]$B_A <= num_drug*threshold_syn){

        drugrank1[i,]$pred_TAG = 'Syn'
        drugrank1[i,]$Score_two_threshold = max(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = min(seedA$Rank, seedB$Rank)
      }else if((drugrank1[i,]$A_B <= num_drug*threshold_syn & (drugrank1[i,]$B_A >= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_ant))|
               ((drugrank1[i,]$A_B >= num_drug*threshold_syn & drugrank1[i,]$A_B <= num_drug*threshold_ant) & drugrank1[i,]$B_A  <= num_drug*threshold_syn)){

        drugrank1[i,]$pred_TAG = 'Syn'
        drugrank1[i,]$Score_two_threshold = max(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = min(seedA$Rank, seedB$Rank)
      }else if(drugrank1[i,]$A_B >= num_drug*threshold_ant & drugrank1[i,]$B_A >= num_drug*threshold_ant){

        drugrank1[i,]$pred_TAG = 'Ant'
        drugrank1[i,]$Score_two_threshold = min(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = max(seedA$Rank, seedB$Rank)
      }else if((drugrank1[i,]$A_B >= num_drug*threshold_ant & (drugrank1[i,]$B_A >= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_ant))|
               ((drugrank1[i,]$A_B >= num_drug*threshold_syn & drugrank1[i,]$A_B <= num_drug*threshold_ant) & drugrank1[i,]$B_A  >= num_drug*threshold_ant)){

        drugrank1[i,]$pred_TAG = 'Ant'
        drugrank1[i,]$Score_two_threshold = min(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = max(seedA$Rank, seedB$Rank)
      }else if((drugrank1[i,]$A_B >= num_drug*threshold_syn & drugrank1[i,]$B_A <= num_drug*threshold_ant)|
               (drugrank1[i,]$A_B <= num_drug*threshold_ant & drugrank1[i,]$B_A >= num_drug*threshold_syn)){

        drugrank1[i,]$pred_TAG = 'Add'
        drugrank1[i,]$Score_two_threshold = mean(mean(seedA$Score), mean(seedB$Score))
        drugrank1[i,]$Rank_two_threshold = mean(round(mean(seedA$Rank),0), round(mean(seedB$Rank),0))
      }else if((drugrank1[i,]$A_B <= num_drug*threshold_syn & drugrank1[i,]$B_A >= num_drug*threshold_ant)|
               (drugrank1[i,]$A_B >= num_drug*threshold_ant & drugrank1[i,]$B_A <= num_drug*threshold_syn)){

        drugrank1[i,]$pred_TAG = 'Add'
        drugrank1[i,]$Score_two_threshold = mean(mean(seedA$Score), mean(seedB$Score))
        drugrank1[i,]$Rank_two_threshold = mean(round(mean(seedA$Rank),0), round(mean(seedB$Rank),0))

      }
    }else if(nrow(seedA)>0 & nrow(seedB)==0){

      drugrank1[i,]$A_B = round(mean(seedA$Rank),0)
      drugrank1[i,]$B_A = round(mean(seedA$Rank),0)
      drugrank1[i,]$A_B_score = mean(seedA$Score)
      drugrank1[i,]$B_A_score = mean(seedB$Score)

      if(drugrank1[i,]$A_B <= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_syn){

        drugrank1[i,]$pred_TAG = 'Syn'
        drugrank1[i,]$Score_two_threshold = max(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = min(seedA$Rank, seedB$Rank)
      }else if((drugrank1[i,]$A_B <= num_drug*threshold_syn & (drugrank1[i,]$B_A  >= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_ant))|
               ((drugrank1[i,]$A_B >= num_drug*threshold_syn & drugrank1[i,]$A_B <= num_drug*threshold_ant) & drugrank1[i,]$B_A  <= num_drug*threshold_syn)){

        drugrank1[i,]$pred_TAG = 'Syn'
        drugrank1[i,]$Score_two_threshold = max(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = min(seedA$Rank, seedB$Rank)
      }else if(drugrank1[i,]$A_B >= num_drug*threshold_ant & drugrank1[i,]$B_A >= num_drug*threshold_ant){

        drugrank1[i,]$pred_TAG = 'Ant'
        drugrank1[i,]$Score_two_threshold = min(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = max(seedA$Rank, seedB$Rank)
      }else if((drugrank1[i,]$A_B >= num_drug*threshold_ant & (drugrank1[i,]$B_A  >= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_ant))|
               ((drugrank1[i,]$A_B >= num_drug*threshold_syn & drugrank1[i,]$A_B <= num_drug*threshold_ant) & drugrank1[i,]$B_A  >= num_drug*threshold_ant)){

        drugrank1[i,]$pred_TAG = 'Ant'
        drugrank1[i,]$Score_two_threshold = min(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = max(seedA$Rank, seedB$Rank)
      }else if((drugrank1[i,]$A_B >= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_ant)|
               (drugrank1[i,]$A_B <= num_drug*threshold_ant & drugrank1[i,]$B_A  >= num_drug*threshold_syn)){

        drugrank1[i,]$pred_TAG = 'Add'
        drugrank1[i,]$Score_two_threshold = mean(seedA$Score)
        drugrank1[i,]$Rank_two_threshold = round(mean(seedA$Rank),0)
      }else if((drugrank1[i,]$A_B <= num_drug*threshold_syn & drugrank1[i,]$B_A  >= num_drug*threshold_ant)|
               (drugrank1[i,]$A_B >= num_drug*threshold_ant & drugrank1[i,]$B_A  <= num_drug*threshold_syn)){

        drugrank1[i,]$pred_TAG = 'Add'
        drugrank1[i,]$Score_two_threshold = mean(seedA$Score)
        drugrank1[i,]$Rank_two_threshold = round(mean(seedA$Rank),0)
      }
    }else if(nrow(seedA)==0 & nrow(seedB)>0){

      drugrank1[i,]$A_B = round(mean(seedB$Rank),0)
      drugrank1[i,]$B_A = round(mean(seedB$Rank),0)
      drugrank1[i,]$A_B_score = mean(seedA$Score)
      drugrank1[i,]$B_A_score = mean(seedB$Score)

      if(drugrank1[i,]$A_B <= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_syn){

        drugrank1[i,]$pred_TAG = 'Syn'
        drugrank1[i,]$Score_two_threshold = max(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = min(seedA$Rank, seedB$Rank)
      }else if((drugrank1[i,]$A_B <= num_drug*threshold_syn & (drugrank1[i,]$B_A  >= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_ant))|
               ((drugrank1[i,]$A_B >= num_drug*threshold_syn & drugrank1[i,]$A_B <= num_drug*threshold_ant) & drugrank1[i,]$B_A  <= num_drug*threshold_syn)){

        drugrank1[i,]$pred_TAG = 'Syn'
        drugrank1[i,]$Score_two_threshold = max(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = min(seedA$Rank, seedB$Rank)
      }else if(drugrank1[i,]$A_B >= num_drug*threshold_ant & drugrank1[i,]$B_A >= num_drug*threshold_ant){

        drugrank1[i,]$pred_TAG = 'Ant'
        drugrank1[i,]$Score_two_threshold = min(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = max(seedA$Rank, seedB$Rank)
      }else if((drugrank1[i,]$A_B >= num_drug*threshold_ant & (drugrank1[i,]$B_A  >= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_ant))|
               ((drugrank1[i,]$A_B >= num_drug*threshold_syn & drugrank1[i,]$A_B <= num_drug*threshold_ant) & drugrank1[i,]$B_A  >= num_drug*threshold_ant)){

        drugrank1[i,]$pred_TAG = 'Ant'
        drugrank1[i,]$Score_two_threshold = min(seedA$Score, seedB$Score)
        drugrank1[i,]$Rank_two_threshold = max(seedA$Rank, seedB$Rank)
      }else if((drugrank1[i,]$A_B >= num_drug*threshold_syn & drugrank1[i,]$B_A  <= num_drug*threshold_ant)|
               (drugrank1[i,]$A_B <= num_drug*threshold_ant & drugrank1[i,]$B_A  >= num_drug*threshold_syn)){

        drugrank1[i,]$pred_TAG = 'Add'
        drugrank1[i,]$Score_two_threshold = mean(seedB$Score)
        drugrank1[i,]$Rank_two_threshold = round(mean(seedB$Rank),0)
      }else if((drugrank1[i,]$A_B <= num_drug*threshold_syn & drugrank1[i,]$B_A  >= num_drug*threshold_ant)|
               (drugrank1[i,]$A_B >= num_drug*threshold_ant & drugrank1[i,]$B_A  <= num_drug*threshold_syn)){

        drugrank1[i,]$pred_TAG = 'Add'
        drugrank1[i,]$Score_two_threshold = mean(seedB$Score)
        drugrank1[i,]$Rank_two_threshold = round(mean(seedB$Rank),0)
      }
    }else{
      print('not predicted')
    }

  }

  drugrank1 = drugrank1[order(drugrank1$Score,decreasing=T),]
  drugrank1$rank2 = 1:nrow(drugrank1)

  if (methods == "sum") {

    drugrank1$Score_sum = drugrank1$A_B_score + drugrank1$B_A_score
    drugrank1 = drugrank1[order(drugrank1$Score_sum, decreasing = T), ]
    drugrank1$Rank_sum = 1:nrow(drugrank1)
    drugpred_res = drugrank1[c('A', 'B' ,'Score_sum', 'Rank_sum')]

  } else if (methods == "two_threshold") {

    drugpred_res = drugrank1[c('A', 'B' ,'Score_two_threshold', 'Rank_two_threshold')]
  }

  if(is.null(outputdir)){

    return(drugpred_res)

  }else{

    if(is.null(filename)){
      write.csv(drugpred_res, paste(outputdir,'drugpred_res.csv', sep = ''), quote = F, row.names = F)
    }else{
      write.csv(drugpred_res, paste(outputdir, filename ,'.csv', sep = ''), quote = F, row.names = F)
    }
    return(drugpred_res)
  }

}


