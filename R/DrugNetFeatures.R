# ATCsim:
# Drug ATC code similarity Calculation

#' Drug ATC Code Distance file preparation
#' @description The function \code{drugATCprep} is used for preparing input
#'   files for calculating drug ATC code distance in Anatomical Therapeutic
#'   Chemical Classification System (ATC system). A drug-ATC file is packed in
#'   the function for newly inputed drug(s) to find its ATC codes
#' @param druglist a list of drugs in \code{data.frame} format
#' @param load_dir path to load or save modelling data files
#' @return \code{drug_atc} a data frame of drug and its corresponding ATC code
#'   in ATC system
#' @export
#'
#' @examples
#' \dontrun{
#' druglist = data.frame(drug = "Sorafenib")
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' drugATCprep(druglist = druglist, load_dir = load_dir)}
#'

drugATCprep <- function(druglist = NULL,
                        load_dir){

  # load(system.file("data", "drug_atc_lib.RData", package = "DComboNet"))
  # load(system.file("data", "druglistlib.RData", package = "DComboNet"))

  if(!is.null(druglist)){
    druglist = data.frame(drugname = union(druglist[,1],druglistlib[,1]))

    names(all_drug_atc) <- c("ATC_code","drugname")
    drug_atc <- merge(druglist, all_drug_atc, by = "drugname")

    dir.create(paste0(load_dir,"drug_net/"))
    dir.create(paste0(load_dir,"drug_net/atc_sim"))
    write.table(drug_atc, paste0(load_dir, "drug_net/atc_sim/drug_atc.txt"), sep = "\t", col.names = T, row.names = F)

  }else{
    druglist = data.frame(drugname = unique(druglistlib[,1]))
    names(all_drug_atc) <- c("ATC_code","drugname")
    drug_atc <- merge(druglist, all_drug_atc, by = "drugname")

    dir.create(paste0(load_dir,"drug_net/"))
    dir.create(paste0(load_dir,"drug_net/atc_sim"))
    write.table(drug_atc, paste0(load_dir,"drug_net/atc_sim/drug_atc.txt"), sep = "\t", col.names = T, row.names = F)

  }
  return(drug_atc[c("drugname","ATC_code")])

}



#' Drug ATC Code Distance Matrix Function
#' @description The function \code{ATC.Dis} is to calculate drug ATC code distance in
#'   Anatomical Therapeutic Chemical Classification System.
#' @param drug_atc a data frame of drug and its corresponding ATC code
#'   in ATC system, generated from function \code{\link{drugATCprep}}
#' @param load_dir path to load or save modelling data files
#' @import igraph
#' @return \code{ATC.Dis} returns ATC code distance matrix, an intermediate file
#'   for computing drug-drug ATC code based similarity
#' @export
#' @seealso \code{\link{drugATCprep}} to prepare input drug-atc file
#' @examples
#'
#' \dontrun{
#' druglist = data.frame(drug = "Sorafenib")
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' drug_atc = drugATCprep(druglist = druglist, load_dir = load_dir)
#' ATC_dis = ATC.Dis(drug_atc = drug_atc, load_dir = load_dir)
#' }
#'

ATC.Dis <-function(drug_atc,
                   load_dir){

  if (requireNamespace("igraph")){

    # load(system.file("data", "atc_net.RData", package = "DComboNet"))
    # ATC <- igraph::read.graph(paste0(load_dir,"/data/ATC_NET/ATC_network.net"),format="pajek")
    # nodes <- read.csv(paste0(load_dir,"/data/ATC_NET/ATC_node.txt"),sep=" ",header=F, stringsAsFactors = F)

    drug_atc = drug_atc[,2]
    atc_freq <- 1/table(unlist(drug_atc))
    drug_atc.index <- rep(0,length(drug_atc))
    for(i in 1:length(drug_atc)){
      if(drug_atc[i] %in% nodes[,2]){
        drug_atc.index[i] <- nodes[nodes[,2]==drug_atc[i],1]
      }
    }

    Shortest.path <- igraph::shortest.paths(ATC,
                                            v=drug_atc.index[drug_atc.index!=0],
                                            to=drug_atc.index[drug_atc.index!=0])
    rownames(Shortest.path) <- drug_atc[drug_atc.index!=0]
    SP <- matrix(0,nrow=length(drug_atc),ncol=length(drug_atc))
    rownames(SP) <- drug_atc
    colnames(SP) <- drug_atc
    for(i in rownames(SP)){
      for(j in rownames(SP)){
        if(i %in% rownames(Shortest.path) && j %in% rownames(Shortest.path))
          SP[which(rownames(SP)==i),which(rownames(SP)==j)] <- atc_freq[i]*atc_freq[j]*exp(-0.25*Shortest.path[which(rownames(Shortest.path)==i),which(rownames(Shortest.path)==j)])
        else
          SP[which(rownames(SP)==i),which(rownames(SP)==j)] <- 0
      }
    }
  }
  return(SP)
}



#' Drug ATC Code-Based Similarity Function
#' @description The function \code{ATC.sim} is used to calculate drug-drug ATC code-based
#'   similarity.
#' @details \code{ATC.sim} integrated \code{\link{ATC.Dis}} and
#'   \code{\link{drugATCprep}} functions. \code{druglist} in \code{data.frame}
#'   should be provided. Their ATC code(s) will be found via \code{drugATCprep}.
#' @details \code{ATC.sim} will compute shortest distance between drugs based on
#'   the shortest distance between ATC codes  within ATC hierarchical network:
#'   \deqn{sim_{ATC}(t_{i},t_{j}) = w(t_{i})*w(t_{j})e^(-\gamma*d(t_{i},t_{j}))}
#'   where t_{i} and t_{j} denote the ATC code i and j, d(t_{i},t_{j}) denotes the
#'   shortest distance between t_{i} and t_{j} and computes with function
#'   \code{\link{ATC.Dis}}; w(t_{i}) and w(t_{j}) represent the weights of the
#'   corresponding ATC codes, and are defined as the inverse of ATC-code
#'   frequencies. \eqn{\gamma} is a pre-defined parameter (set as 0.25). For
#'   drugs with multiple ATC codes, average ATC code similarities will be
#'   considered.
#' @details ATC code based similarity between two drugs will be calculated and
#'   saved in \code{load_dir} path provided by user.
#'
#' @param druglist a list of drugs in \code{data frame} format
#' @param load_dir path to load or save modeling data files
#' @import compiler
#' @return \code{ATC.sim} returns \code{ATCsim.matrix}, drug-drug ATC code based
#'   similarity matrix. Row name and column name are drug names, values
#'   represent drug-drug ATC code based similarity. Matrix will be saved in
#'   provided path \code{load_dir} under \emph{"/drug_net/atc_sim/"} folder, under the
#'   name of "ATCsim.txt"
#' @export
#' @seealso \code{\link{drugATCprep}} to prepare input drug_atc file,
#'   \code{\link{ATC.Dis}} for computing intermediated ATC distance matrix
#'
#' @examples
#'
#' \dontrun{
#' druglist = data.frame(drug = "Sorafenib")
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' ATC.sim(druglist = druglist, load_dir = load_dir)
#' }
#'


ATC.sim <- function(druglist = NULL,
                    load_dir){

  if(requireNamespace("compiler")){

    c.drugATCprep = compiler::cmpfun(drugATCprep)
    c.ATC.Dis = compiler::cmpfun(ATC.Dis)
    drug_atc = c.drugATCprep(druglist,load_dir)
    DT.M <- as.matrix(table(drug_atc))
    ATCD.M <- c.ATC.Dis(drug_atc,load_dir)
    N.drug <- nrow(DT.M)
    ATCsim.matrix <- matrix(0,nrow=N.drug,ncol=N.drug)
    rownames(ATCsim.matrix) <- rownames(DT.M)
    colnames(ATCsim.matrix) <- rownames(DT.M)
    for(i in 1:N.drug){
      for(j in i:N.drug){
        ATCsim.matrix[i,j] <- max(ATCD.M[DT.M[i,]>0,DT.M[j,]>0])
        ATCsim.matrix[j,i] <- ATCsim.matrix[i,j]
      }
    }
    write.table(ATCsim.matrix,paste0(load_dir,"drug_net/atc_sim/ATCsim.txt"),sep="\t",quote = F)
    return(ATCsim.matrix)
  }
}

#' Tanimoto Coefficient between Drug Chemical Structure fingerprints
#'
#' @description The function \code{FP.tanimoro} is used for computing Tanimoto
#'   Coefficient between drug chemical structure fingerprints. The Chemical
#'   fingerprints can be generated from multiple softwares. Here we use PEDAL
#'   generated fingerprint file.
#'
#' @param fingerprint Drug fingerprint matrix (binary 0-1 matrix) generated via
#'   PEDAL, rows denote drugs and columns denote particular substructures in the
#'   molecule.
#'
#' @return tanimoto_score.matrix
#' @export
#'
#' @examples
#' \dontrun{
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' fingerprint = read.table(paste0(load_dir,
#' "/data/fingerprints/fingerprints.csv"), sep = ",", header = TRUE,
#' stringsAsFactors = FALSE)
#' FP.tanimoro(fingerprint = fingerprint)
#' }
#'

FP.tanimoro <- function(fingerprint){
  tanimoto_score.matrix <- matrix(0,nrow=nrow(fingerprint),ncol=nrow(fingerprint))
  row.names(tanimoto_score.matrix) <- rownames(fingerprint)
  colnames(tanimoto_score.matrix) <- rownames(fingerprint)
  for(i in 1:nrow(fingerprint)){
    for(j in i:nrow(fingerprint)){
      fp1 <- as.numeric(unlist(fingerprint[i,]))
      fp2 <- as.numeric(unlist(fingerprint[j,]))
      A = sum((fp1 | fp2) & !fp2)
      B = sum((fp1 | fp2) & !fp1)
      C = sum(fp1 & fp2)
      tanimoto_score.matrix[i,j] = C/(A + B + C)
      tanimoto_score.matrix[j,i] = tanimoto_score.matrix[i,j]
    }
  }
  return(tanimoto_score.matrix)
}


#' Drug Chemical Fingerprint Similarity Function
#' @description The function \code{structure.sim} is for computing drug chemical
#'   structure fingerprint based similarity.
#'@details Three kinds of chemical fingerprint based similarity will be computed
#'  and integrated into one score. Tanimoto index is used to calculate the
#'  similarity matrix. The Chemical fingerprints can be generated from multiple
#'  softwares, here we use PEDAL to generate three kinds fingerprint matrix,
#'  including CDK fingerprint, pubchem fingerprint and MACCS fingerprint. The
#'  average score of these three fingerprint based similarity is computed as
#'  final drug-drug chemical fingerprint-based similarity. Four similarity
#'  matrices will be calculated and saved to path provided by user
#'  (\code{load_dir}) under \emph{"/data/fingerprints/"} folder.
#' @param CDK_FP Drug fingerprint matrix generated from PEDAL
#' @param pubchemFP Drug Pubchem fingerprint matrix generated from PEDAL
#' @param MACCS_FP Drug MACCS fingerprint matrix generated from PEDAL
#' @param load_dir path to load and save model data file(s)
#' @import compiler
#' @return \code{structure_sim} drug-drug chemical fingerprint-based similarity
#'   matrix, rows and columns represent drug names
#' @export
#'
#' @examples
#' \dontrun{
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' structure.sim(load_dir = load_dir)
#' }
#'
#'

structure.sim <- function(CDK_FP=NULL,
                          pubchemFP=NULL,
                          MACCS_FP=NULL,
                          load_dir){

  if (requireNamespace("compiler")){

    c.FP.tanimoro = compiler::cmpfun(FP.tanimoro)

    # load(system.file("data", "fingerprint_lib.RData", package = "DComboNet"))
    # CDK_FP_lib = read.csv(paste0(load_dir,"/data/fingerprints/fingerprints.csv"),stringsAsFactors = F)
    # pubchemFP_lib = read.csv(paste0(load_dir,"/data/fingerprints/pubchem_fingerprints.csv"),stringsAsFactors = F)
    # MACCS_FP_lib = read.csv(paste0(load_dir,"/data/fingerprints/MACCS_fingerprints.csv"),stringsAsFactors = F)

    if((!is.null(CDK_FP))& (!is.null(pubchemFP))& (!is.null(MACCS_FP))){

      CDK_FP = unique(rbind(CDK_FP_lib,CDK_FP))
      pubchemFP = unique(rbind(pubchemFP_lib,pubchemFP))
      MACCS_FP = unique(rbind(MACCS_FP_lib,MACCS_FP))

    }else{

      CDK_FP = CDK_FP_lib
      pubchemFP = pubchemFP_lib
      MACCS_FP = MACCS_FP_lib

    }

    rownames(CDK_FP) = CDK_FP[,1]
    rownames(pubchemFP) = pubchemFP[,1]
    rownames(MACCS_FP) = MACCS_FP[,1]

    CDK_FP_tanimoto_dis.matrix <- c.FP.tanimoro(CDK_FP[-1])
    pubchemFP_tanimoto_dis.matrix <- c.FP.tanimoro(pubchemFP[-1])
    MACCS_FP_tanimoto_dis.matrix <- c.FP.tanimoro(MACCS_FP[-1])
    structure_sim <- (CDK_FP_tanimoto_dis.matrix + pubchemFP_tanimoto_dis.matrix + MACCS_FP_tanimoto_dis.matrix)/3
    dir.create(paste0(load_dir,"drug_net/structure_sim"))

    write.table(CDK_FP_tanimoto_dis.matrix, paste0(load_dir,"drug_net/structure_sim/CDK_FP_tanimoto.txt"), sep = "\t", col.names = T, row.names  = T, quote = F)
    write.table(pubchemFP_tanimoto_dis.matrix, paste0(load_dir,"drug_net/structure_sim/pubchemFP_tanimoto.txt"), sep = "\t", col.names = T, row.names = T, quote = F)
    write.table(MACCS_FP_tanimoto_dis.matrix, paste0(load_dir,"drug_net/structure_sim/MACCS_FP_tanimoto.txt"), sep = "\t", col.names = T, row.names =T, quote = F)
    write.table(structure_sim, paste0(load_dir,"drug_net/structure_sim/structure_sim.txt"), sep = "\t", col.names = T, row.names = T, quote = F)
    return(structure_sim)
  }
}



#' Drug side effect similarity file preparation
#' @description The function \code{drugSEprep} is used to prepare input file for
#'   function \code{SE.tanimoto}
#' @details \code{drugSEprep} provides drug-side effect file download from
#'   \code{SIDER4} database, where newly inputed drugs can be matched to their
#'   known side effect terms through "drugbank ID".
#' @param druglist drug and its corresponding "drugbank ID" in \code{data
#'   frame}
#' @param load_dir path to load or save modeling data files
#' @importFrom utils read.delim
#' @return \code{drug_SE} a data frame contains drugs and their side effect
#'   terms. If druglist contains newly inputed drug, their drug-sideeffect
#'   dataframe will be merged with provided data
#' @export
#' @examples
#'
#' \dontrun{
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' druglist = NULL
#' # or druglist = read.table("druglist_example.csv", sep = ",", header = TRUE,
#' # stringsAsFactors = FALSE)
#' drugSEprep(druglist = druglist, load_dir = load_dir)
#' }
#'

drugSEprep <- function(druglist = NULL,
                       load_dir){

  # all_side_effect = read.delim(paste0(load_dir,"data/side-effects.tsv"),sep="\t",header = T,stringsAsFactors = F)
  # drugselib = read.csv(paste0(load_dir,"data/drug_se.txt"),sep="\t",header=T, stringsAsFactors = F)

  # load(system.file("data", "drug_sideeffect.RData", package = "DComboNet"))

  if(is.null(druglist)){
    drug_SE = drugselib

  }else{
    names(druglist) = c("Drug","drugbank_id")
    drug_SE = merge(druglist,all_side_effect,by = "drugbank_id")
    drug_SE=drug_SE[c("Drug","side_effect_name")]
    drug_SE = rbind(drugselib,drug_SE)
  }

  dir.create(paste0(load_dir,"/drug_net/sideeffect_sim"))
  write.table(drug_SE, paste0(load_dir,"/drug_net/sideeffect_sim/drug_SE.txt"), sep = "\t", col.names = T, row.names =F, quote = F)
  return(drug_SE)
}


#' Tanimoto Coefficient between Drug side effect
#' @description The function \code{SE.sim} is used to generate drug side effect
#'   based similarity matrix
#' @details drug-side effect dataframe will be converted into 0-1 binary matrix
#'   where rows denote drugs columns denote side effect terms, 1 represent drug
#'   has corresponding side effect term extracted from \code{SIDER4} database
#'   while 0 represent does not has.  Tanimoto coefficient was considered to
#'   compute side effect-based similarity.
#' @param druglist a dataframe of drug and corresponding "drugbank ID"
#' @param load_dir path to load and save modeling data files
#' @return \code{SE.tanimoto.matrix} Drug similarity matrix based on drug side
#'   effects Tanimoto Coefficient.
#' @export
#' @examples
#'
#' \dontrun{
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' druglist = NULL
#' # or druglist = read.table("druglist_example.csv", sep = ",", header = TRUE,
#' # stringsAsFactors = FALSE)
#' SE.sim(druglist = druglist, load_dir = load_dir)
#' }
#'

SE.sim <- function(druglist=NULL, load_dir){

  if (requireNamespace("compiler")){

    c.drugSEprep = compiler::cmpfun(drugSEprep)

    # load(system.file("data", "drug_sideeffect.RData", package = "DComboNet"))

    if(!is.null(druglist)){

      drug_se = c.drugSEprep(druglist,load_dir)

    }else{

      drug_se = drugselib #read.csv(paste0(load_dir,"data/drug_se.txt"),sep="\t",header=T, stringsAsFactors = F)
    }
    SE.matrix <- as.matrix(table(drug_se))
    SE.tanimoto.matrix <- matrix(0,nrow=nrow(SE.matrix),ncol=nrow(SE.matrix))
    row.names(SE.tanimoto.matrix) <- rownames(SE.matrix)
    colnames(SE.tanimoto.matrix) <- rownames(SE.matrix)
    for(i in 1:nrow(SE.matrix)){
      for(j in i:nrow(SE.matrix)){
        fp1 <- as.numeric(unlist(SE.matrix[i,]))
        fp2 <- as.numeric(unlist(SE.matrix[j,]))
        A = sum((fp1 | fp2) & !fp2)
        B = sum((fp1 | fp2) & !fp1)
        C = sum(fp1 & fp2)
        SE.tanimoto.matrix[i,j] = C/(A + B + C)
        SE.tanimoto.matrix[j,i] = SE.tanimoto.matrix[i,j]
      }
    }
    write.table(SE.tanimoto.matrix, paste0(load_dir,"drug_net/sideeffect_sim/sideeffect_sim.txt"), sep = "\t", col.names = T, row.names = T, quote = F)
    return(SE.tanimoto.matrix)
  }
}


#' Drug drug interaction network edge weight: Drug-Drug similarity
#' @description The function \code{DrugNetFeature} is to calculate drug ATC
#'   code-based similarity, chemical structure fingerprint-based similarity and
#'   drug side effect-based similarity and integrate together as the final
#'   pharmacological similarity.
#' @details \code{DrugNetFeature} integrates calculation of three drug-drug
#'   similarity functions, \code{druglist} included "drugbank ID" is
#'   required to provide for file preparation. drug ATC code-based similarity,
#'   chemical structure fingerprint-based similarity and drug side effect-based
#'   similarity will be computed seperately.
#' @details \code{NA} generated due to lack of data will be imputed by the \code{mean}
#'   of each similarity score. Three simiarity score will be integrated into
#'   phamacological score and the score of known combinations will be replaced
#'   as \code{1}. Final \code{feature} table will be saved in provided
#'   \code{load_dir} path under \emph{/drug_net/} folder named with
#'   "features.csv" for level one model or "features_\emph{cellline}.csv" if
#'   level two model is called and \code{cellline} is provided.
#'
#'
#' @param druglist drug and its corresponding "drugbank ID" in \code{data
#'   frame}
#' @param cellline cancer cell lines name, provide only when level two model is
#'   called (\code{model = "L2"}).
#' @param model define which level model should be provided. Options are level
#'   one model ("\code{L1}") and level two model ("\code{L2}")
#' @param CDK_FP Drug fingerprint matrix generated from PEDAL
#' @param pubchemFP Drug Pubchem fingerprint matrix generated from PEDAL
#' @param MACCS_FP Drug MACCS fingerprint matrix generated from PEDAL
#' @param load_dir path to load and save modeling data files
#' @import reshape2
#' @import compiler
#' @return \code{features} table that includes the three drug similarity values
#'   and the integrated pharmacological score
#' @export
#'
#' @examples
#'
# \dontrun{
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' druglist = NULL
#' # or druglist = read.table("druglist_example.csv", sep = ",", header = TRUE,
#' # stringsAsFactors = FALSE)
#' # For level one model
#' cellline = NULL
#' model = "L1"
#'
#'
#' # For level two model
#' cellline = "HEPG2"
#' model = "L2"
#'
#'
#' # If no newly inputed drug(s): druglist = NULL
#' CDK_FP=NULL
#' pubchemFP=NULL
#' MACCS_FP=NULL
#' # Otherwise, please input PaDEL generated three fingerprint files
#' DrugNetFeature(druglist = druglist, cellline = cellline, model = model,
#' CDK_FP = CDK_FP, pubchemFP = pubchemFP, MACCS_FP = MACCS_FP, load_dir =
#' load_dir)
#' }
#'
#' @seealso \code{\link{ATC.sim}} for ATC code based similarity,
#'   \code{\link{structure.sim}} for chemical fingerprint based similarity,
#'   \code{\link{SE.sim}} for side effect based similarity.



DrugNetFeature <- function(druglist = NULL,
                           cellline = NULL,
                           model = c("L1","L2"),
                           CDK_FP=NULL,
                           pubchemFP=NULL,
                           MACCS_FP=NULL,
                           load_dir){

  if(requireNamespace(c("reshape2","compiler"))){

    # load(system.file("data", "drugpair_lib.RData", package = "DComboNet"))
    # load(system.file("data", "druglistlib.RData", package = "DComboNet"))

    c.ATC.sim = compiler::cmpfun(ATC.sim)
    c.structure.sim = compiler::cmpfun(structure.sim)
    c.SE.sim = compiler::cmpfun(SE.sim)

    # druglistlib = read.csv(paste0(load_dir,"data/druglist.csv"), header = T, stringsAsFactors = F)
    # drugpair_lib = read.csv(paste0(load_dir,"data/drugpair_base.csv"),stringsAsFactors = F)

    drugpair_P = drugpair_lib[drugpair_lib$TAG == "P",]
    drugpair_P$label = paste0(drugpair_P$A, "_",drugpair_P$B)
    drugpair_N = drugpair_lib[drugpair_lib$TAG == "N",]


    if(!is.null(druglist)){

      druglib_sum = union(druglist[,1], druglistlib[,1])

      dp_tmp = data.frame(A = rep(druglib_sum[1],length(druglib_sum)),B = druglib_sum,TAG = "N")

      for(i in 2:nrow(druglist)){
        #if(!(druglist[i,1]%in% druglib_sum)){

          dp_tmp2 =data.frame(A = rep(druglib_sum[i],length(druglib_sum)),B = druglib_sum,TAG = "N")
          dp_tmp = rbind(dp_tmp, dp_tmp2)

       # }
      }

      drugpair = unique(dp_tmp)
      drugpair$label = paste0(drugpair$A,"_", drugpair$B)
      drugpair = merge(drugpair, drugpair_P, by=c("A","B","label"), all.x = T)
      drugpair = unique(drugpair[-c(3,4)])
      names(drugpair) = c("A","B","TAG")
      drugpair[is.na(drugpair$TAG),]$TAG = "N"
      drugpair = unique(rbind(drugpair, drugpair_N))

      print("----- check if the input drugs exist in the basic drugnet -----")
      print(paste0("Basic drug number: ",nrow(druglistlib)))
      print(paste0("Input drug number: ",nrow(druglist)))

      if(length(intersect(druglistlib[,1],druglist[,1]))!=0){
        print(paste0("Overlaping drug number: ",length(intersect(druglistlib[,1],druglist[,1]))))
      }else{
        print("None of the input drugs is in the basic drugnet.")
      }

      print(paste0("Over all: ",length(druglib_sum)," drugs"))

      ATCsim_matrix = c.ATC.sim(druglist,load_dir)
      STsim_matrix = c.structure.sim(CDK_FP = CDK_FP,
                                     pubchemFP = pubchemFP,
                                     MACCS_FP = MACCS_FP,
                                     load_dir)
      SEsim_matrix = c.SE.sim(druglist, load_dir)
    }else{
      print("Prediction will run within our database ~~~")
      ATCsim_matrix = c.ATC.sim(druglist=NULL,load_dir)
      STsim_matrix = c.structure.sim(CDK_FP=NULL,
                                     pubchemFP=NULL,
                                     MACCS_FP=NULL,
                                     load_dir)
      SEsim_matrix = c.SE.sim(druglist=NULL, load_dir)

    }

    colnames(ATCsim_matrix) <- rownames(ATCsim_matrix)
    colnames(STsim_matrix) <- rownames(STsim_matrix)
    colnames(SEsim_matrix) <- rownames(SEsim_matrix)
    ATCsim_matrix[!upper.tri(ATCsim_matrix)] <- -1
    STsim_matrix[!upper.tri(STsim_matrix)] <- -1
    SEsim_matrix[!upper.tri(SEsim_matrix)] <- -1
    ATCsim_matrix_melt <- reshape2::melt(ATCsim_matrix); names(ATCsim_matrix_melt) <- c("A","B","ATCsim")
    STsim_matrix_melt <- reshape2::melt(STsim_matrix); names(STsim_matrix_melt) <- c("A","B","STsim")
    SEsim_matrix_melt <- reshape2::melt(SEsim_matrix); names(SEsim_matrix_melt) <- c("A","B","SEsim")
    ATCsim_matrix_melt = ATCsim_matrix_melt[ATCsim_matrix_melt$ATCsim != -1,]
    STsim_matrix_melt = STsim_matrix_melt[STsim_matrix_melt$STsim != -1,]
    SEsim_matrix_melt = SEsim_matrix_melt[SEsim_matrix_melt$SEsim != -1,]


    if(nrow(ATCsim_matrix_melt[ATCsim_matrix_melt$ATCsim == -Inf,])!=0){ATCsim_matrix_melt[ATCsim_matrix_melt$ATCsim == -Inf,]$ATCsim = 1/nrow(ATCsim_matrix_melt)}

    features <- unique(data.frame(merge(merge(merge(drugpair, STsim_matrix_melt, by = c("A","B"), all.x = T),
                                              ATCsim_matrix_melt, by = c("A","B"), all.x = T),
                                        SEsim_matrix_melt, by = c("A","B"), all.x = T)))


    if(nrow(features[is.na(features$ATCsim),])!=0){features[is.na(features$ATCsim),]$ATCsim = mean(features[!is.na(features$ATCsim),]$ATCsim)}
    if(nrow(features[is.na(features$STsim),])!=0){features[is.na(features$STsim),]$STsim =  mean(features[!is.na(features$STsim),]$STsim)}
    if(nrow(features[is.na(features$SEsim),])!=0){features[is.na(features$SEsim),]$SEsim = mean(features[!is.na(features$SEsim),]$SEsim)}

    if(nrow(features[is.na(features$ATCsim),])!=0){features[is.na(features$ATCsim),]$ATCsim = mean(features[!is.na(features$ATCsim),]$ATCsim)}
    if(nrow(features[is.na(features$STsim),])!=0){features[is.na(features$STsim),]$STsim = 0.1}
    if(nrow(features[is.na(features$SEsim),])!=0){features[is.na(features$SEsim),]$SEsim = mean(features[!is.na(features$SEsim),]$SEsim)}

    features = unique(features[c("A","B","ATCsim","STsim","SEsim","TAG")])
    score = 1-(1-features[,3])*(1-features[,4])*(1-features[,5])
    features$integrated_score = score
    features$integrated_score2 = score

    features = unique(features[c("A","B","ATCsim","STsim","SEsim","integrated_score","integrated_score2","TAG")])

    if(nrow(features[features$TAG == "P",]) != 0 ){features[features$TAG == "P",]$integrated_score = 1}
    features = features[features$integrated_score2 >= 0.2,]

    if(model == "L1" | is.null(cellline)){
      write.table(features,paste0(load_dir,"/drug_net/features.csv"),sep = ",",row.names = F, quote = F)
    }else if(model == "L2" | (is.null(cellline)==FALSE)){
      write.table(features,paste0(load_dir,"/drug_net/features_",cellline,".csv"),sep = ",",row.names = F, quote = F)
    }
    return(features)
  }
}

