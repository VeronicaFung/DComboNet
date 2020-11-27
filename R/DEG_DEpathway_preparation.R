#' Drug induced differentially expressed gene/pathway list preparation
#' @description The function \code{DEG_DEP_preparation} is to prepare drug
#'   induced differentially expressed gene table and pathway table for level two
#'   model(\code{model = "L2"}) from LINCS database. Differentially expressed
#'   genes were selected by functions \code{lmFIt} and \code{eBayes} in the
#'   \code{Limma} package. Differential regulated pathways were obtained by
#'   \code{gsva} function from \code{GSVA} package.
#'
#' @param druglist drug and its corresponding "drugbank ID" in
#'   \code{data.frame}
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
#'   should be taken as input here
#' @param core numeric, the number of core used for parallel running \code{gsva}
#' @param load_dir path to load or save modeling data files
#' @import prada
#' @import rhdf5
#' @import GSVA
#' @import limma
#' @importFrom stats model.matrix p.adjust
#'
#' @examples
#' druglist = data.frame(Druglist = c("Sorafenib", "Sunitinib"))
#' dataset = "92742" # which dataset to use
#' cellline = "HEPG2"
#' t = 6
#' \dontrun{
#' # access the number of available cores to run GSVA
#' library(parallel)
#' num_cores<-detectCores(logical=F)
#' }
#' num_cores = 2
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' \dontrun{
#' DEG_DEP_preparation(druglist = druglist, dataset = dataset, cellline =
#' cellline, treatment_time = t, core = num_cores, load_dir = load_dir)
#' }
#'


DEG_DEP_preparation <- function(druglist,
                                dataset,
                                cellline,
                                treatment_time,
                                core,
                                load_dir){

  if (requireNamespace("prada", "rhdf5", "GSVA", "limma","Hmisc")){

    source(paste0(load_dir, "/LINCS_data/l1ktools-master/R/cmap/io.R"))
    options(stringsAsFactors = F)
    options(quote=NULL)
    # druglist = read.csv("druglist.csv")
    druglist = druglist[1]
    names(druglist) = "Druglist"
    if(dataset =="70138"){
      inst = read.csv(paste0(load_dir, "/LINCS_data/GSE", dataset, "/file/GSE", dataset, "_Broad_LINCS_inst_info.txt"),header=T, sep="\t")

      drug_inst = inst[tolower(inst$pert_iname) %in% tolower(c(druglist$Druglist,"DMSO")),]

      geneinfo = read.csv(paste0(load_dir, "/LINCS_data/GSE", dataset, "/file/GSE92742_Broad_LINCS_gene_info.txt"),sep="\t",header = T)
      load(paste0(load_dir, "/LINCS_data/GeneSets_kegg_go.RData"))


      dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3"))

      drug_inst_selected = drug_inst[drug_inst$pert_time == treatment_time & drug_inst$cell_id == cell,]
      drug_inst_selected2 <- drug_inst_selected[drug_inst_selected$pert_dose >= 5 ,]

      if(nrow(drug_inst_selected2)!=0){

        dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline))
        dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time))
        dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time,"/all_geneset/"))
        dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time,"/all_genes/"))

        DMSO <- drug_inst[drug_inst$pert_time == treatment_time & drug_inst$cell_id == cellline & drug_inst$pert_dose==-666,]

        for(drug in unique( drug_inst_selected2$pert_iname)){

          if(drug != "DMSO"){

            # print(paste(cellline,treatment_time,drug))

            d = drug_inst_selected2[drug_inst_selected2$pert_iname == drug,]
            DMSO_match = DMSO[DMSO$det_plate %in% d$det_plate,]

            treatment_design = data.frame(treatment = c(rep(1,length(d$inst_id)),
                                                        rep(0,length(DMSO_match$inst_id))),
                                          treatment_ID =  c(d$inst_id,
                                                            DMSO_match$inst_id ),
                                          treatment_dose = c(d$pert_dose,
                                                             DMSO_match$pert_dose))

            row.names(treatment_design) = treatment_design$treatment_ID

            drug_70138_3h_L3 <- parse.gctx(paste0(load_dir, "/LINCS_data/GSE", dataset, "/GSE", dataset, "_Broad_LINCS_Level3_INF_mlr12k_n345976x12328.gctx"),
                                           cid = treatment_design$treatment_ID )
            mat <- as.data.frame(drug_70138_3h_L3@mat)
            mat$pr_gene_id <- rownames(mat)
            mat <- merge(mat, geneinfo, by = "pr_gene_id")
            rownames(mat) <- mat$pr_gene_symbol
            mat <- as.matrix(mat[c(treatment_design$treatment_ID)])
            design = cbind(data.frame(drug = rep(1, nrow(treatment_design))),treatment_design[1] )

            gsva_es <- GSVA::gsva(mat, gs.kegg, mx.diff = 1,parallel.sz = core)
            ## fit the same linear model now to the GSVA enrichment scores
            gsvafit <- limma::lmFit(gsva_es, design)
            ## estimate moderated t-statistics
            gsvafit <- limma::eBayes(gsvafit)
            ## set1 is differentially expressed
            #topTable(gsvafit, coef="treatment")
            gsva_allGeneSets <- limma::topTable(gsvafit, coef="treatment", number=Inf)
            gsva_allGeneSets$P.Value.fdr = p.adjust(gsva_allGeneSets$P.Value, method = "fdr")


            fit <- limma::lmFit(mat , design)
            fit <- limma::eBayes(fit)
            #topTable(fit, coef="treatment")
            allgenes <- limma::topTable(fit, coef="treatment", number=Inf)

            write.csv(allgenes, paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time,"/all_genes/allgenes_",drug,".csv"), quote = F)
            write.csv(gsva_allGeneSets, paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time,"/all_geneset/gsva_GeneSets_",drug,".csv"), quote = F)
          }
        }
      }
    }else if( dataset == "92742" ){

      inst = read.csv(paste0(load_dir, "LINCS_data/GSE", dataset, "/file/GSE", dataset, "_Broad_LINCS_inst_info.txt"),header=T, sep="\t")
      drug_inst = inst[tolower(inst$pert_iname) %in% tolower(c(druglist$Druglist,"DMSO")),]

      geneinfo = read.csv(paste0(load_dir, "/LINCS_data/GSE70138/file/GSE92742_Broad_LINCS_gene_info.txt"), sep="\t", header = T)
      load(paste0(load_dir, "/LINCS_data/GeneSets_kegg_go.RData"))

      dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3"))

      drug_inst[drug_inst$pert_time == "24|4",]$pert_time = "24"
      drug_inst[drug_inst$pert_time == "6|4",]$pert_time = "6"

      drug_inst_selected = drug_inst[drug_inst$pert_time == treatment_time & drug_inst$cell_id == cell,]
      drug_inst_selected$pert_dose =  as.numeric(drug_inst_selected$pert_dose)
      drug_inst_selected2 <- drug_inst_selected[drug_inst_selected$pert_dose >= 5 ,]

      if(nrow(drug_inst_selected2)!=0){

        dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline))
        dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time))
        dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time,"/all_geneset/"))
        dir.create(paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time,"/all_genes/"))


        DMSO <- drug_inst[drug_inst$pert_time == treatment_time & drug_inst$cell_id == cell & drug_inst$pert_dose==-666,]

        for(drug in unique(drug_inst_selected2$pert_iname)){

          if(drug != "DMSO"){

            d = drug_inst_selected[drug_inst_selected$pert_iname == drug,]
            d = d[d$pert_dose>= 5,]
            DMSO_match = DMSO[DMSO$rna_plate %in% d$rna_plate,]
            treatment_design = data.frame(treatment = c(rep(1,length(d$distil_id)),
                                                        rep(0,length(DMSO_match$distil_id))),
                                          treatment_ID =  c(d$distil_id,
                                                            DMSO_match$distil_id ),
                                          treatment_dose = c(d$pert_dose,
                                                             DMSO_match$pert_dose))
            row.names(treatment_design) = treatment_design$treatment_ID

            if(length(unique(treatment_design$treatment))>1){

              drug_92742_3h_L3 <- parse.gctx(paste0(load_dir, "/LINCS_data/GSE", dataset, "/GSE", dataset, "_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx"),
                                             cid = treatment_design$treatment_ID )

              mat <- as.data.frame(drug_92742_3h_L3@mat)
              mat$pr_gene_id <- rownames(mat)
              mat <- merge(mat, geneinfo, by = "pr_gene_id")
              rownames(mat) <- mat$pr_gene_symbol
              mat <- as.matrix(mat[c(treatment_design$treatment_ID)])
              design = model.matrix(~factor(treatment_design$treatment),levels=c(1,2))
              colnames(design) = c("control","treatment")
              gsva_es <- GSVA::gsva(mat, gs.kegg, mx.diff=1, parallel.sz = core)

              ## fit the same linear model now to the GSVA enrichment scores
              gsvafit <- limma::lmFit(gsva_es, design)
              ## estimate moderated t-statistics
              gsvafit <- limma::eBayes(gsvafit)
              ## set1 is differentially expressed
              #topTable(gsvafit, coef="treatment")
              gsva_allGeneSets <- limma::topTable(gsvafit, coef="treatment", number=Inf)
              gsva_allGeneSets$P.Value.fdr = p.adjust(gsva_allGeneSets$P.Value, method = "fdr")


              fit <- limma::lmFit(mat , design)
              fit <- limma::eBayes(fit)
              #topTable(fit, coef="treatment")
              allgenes <- limma::topTable(fit, coef="treatment", number=Inf)


              write.csv(allgenes, paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time,"/all_genes/allgenes_",drug,".csv"), quote = F)
              write.csv(gsva_allGeneSets, paste0(load_dir, "/LINCS_data/gsva_DComboNet_GSE", dataset, "_L3/",cellline,"/",treatment_time,"/all_geneset/gsva_GeneSets_",drug,".csv"), quote = F)


            }
          }
        }
      }
    }
  }
}


#' Drug induced differentially expressed gene list
#' preparation
#' @description The function \code{drugDEG_preparation} is to prepare drug
#'   induced differentially expressed gene table and pathway table for level two
#'   model(\code{model = "L2"}). If drug or cancer cell line you are interested
#'   in are not included in our dataset but contained in LINCS datasets
#'   ("GSE70138" or "GSE92742"), function \code{DEG_DEP_preparation} provides
#'   integrated analysis tools to generate matching files.
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
#'   should be taken as input here
#' @param foldchange numeric, fold change of gene expression between
#'   drug-treated cancer cell line and DMSO treated cell line
#' @param pvalue numeric, p-values corresponding to the t-statistics
#' @param load_dir path to load or save modeling data files
#' @import Hmisc
#' @importFrom stats model.matrix p.adjust
#'
#' @examples
#' dataset = "92742" # which dataset to use
#' cellline = "HEPG2"
#' t = 6
#' FC = 0.5
#' p = 0.05
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' \dontrun{
#' drugDEG_preparation(cellline = cellline, dataset = dataset, treatment_time =
#' t, foldchange = FC, pvalue = p, load_dir = load_dir)
#' }
#'

drugDEG_preparation = function(cellline,
                               dataset,
                               treatment_time,
                               foldchange,
                               pvalue,
                               load_dir){

  if (requireNamespace("Hmisc")){

    celllines = list.files(paste0(load_dir, "/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/"))

    a = list.files(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"/",treatment_time,"/all_genes/"))
    path = paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"/",treatment_time,"/all_genes/",a)
    n=length(a)
    druglist = gsub("allgenes_","", gsub(".csv","",a))
    druglist = Hmisc::capitalize(druglist)

    dDEp = read.table(path[1], sep = ",", header = T, stringsAsFactors = F)
    # dDEp = dDEp[dDEp$P.Value <=0.05,]
    dDEp = dDEp[abs(dDEp$logFC)>=foldchange & dDEp$P.Value <=pvalue,]
    if(nrow(dDEp) !=0){
      dDEp = data.frame(drug = druglist[1], DEG = dDEp$X, logFC=dDEp$logFC)
    }

    for(i in 2:n){
      dDEp2 = try(  read.table(path[i], sep = ",", header = T, stringsAsFactors = F), silent = TRUE)
      dDEp2 = dDEp2[abs(dDEp2$logFC)>=foldchange & dDEp2$P.Value <=pvalue,]
      # dDEp2 = dDEp2[ dDEp2$P.Value <=0.05,]
      if(nrow(dDEp2) !=0){
        dDEp2 = data.frame(drug = druglist[i], DEG =dDEp2$X, logFC=dDEp2$logFC)
        dDEp = rbind(dDEp, dDEp2)
      }
    }

    write.csv(dDEp, paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEG_",treatment_time,"h.csv"),quote=F, row.names = F)

  }
}


#' Drug induced pathway activation change list preparation
#' @description The function \code{drugDEP_preparation} is to prepare drug
#'   induced differentially expressed gene table and pathway table for level two
#'   model(\code{model = "L2"}). If drug or cancer cell line you are interested
#'   in are not included in our dataset but contained in LINCS datasets
#'   ("GSE70138" or "GSE92742"), function \code{DEG_DEP_preparation} provides
#'   integrated analysis tools to generate matching files.
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
#'   should be taken as input here
#' @param foldchange numeric, fold change of gene expression between drug-treated
#'   cancer cell line and DMSO treated cell line
#' @param pvalue numeric, p-values corresponding to the t-statistics
#' @param load_dir path to load or save modeling data files
#' @import Hmisc
#' @importFrom stats model.matrix p.adjust
#'
#' @examples
#' dataset = "92742" # which dataset to use
#' cellline = "HEPG2"
#' t = 6
#' FC = 0.5
#' p = 0.05
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' \dontrun{
#' drugDEP_preparation(cellline = cellline, dataset = dataset, treatment_time =
#' t, foldchange = FC, pvalue = p, load_dir = load_dir)
#'}
#'

drugDEP_preparation = function(cellline,
                               dataset,
                               treatment_time,
                               foldchange,
                               pvalue,
                               load_dir){

  if (requireNamespace("Hmisc")){

    # genepathway = read.table(paste0(load_dir,"/pathway/KEGG_ID_gene.txt"),sep="\t",header = T,stringsAsFactors = F)
    genepathway = unique(genepathway[c(1,2)])
    celllines = list.files(paste0(load_dir,"/LINCS/pathway_enrich_gsva_GSE",dataset,"/"))



    a = list.files(paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"/",treatment_time,"/all_geneset/"))
    path = paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"/",treatment_time,"/all_geneset/",a)
    n=length(a)
    druglist = gsub("gsva_GeneSets_","",
                    gsub(".csv","",a))
    druglist = Hmisc::capitalize(druglist)

    dDEp = read.table(path[1], sep = ",", header = T, stringsAsFactors = F)
    if(nrow(dDEp)!=0 ){
      dDEp = dDEp[dDEp$P.Value <= pvalue,]
      dDEp = data.frame(drug = druglist[1], pathway_name = dDEp$X, logFC=dDEp$logFC)
    }


    for(i in 2:n){
      dDEp2 = try(  read.table(path[i], sep = ",", header = T, stringsAsFactors = F), silent = TRUE)
      dDEp2 = dDEp2[dDEp2$P.Value <= pvalue,]
      if(nrow(dDEp2)!=0 ){

        dDEp2 = data.frame(drug = druglist[i], pathway_name =dDEp2$X, logFC=dDEp2$logFC)
        dDEp = rbind(dDEp, dDEp2)

      }
    }

    dDEp2 = merge(dDEp, genepathway, by = "pathway_name")
    write.csv(dDEp2, paste0(load_dir,"/LINCS_data/pathway_enrich_gsva_GSE",dataset,"/",cellline,"_DEgene_",treatment_time,"h.csv"),quote=F, row.names = F)
  }
}

