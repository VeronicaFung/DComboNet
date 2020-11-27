#' druglistlib
#'
#' all drugs in pre-constructed network
#' @docType data
#' @name druglistlib
NULL

#' drugpair_lib
#'
#' all drug pairs in pre-constructed network, including known drug combinations
#' marked as \code{P} and randomly reshuffled pairs marked as \code{N} in "TAG"
#' column.
#' @docType data
#' @name drugpair_lib
NULL

#' Hierarchical network of ATC code system
#'
#' ATC network, igraph object, used to calculate drug ATC based similarity.
#' @docType data
#' @name ATC

NULL

#' Nodes in ATC network
#'
#' nodes in ATC network, used to calculate drug ATC based similarity.
#' @docType data
#' @name nodes

NULL

#' all_drug_atc
#'
#' ATC codes in ATC network and their drug names
#' @docType data
#' @name all_drug_atc
NULL

#' CDK fingerprints
#'
#' CDK fingerprint table for drugs in pre-constructed network generated via
#' PaDel software
#' @docType data
#' @name CDK_FP_lib
NULL

#' MACCS fingerprints
#'
#' MACCS fingerprint table for drugs in pre-constructed network generated via
#' PaDel software
#' @docType data
#' @name MACCS_FP_lib
NULL

#' PubChem fingerprints
#'
#' PubChem fingerprint table for drugs in pre-constructed network generated via
#' PaDel software
#' @docType data
#' @name pubchemFP_lib
NULL


#' Drug side effect database
#'
#' all drugs and their corresponding side effect terms extracted from SIDER(v4)
#' database
#' @docType data
#' @name all_side_effect
NULL

#' Drugs and their side effect terms
#'
#' all drugs in pre-constructed network and their corresponding side effect
#' terms extracted from SIDER(v4) database
#' @docType data
#' @name drugselib
NULL

#' Drugs and their target genes
#'
#' Drugs in pre-constructed network and their target genes
#' @docType data
#' @name dt_lib
NULL

#' Drugs and their related genes
#'
#' Drugs in pre-constructed network and their related genes extracted from IPA
#' (Ingenuity® Pathway Analysis) tool.
#' @docType data
#' @name dg_lib
NULL


#' Protein-protein interaction network
#'
#' Protein-protein interaction network table \code{PPI_table}, extracted from inBio map. This
#' table contains protein-protein interactions and the confidence scores
#' supporting the existence of edges.
#' @docType data
#' @name PPI_table
NULL


#' Cancer related genes
#'
#' Cancer related gene list, extracted from KEGG Cancer related pathways
#' @docType data
#' @name cancer_gene
NULL

#' Cancer cell lines in CCLE database
#'
#' Avaliable cancer cell lines in CCLE database, checking if user inputs cell
#' line name in including in CCLE so Level 2 model can be called without
#' requiring extra cancer specific expressed gene table
#' @docType data
#' @name cellline_lib
NULL


#' Cancer cell lines in LINCS database
#'
#' Avaliable cancer cell lines in LINCS database, checking if user inputs cell
#' line name in including in LINCS so Level 2 model can be called without
#' requiring extra drug induced differentially expressed gene table
#' @docType data
#' @name celllineDEG_lib
NULL



#' KEGG pathway and genes within pathways
#'
#' Table of KEGG pathway and genes within pathways. The table contains three
#' columns, first column denotes the KEGG pathway ID, second column denotes the
#' names of these pathways (all letter in uppercase), third column denotes genes
#' within these pathways.
#' @docType data
#' @name genepathway
NULL

#' KEGG pathway-pathway interaction network
#'
#' Table of KEGG pathway-pathway interaction network downloaded from publication
#' (PMID: 26687590). The first two columns denote the interaction between two
#' pathway nodes; the third column denotes the interaction types.
#' @docType data
#' @name WWI_all
NULL

#' Cancer cell lines in CCLE database
#'
#' Avaliable cancer cell lines in CCLE database, checking if user inputs cell
#' line name in including in CCLE so Level 2 model can be called without extra
#' drug induced differentially activated pathway table inputed
#' @docType data
#' @name celllineDEP_lib
NULL

