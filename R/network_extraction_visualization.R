#' Extract predicted drug combination subnetwork
#'
#' @description Function \code{network_extract} is to extract subnetwork of
#'   predicted drug combinations from drug-gene/pathway network to help infer
#'   possible mechanism of drug combinations and visualization. Note that this
#'   function only works after user did prediction via \code{DComboNet()} and
#'   drug/gene/pathway rank tables have been generated and saved in fetchable
#'   path provided in \code{resultdir}.
#' @param drugseed drug that user is interested in perdicting combinable drug(s)
#' @param drugcandidate predicted combinable drug for drugseed
#' @param drugtarget drug-target gene interaction file corresponding to
#'   \code{drugseed} and \code{drugcandidate} in \code{data.frame} format
#'   (otherwise \code{as.data.frame} can be used to transform the format of
#'   file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol).
#' @param generank_filter numeric, to select genes according to the rank,
#'   default value is 0.01 meaning only genes ranked on top 1% will be kept in
#'   subnetwork
#' @param pathwayrank_filter umeric, to select pathways according to the rank,
#'   default value is 0.1 meaning only pathways ranked on top 10% will be kept in
#'   subnetwork
#' @param model choose model type, two models here to choose: "L1", "L2"
#' @param cellline cancer cell lines name.
#' @param load_dir path to load or save modelin data files
#' @param resultdir path to save result files
#' @import igraph
#' @return Predicted drug combination subnetwork table extracted from
#'   drug-gene/pathway network constructed based on user requirements
#'
#'
#' @examples
#' drugseed = "Sorafenib"
#' drugcandidate = "Vorinostat"
#' drugtarget = data.frame(Drug = c(rep(drugseed,10),rep(drugcandidate,5)),target
#' = c("BRAF", "FGFR1", "FLT1", "FLT3", "FLT4", "KDR", "KIT", "PDGFRB", "RAF1",
#' "RET", "HDAC1", "HDAC2", "HDAC3", "HDAC6", "HDAC8"))
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' resultdir = "G:/lab/DCcomboNet/Rpackage/tryout_result/"
#' \dontrun{network_extract(drugseed = drugseed, drugcandidate = drugcandidate,
#' drugtarget = drugtarget, generank_filter = 0.01, pathwayrank_filter = 0.1,
#' model = "L1", load_dir = load_dir, resultdir = resultdir)}
#' @export
#'
network_extract <- function(drugseed,
                            drugcandidate,
                            drugtarget = NULL,
                            generank_filter = 0.01,
                            pathwayrank_filter = 0.1,
                            model = c('L1','L2'),
                            cellline = NULL,
                            load_dir,
                            resultdir){

  if (requireNamespace('igraph','load_all')){
    # dt_lib = read.csv(paste0(load_dir, '/data/drugtarget.csv'), sep = ',',header = T, stringsAsFactors = F)
    devtools::load_all()

    if(is.null(drugtarget)){
      drugtarget = dt_lib
    }else{
      names(drugtarget) = names(dt_lib)
      drugtarget = unique(rbind(dt_lib, drugtarget))
    }

    if(model == 'L1'){
      load(paste0(resultdir,model,'_result/potential_net/net_data.Rdata'))
    }else if(model == 'L2'){
      load(paste0(resultdir,model,'_result/potential_net/',cellline,'_net_data.Rdata'))
    }
    names(genenet) = c('N1','N2')
    names(genepathwaynet) = c('N1','N2')
    multiplex_graph = igraph::graph_from_data_frame(multiplex_df)
    igraph::E(multiplex_graph)$weight = multiplex_df$value
    multiplex_graph <- igraph::simplify(multiplex_graph, remove.multiple = T, remove.loops = T)
    node_mn = igraph::V(multiplex_graph)
    igraph::V(multiplex_graph)$names = node_mn


    generank = read.csv(paste0(resultdir,model,'_result/generank/',drugseed,'_rank.csv'))
    generank$rank = 1:nrow(generank)
    pathwayrank = read.csv(paste0(resultdir,model,'_result/pathwayrank/',drugseed, '_rank.csv'))
    pathwayrank$rank = 1:nrow(pathwayrank)


    dp_subN1 = dnet::dNetInduce(multiplex_graph, nodes_query=c(drugseed,drugcandidate), knn=1, remove.loops=T, largest.comp=T)
    adj1 = igraph::get.adjacency(dp_subN1)
    df1 = reshape2::melt(as.matrix(adj1))
    df1 = df1[df1$value!=0,]
    df1 = df1[df1$Var1%in%c(drugseed, drugcandidate) | df1$Var2%in%c(drugseed,drugcandidate),]
    names(df1) = c('GeneNames','Var2','Value')
    df1_generank = merge(generank,df1,by='GeneNames')
    names(df1) = c('PathwayID','Var2','Value')
    df1_pathwayrank = merge(pathwayrank,df1,by='PathwayID')
    df1_generank2 = df1_generank[df1_generank$rank<=round( generank_filter * nrow(generank),0),]
    df1_pathwayrank2 = df1_pathwayrank[df1_pathwayrank$rank<=round(pathwayrank_filter * nrow(df1_pathwayrank),0),]
    df1_generank_f = df1_generank2[c(4,1)];names(df1_generank_f) = c('N1','N2')
    dt =  drugtarget[drugtarget$Drug %in% c(drugseed, drugcandidate),];names(dt ) = c('N1','N2')
    df1_generank_f = rbind(df1_generank_f,dt)
    df1_pathwayrank_f = df1_pathwayrank2[c(4,1)];names(df1_pathwayrank_f) = c('N1','N2')

    df1_gene=unique(df1_generank_f[2])
    gene_connection = merge(df1_gene,genenet,by='N2')
    names(df1_gene)='N1'
    gene_connection = merge(df1_gene,gene_connection,by='N1')


    df1_pathway=unique(df1_pathwayrank_f[2])
    dp_connection = merge(df1_pathway,genepathwaynet,by='N2')
    names(df1_pathway)='N1'
    dp_connection = merge(df1_pathway,dp_connection,by='N1')


    gp_connection1 = merge(unique(gene_connection[1]),genepathwaynet,by='N1')
    nrow(gp_connection1)
    gp_connection2 = merge(unique(gene_connection[2]),genepathwaynet,by='N2')
    nrow(gp_connection2)

    gp_connection_final = gp_connection2[gp_connection2$N1%in%df1_pathwayrank_f$N2, ]

    names(gene_connection)=c('N1','N2')
    names(gp_connection_final)=c('N1','N2')

    df1_subN2 = unique(rbind(df1_generank_f,
                             df1_pathwayrank_f,
                             gene_connection[c('N1','N2')],
                             gp_connection_final[c('N1','N2')]))

    df1_generank2 = df1_generank2[c(4,1,2,3)];names(df1_generank2)=c('N1','N2','Gene.Score','Gene.rank')
    df1_pathwayrank2 = df1_pathwayrank2[c(4,1,2,3)];names(df1_pathwayrank2)=c('N1','N2','Pathway.Score','Pathway.rank')
    network_file = merge(df1_subN2,df1_generank2,by=c('N1','N2'),all.x = T)
    network_file = merge(network_file,df1_pathwayrank2,by=c('N1','N2'),all.x = T)

    dir.create(paste0(resultdir, model,'_result/infered_net/'))

    if(model == "L1"){
      write.csv(network_file,paste0(resultdir, model,'_result/infered_net/', drugseed, '_', drugcandidate,'_network_file.csv'),quote = F,row.names = F)
    }else if(model == "L2"){
      write.csv(network_file,paste0(resultdir, model,'_result/infered_net/',cellline,'_', drugseed, '_', drugcandidate,'_network_file.csv'),quote = F,row.names = F)

    }
  }
}


#' Predicted drug combination subnetwork visualization
#'
#' @description Function \code{network_visualization} is to visualize subnetwork
#'   of predicted drug combinations from drug-gene/pathway network to help infer
#'   possible mechanism of drug combinations. Note that subnetwork table should
#'   be prepared via function \code{network_extract}.
#' @param drugseed drug that user is interested in perdicting combinable drug(s)
#' @param drugcandidate predicted combinable drug for drugseed
#' @param drugtarget drug-target gene interaction file in \code{data.frame}
#'   format (otherwise \code{as.data.frame} can be used to transform the format
#'   of file). Target gene names should convert to official gene name (i.e. HGNC
#'   gene symbol)
#' @param model choose model type, two models here to choose: "L1", "L2"
#' @param cellline cancer cell lines name.
#' @param load_dir path to load or save modelin data files
#' @param resultdir path to save result files
#' @return Interactive subnetwork graph of predicted drug combination. An
#'   \code{.graphml} network file will be saved in provided path and can be
#'   easily import into \code{Cytoscape} for further analysis.
#'
#' @examples
#' drugseed = "Sorafenib"
#' drugcandidate = "Vorinostat"
#' drugtarget = data.frame(Drug = c(rep(drugseed,10),rep(drugcandidate,5)),target
#' = c("BRAF", "FGFR1", "FLT1", "FLT3", "FLT4", "KDR", "KIT", "PDGFRB", "RAF1",
#' "RET", "HDAC1", "HDAC2", "HDAC3", "HDAC6", "HDAC8"))
#'
#' load_dir = "G:/lab/DCcomboNet/Rpackage/input_data/"
#' resultdir = "G:/lab/DCcomboNet/Rpackage/tryout_result/"
#'
#' \dontrun{network_visualization(drugseed = drugseed, drugcandidate =
#' drugcandidate, drugtarget = drugtarget, model = "L2", load_dir = load_dir,
#' resultdir = resultdir)}
#'
#' @export
#'

network_visualization <- function(drugseed,
                                  drugcandidate,
                                  drugtarget = NULL,
                                  model = c('L1','L2'),
                                  cellline = NULL,
                                  load_dir,
                                  resultdir){

  if (requireNamespace('devtools','igraph','visNetwork')){

    # dt_lib = read.csv(paste0(load_dir, '/data/drugtarget.csv'), sep = ',',header = T, stringsAsFactors = F)

    devtools::load_all()

    if(is.null(drugtarget)){
      drugtarget = dt_lib
    }else{
      names(drugtarget) = names(dt_lib)
      drugtarget = unique(rbind(dt_lib, drugtarget))
    }

    if(model == "L1"){
      df = read.table(paste0(resultdir, model,'_result/infered_net/', drugseed, '_', drugcandidate,'_network_file.csv'), sep = ',', header = T, stringsAsFactors = F)
    }else if(model == "L2"){
      df = read.table(paste0(resultdir, model,'_result/infered_net/',cellline,'_', drugseed, '_', drugcandidate,'_network_file.csv'), sep = ',', header = T, stringsAsFactors = F)
    }

    # df <- df[order(df$combined_score, decreasing = T),]
    nodes = union(df[,1],df[,2])
    net = igraph::graph_from_data_frame(d=df, vertices = unique(nodes), directed = F)
    net = igraph::simplify(net, remove.multiple = TRUE, remove.loops = TRUE)

    edge.color <- colorRampPalette(c('#D6D6D6','#383838'), alpha=TRUE)
    igraph::E(net)$color <- edge.color(igraph::ecount(net))

    drugseed_target = drugtarget[drugtarget$Drug == drugseed, ]$Target
    drugcandidate_target = drugtarget[drugtarget$Drug == drugcandidate, ]$Target

    igraph::V(net)$color = '#E7F3FD'
    igraph::V(net)[which(igraph::V(net)$name == drugseed)]$color = '#89D0F5'
    igraph::V(net)[which(igraph::V(net)$name == drugcandidate)]$color = '#FF9900'
    igraph::V(net)[which(igraph::V(net)$name %in% drugseed_target)]$color = '#99FFFF'
    igraph::V(net)[which(igraph::V(net)$name %in% drugcandidate_target)]$color = '#FFCC99'
    igraph::V(net)$shape = 'box'
    igraph::V(net)[which(igraph::V(net)$name %in% c(drugseed, drugcandidate))]$shape = 'ellipse'
    igraph::V(net)$group = 'other-genes/pathways'
    igraph::V(net)[which(igraph::V(net)$name %in% drugseed)]$group = 'drugseed'
    igraph::V(net)[which(igraph::V(net)$name %in% drugcandidate)]$group = 'drugcandidate'
    igraph::V(net)[which(igraph::V(net)$name %in% drugtarget[drugtarget$Drug %in% drugseed,]$Target)]$group = 'drugseed-target-genes'
    igraph::V(net)[which(igraph::V(net)$name %in% drugtarget[drugtarget$Drug %in% drugcandidate,]$Target)]$group = 'drugcandidate-target-genes'
    igraph::V(net)$font.size = 40

    # tkid <- igraph::tkplot(net) #tkid is the id of the tkplot that will open
    # l <- tkplot.getcoords(18) # grab the coordinates from tkplot
    igraph::write_graph(net, paste0(resultdir, model,'_result/infered_net/', drugseed, '_', drugcandidate,'subnetwork.graphml'), format = "graphml")


    data <- toVisNetworkData(net)

    visNetwork(nodes = data$nodes, edges = data$edges)%>%
      visNodes(size = 40)%>%
      # visHierarchicalLayout()%>%
      visIgraphLayout(type = 'full',randomSeed = 123) %>%
      visOptions(highlightNearest = list(enabled = TRUE,  hideColor = 'lightgrey', hover = T),
                 nodesIdSelection =list(enabled = TRUE), selectedBy = "group") %>%
      # visConfigure(enabled = TRUE)#%>%
      addFontAwesome() %>%
      visGroups(groupname = "drugseed", color = '#89D0F5')%>%
      visGroups(groupname = "drugcandidate", color = '#FF9900')%>%
      visGroups(groupname = "drugseed-target-genes", color = '#99FFFF', shape = 'box')%>%
      visGroups(groupname = "drugcandidate-target-genes", color = '#FFCC99', shape = 'box')%>%
      visGroups(groupname = "other-genes/pathways", color = '#E7F3FD', shape = 'box')%>%
      visLegend() %>%
      visInteraction(navigationButtons = TRUE) %>%
      visOptions(manipulation = TRUE)  %>%
      visSave(file = paste0(resultdir, model,'_result/infered_net/network.html'))


    visNetwork(nodes = data$nodes, edges = data$edges)%>%
      visNodes(size = 40)%>%
      visIgraphLayout(type = 'full',randomSeed = 123) %>%
      visOptions(highlightNearest = list(enabled = TRUE,  hideColor = 'lightgrey', hover = T),
                 nodesIdSelection = TRUE, selectedBy = "group") %>%
      # visConfigure(enabled = TRUE)#%>%
      addFontAwesome() %>%
      visGroups(groupname = "drugseed", color = '#89D0F5')%>%
      visGroups(groupname = "drugcandidate", color = '#FF9900')%>%
      visGroups(groupname = "drugseed-target-genes", color = '#99FFFF', shape = 'box')%>%
      visGroups(groupname = "drugcandidate-target-genes", color = '#FFCC99', shape = 'box')%>%
      visGroups(groupname = "other-genes/pathways", color = '#E7F3FD', shape = 'box')%>%
      visLegend() %>%
      visInteraction(navigationButtons = TRUE) #%>%
    # visOptions(manipulation = TRUE) # %>%
    # visSave(file = paste0(resultdir, model,'_result/infered_net/network.html'))


  }
}
