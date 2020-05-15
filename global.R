library(shinydashboard)
library(shinycssloaders)
library(DT)
library(tidyverse)
library(dendextend)
library(visNetwork)
library(heatmaply)
library(Matrix)
library(shinyjs)
library(shinyBS)

load("data/cerinaDataFinal.rda")
rna_norm_counts <- t(t(rna_counts)/size_factors$rna)
mirna_norm_counts <- t(t(mirna_counts)/size_factors$mirna)
circ_norm_counts <- t(t(circ_counts)/size_factors$rna)
rna_cpm_counts <- t(t(rna_counts)/lib_sizes$rna)*1e06
mirna_cpm_counts <- t(t(mirna_counts)/lib_sizes$mirna)*1e06
circ_srpbm_counts <- t(t(circ_counts)/lib_sizes$rna)*1e09
tissue_colors <- gplots::col2hex(c("burlywood3","blue","brown2","chartreuse3","darkgoldenrod2",
                                   "darkorchid","darkorange3","hotpink","palegreen2","slateblue2","gray45"))
names(tissue_colors) <- levels(design$Tissue)

pathway_colors <- c('#e6194b','#3cb44b','#ffe119','#4363d8','#f58231',
                    '#911eb4','#46f0f0','#f032e6','#bcf60c','#fabebe', 
                    '#008080','#e6beff','#9a6324','#fffac8','#800000', 
                    '#aaffc3','#808000','#ffd8b1','#000075','#808080')

infoPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
      )
    ),
    tags$a(
      href = "#", class = "btn btn-mini", `data-toggle` = "popover",`data-html`="true",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
      
      tags$i(class="fa fa-info-circle")
    )
  )
}


# pathway analysis
pathwayAnalysis <- function(gene_list, pathways, fdr.adjusted = TRUE, pval.thresh = .05){
  result <- list()
  for(i in 1:length(pathways)){
    tab <- pvals <- genes <- list()
    n_gene_list <- length(gene_list)
    n_pathway <- length(pathways[[i]])
    n_common <- length(intersect(gene_list, pathways[[i]]))
    n_genome <- 20000
    tab <- matrix(c(n_common,n_pathway,n_gene_list-n_common,n_genome-n_pathway), ncol = 2, byrow = TRUE)
    pvals <- fisher.test(tab, alternative = "g")$p.value
    genes <- intersect(gene_list, pathways[[i]])
    if(length(genes) == 0){
      genes <- NA
    }
    result[[i]] <- data.frame(Function = names(pathways)[i], pvalue = unlist(pvals), Genes = "")
    result[[i]]$Genes <- as.character(result[[i]]$Genes)
    for(k in 1:nrow(result[[i]])){
      result[[i]]$Genes[k] <- list(genes)
    }
    result[[i]]$Count <- length(genes)
    result[[i]]$Function.Size <- n_pathway
  }
  pathway_results <- do.call(rbind, result)
  pathway_results$FDR <- 0
  pathway_results$FDR <- p.adjust(pathway_results$pvalue, "fdr")
  pathway_results <- pathway_results[order(pathway_results$pvalue, decreasing = FALSE),]
  pathway_results$Function <- as.character(pathway_results$Function)
  if(fdr.adjusted){
    pathway_results <- pathway_results[pathway_results$FDR <= pval.thresh,]
  } else{
    pathway_results <- pathway_results[pathway_results$pvalue <= pval.thresh,] 
  }
  pathway_results$Function <- as.character(pathway_results$Function)
  pathway_results$Genes <- as.character(pathway_results$Genes)
  if(length(which(pathway_results$Genes == "NA")) > 0){
    pathway_results$Count[pathway_results$Genes == "NA"] <- 0 
  }
  pathway_results <- pathway_results[,c("Function", "Genes", "Count", "Function.Size", "pvalue", "FDR")]
  return(pathway_results)
}


annotationIdentifier <- function(id, id.map, default.ids, id.group = NULL){
  if(is.null(id.group)){
    id.first <- id[1]
    if(length(grep("hsa-miR", id.first, fixed = TRUE)) > 0 | length(grep("hsa-mir", id.first, fixed = TRUE)) > 0 | 
       length(grep("hsa-let", id.first, fixed = TRUE)) > 0 | id.first %in% default.ids){
      id.converted <- id
    } else if(length(grep("hsa_circ", id.first)) > 0){
      id.converted <- id.map$circ_ids$circRNA[match(id, id.map$circ_ids$circBase_ID, nomatch = 0)]
    } else if(length(grep("hsa-",id.first, fixed = TRUE)) > 0){
      id.converted <- id.map$circ_ids$circRNA[match(id, id.map$circ_ids$circAtlas_ID, nomatch = 0)]
    } else if(length(grep("|", id.first, fixed = TRUE)) > 0){
      id.converted <- id.map$circ_ids$circRNA[match(id, id.map$circ_ids$ciri2_ID, nomatch = 0)]
    } else{
      id.converted <- id.map$rna_ids$gene_id[match(id, id.map$rna_ids$gene_symbol, nomatch = 0)]
    }
  } else if(id.group == "circ"){
    id.first <- id[1]
    if(id.first %in% default.ids){
      id.converted <- id
    } else if(length(grep("hsa_circ", id.first)) > 0){
      id.converted <- id.map$circRNA[match(id, id.map$circBase_ID, nomatch = 0)]
    } else if(length(grep("hsa-",id.first, fixed = TRUE)) > 0){
      id.converted <- id.map$circRNA[match(id, id.map$circAtlas_ID, nomatch = 0)]
    } else if(length(grep("|", id.first, fixed = TRUE)) > 0){
      id.converted <- id.map$circRNA[match(id, id.map$ciri2_ID, nomatch = 0)]
    } else{
      id.converted <- id.map$circRNA[match(id, id.map$gene_symbol, nomatch = 0)]
    }
  }
  return(id.converted)
}


circRnaFuncAnalysis <- function(circ.mat, mirna.mat, circ.id, mirnas.selected = NULL, nrep = 10000, seed = 1){
  withProgress(message = 'Running permutation test...', min = 0, max = 1, value = .4, {
    set.seed(seed)
    circ.scores <- circ.mat[,circ.id]
    if(!is.null(mirnas.selected)){
      circ.scores[!names(circ.scores) %in% mirnas.selected] <- 0
    }
    mirnas <- names(circ.scores)[circ.scores > 0]
    rownames(circ.scores) <- NULL
    scores.observed <- mirna.mat %*% circ.scores
    n <- length(circ.scores)
    empty.list <- vector("list", nrep)
    circ.scores.list <- lapply(empty.list, function(x) sapply(n, function(x) sample(circ.scores, x)))
    circ.scores.mat <- Matrix(do.call(cbind, circ.scores.list), sparse = TRUE)
    scores.null <- mirna.mat %*% circ.scores.mat
    p <- (Matrix::rowSums(scores.null-scores.observed[,1] >= 0)+1)/(nrep+1)
    scores <- data.frame(Score = cume_dist(-log10(p)))
  })
  return(list(p = p, mirnas = mirnas, scores = scores))
}

cleanFun <- function(htmlString) {
  return(gsub("<.*?>", "", htmlString))
}

