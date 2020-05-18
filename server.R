options(shiny.maxRequestSize=1000*1024^2)
server <- function(input, output, session) {
  plotDendrogramObject <- reactive({
    p <- vector("list", 3)
    for(i in 1:length(p)){
      p[[i]] <- ggplot(as.ggdend(hc[[i]] %>% set('labels_cex', .9))) + 
        ylim(-.3*max(get_branches_heights(hc[[i]])), max(get_branches_heights(hc[[i]])))
    }
    names(p) <- names(hc)
    return(p)
  })
  
  output$rnaDendrogram <- renderPlot({
    plotDendrogramObject()$rna
  }, height = 600)
  
  output$mirnaDendrogram <- renderPlot({
    plotDendrogramObject()$mirna
  }, height = 600)
  
  output$circDendrogram <- renderPlot({
    plotDendrogramObject()$circ
  }, height = 600)
  
  circAnnots <- reactive({
    circAnnots <- annots[annots$circRNA %in% rownames(circ_norm_counts),]
    circAnnots <- c(as.character(circAnnots$circRNA), as.character(circAnnots$ciri2_ID), 
                    as.character(circAnnots$circAtlas_ID), as.character(circAnnots$circBase_ID), 
                    as.character(circAnnots$gene_symbol[-which(circAnnots$gene_symbol == "n/a")]))
    circAnnots <- unique(circAnnots[-which(is.na(circAnnots))])
    return(circAnnots)
  })
  
  correlationAxisVariables <- reactive({
    circAnnots <- circAnnots()[-which(circAnnots() %in% annots$gene_symbol)]
    vars <- unique(c(circAnnots, rownames(mirna_norm_counts), rownames(rna_norm_counts), as.character(rna.norm.summary$gene_symbol)))
    vars
  })
  
  circ_current_selection <- reactiveVal(NULL)
  
  # now store your current selection in the reactive value
  observeEvent(input$selectCirc, {
    circ_current_selection(input$selectCirc)
  })
  
  observe({
    updateSelectizeInput(session, "selectGene", choices = c(rownames(rna_counts), unique(as.character(rna.norm.summary$gene_symbol))), 
                         selected = "PIK3R3", options = list(maxItems = 2), server = TRUE)
    updateSelectizeInput(session, "selectMiRna", choices = c(rownames(mirna_counts)), 
                         selected = "hsa-miR-7-5p", options = list(maxItems = 2), server = TRUE)
    updateSelectizeInput(session, "selectCircRna", choices = circAnnots(), 
                         selected = "CDR1-AS", options = list(maxItems = 2), server = TRUE)
    updateSelectizeInput(session, "xaxis", choices = correlationAxisVariables(), 
                         selected = "chrX:139865339-139866824", options = list(maxItems = 2), server = TRUE)
    updateSelectizeInput(session, "yaxis", choices = correlationAxisVariables(), 
                         selected = "PIK3R3", options = list(maxItems = 2), server = TRUE)
    updateSelectizeInput(session, "selectMirnaNetwork", choices = unique(scores_circ_mirna$miRNA), 
                         selected = "hsa-miR-7-5p", options = list(maxItems = 2), server = TRUE)
    updateSelectizeInput(session, "circForAnalysis", choices = circAnnots(), 
                         selected = "CDR1-AS", options = list(maxItems = 2), server = TRUE)
  })
  
  expressionProfile <- reactive({
    if(input$expressionProfileTabs == "RNA"){
      profile <- rna.norm.summary
      if(!input$rnaNormCounts){
        profile$Expression <- profile$Mean.CPM
        ylabel <- "CPM"
        if(input$rnaLog2Scale){
          profile$Expression <- profile$Mean.Log2.CPM
        }
      } else{
        profile$Expression <- profile$Mean.Expression
        ylabel <- "Normalized Expression"
        if(input$rnaLog2Scale){
          profile$Expression <- profile$Mean.Log2.Expression
        }
      }
      selectVar <- input$selectGene
      if(length(grep("ENSG", selectVar)) > 0){
        var <- sym("gene_id")
      } else {
        var <- sym("gene_symbol")
      }
    } else if(input$expressionProfileTabs == "miRNA"){
      profile <- mirna.norm.summary
      if(!input$mirnaNormCounts){
        profile$Expression <- profile$Mean.CPM
        ylabel <- "CPM"
        if(input$mirnaLog2Scale){
          profile$Expression <- profile$Mean.Log2.CPM
        }
      } else{
        profile$Expression <- profile$Mean.Expression
        ylabel <- "Normalized Expression"
        if(input$mirnaLog2Scale){
          profile$Expression <- profile$Mean.Log2.Expression
        }
      }
      selectVar <- input$selectMiRna
      var <- sym("miRNA")
    } else {
      profile <- circ.norm.summary %>% left_join(annots, by = "circRNA") %>%
        select(Tissue, circRNA, ciri2_ID, circAtlas_ID, circBase_ID, gene_symbol, Mean.Expression, Mean.Log2.Expression, Mean.SRPBM, Mean.Log2.SRPBM)
      if(!input$circNormCounts){
        profile$Expression <- profile$Mean.SRPBM
        ylabel <- "SRPBM"
        if(input$circLog2Scale){
          profile$Expression <- profile$Mean.Log2.SRPBM
        }
      } else{
        profile$Expression <- profile$Mean.Expression
        ylabel <- "Normalized Expression"
        if(input$circLog2Scale){
          profile$Expression <- profile$Mean.Log2.Expression
        }
      }
      selectVar <- input$selectCircRna
      if(length(grep("hsa_circ", selectVar)) > 0){
        var <- sym("circBase_ID")
      } else if(length(grep("hsa-", selectVar, fixed = TRUE)) > 0) {
        var <- sym("circAtlas_ID")
      } else if(length(grep("|", selectVar, fixed = TRUE)) > 0){
        var <- sym("ciri2_ID")
      } else if(length(grep("chr", selectVar)) > 0){
        var <- sym("circRNA")
      } else{
        var <- sym("gene_symbol")
      }
    }
    if(selectVar == ""){return(NULL)}
    return(list(profile = profile, selectVar = selectVar, var = var, ylabel = ylabel))
  })
  
  expressionProfileData <- reactive({
    if(is.null(expressionProfile())){return(NULL)}
    dat <- filter(expressionProfile()$profile, !!expressionProfile()$var == expressionProfile()$selectVar)
    return(dat)
  })
  
  output$multipleMapping <- renderUI({
    if(length(unique(expressionProfileData()$circRNA)) > 0 & input$selectCircRna %in% annots$gene_symbol){
      return(
        selectizeInput("selectCirc", "circRNA by Position", choices = unique(expressionProfileData()$circRNA), selected = circ_current_selection(), multiple = FALSE)
      )
    } else{
      return(NULL)
    }
  })
  
  expressionProfileObject <- reactive({
    if(is.null(expressionProfileData())){return(NULL)}
    ylabel <- expressionProfile()$ylabel
    if(isTruthy(input$selectCirc) & length(unique(expressionProfileData()$circRNA)) > 1){
      p <- ggplot(filter(expressionProfile()$profile, !!sym("circRNA") == input$selectCirc), aes(x = Tissue, y = Expression)) + 
        stat_summary(geom = "bar", aes(fill = Tissue), width = .4) + 
        scale_fill_manual(values = tissue_colors) + ylab(ylabel) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
              legend.position = "none")
    } else{
      p <- ggplot(expressionProfileData(), aes(x = Tissue, y = Expression)) + 
        stat_summary(geom = "bar", aes(fill = Tissue), width = .4) + 
        scale_fill_manual(values = tissue_colors) + ylab(ylabel) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
              legend.position = "none") 
    }
    return(p)
  })
  
  output$rnaExpressionProfilePlot <- renderPlot({
    expressionProfileObject()
  }, height = 500)
  
  output$mirnaExpressionProfilePlot <- renderPlot({
    expressionProfileObject()
  }, height = 500)
  
  output$circExpressionProfilePlot <- renderPlot({
    expressionProfileObject()
  }, height = 500)
  
  correlationData <- reactive({
    if(!input$normCounts){
      rna_norm_counts <- rna_cpm_counts
      mirna_norm_counts <- mirna_cpm_counts
      circ_norm_counts <- circ_srpbm_counts
    }
    rna_ids <- rna.norm.summary[match(rownames(rna_norm_counts), rna.norm.summary$gene_id, nomatch = 0),c("gene_id","gene_symbol")]
    mirna_ids <- mirna.norm.summary[match(rownames(mirna_norm_counts), mirna.norm.summary$miRNA, nomatch = 0),c("miRNA","miRNA")]
    circ_ids <- annots[match(rownames(circ_norm_counts), annots$circRNA, nomatch = 0),c("circRNA","ciri2_ID","circAtlas_ID","circBase_ID")]
    norm_data <- data.frame(t(rbind(log2(rna_norm_counts + .1), log2(mirna_norm_counts + .1), log2(circ_norm_counts + .1))), check.names = FALSE)
    norm_data$Tissue <- design$Tissue
    ids <- list(rna_ids = rna_ids, mirna_ids = mirna_ids, circ_ids = circ_ids)
    return(list(data = norm_data, ids = ids))
  })
  
  correlationPlot <- reactive({
    ids <- correlationData()$ids
    dat <- correlationData()$data
    colnames(dat)[-ncol(dat)] <- make.unique(colnames(dat)[-ncol(dat)])
    xaxis <- input$xaxis
    yaxis <- input$yaxis
    x_name <- annotationIdentifier(id = xaxis, id.map = ids, default.ids = colnames(dat))
    y_name <- annotationIdentifier(id = yaxis, id.map = ids, default.ids = colnames(dat))
    corr <- c(pearson = round(cor(dat[,x_name], dat[,y_name]), 2), 
              spearman = round(cor(dat[,x_name], dat[,y_name], method = "spearman"), 2))
    p <- ggplot(dat, aes_(as.name(x_name), as.name(y_name))) + geom_point(aes(colour = Tissue), size = 5) + 
      scale_color_manual(values = tissue_colors) + xlab(input$xaxis) + ylab(input$yaxis)
    p <- p + geom_smooth(method = lm, se = TRUE) + 
      ggtitle(paste0("pearson r = ", corr["pearson"], ", spearman rho = ", corr["spearman"])) +
      theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12), legend.title = element_text(size = 14), 
            legend.text = element_text(size = 12), title = element_text(size = 14))
    return(list(p = p, corr = corr))
  })
  
  output$correlationScatterPlot <- renderPlot({
    correlationPlot()$p
  }, height = 700)
  
  tissueSpecificSelection <- reactiveVal(NULL)
  
  circRnas <- reactive({
    if(is.null(tissueSpecificSelection())){
      return(NULL)
    } else if(tissueSpecificSelection() == "unspecific"){
      choices <- unique(scores_circ_mirna$circRNA[!scores_circ_mirna$circRNA %in% 
                                                    tissue_specific_summary$circRNA[tissue_specific_summary$Tissue == tissueSpecificSelection()]])
    } else {
      choices <- unique(scores_circ_mirna$circRNA[scores_circ_mirna$circRNA %in% 
                                                    tissue_specific_summary$circRNA[tissue_specific_summary$Tissue == tissueSpecificSelection()]])
    }
    return(choices)
  })
  
  output$correlationThresholds <- renderUI({
    if(input$spearmanFilter){
      return(
        numericInput("mRnaCircRnaCorrThresh", "circRNA-mRNA spearman rho",
                     min = 0, max = 1, value = .5, step = .05)
      )
    } else{
      return(NULL)
    }
  })
  
  mirnaNetworkData <- eventReactive(input$plotMiRnaNetwork,{
    if(!input$mirnaNetworkCircNormCounts){
      circ_norm_counts <- circ_srpbm_counts
      circ.norm.summary$Mean.Expression <- circ.norm.summary$Mean.SRPBM
    }
    circs <- unique(circ.norm.summary$circRNA[circ.norm.summary$Mean.Expression >= input$filterByExpression])
    edges <- scores_circ_mirna[scores_circ_mirna$miRNA %in% input$selectMirnaNetwork,] %>%
      select(miRNA, circRNA, MRE, MRE.per.kb, Score) %>%
      rename(from = miRNA, to = circRNA) %>%
      filter(to %in% circs & Score >= input$filterByScore)
    corr_mat <- data.frame(cor(t(circ_norm_counts[rownames(circ_norm_counts) %in% edges$to,])), check.names = FALSE)
    corr_mat$circRNA <- rownames(corr_mat)
    corr_scores <- reshape2::melt(corr_mat)
    corr_scores <- corr_scores[!corr_scores$circRNA == corr_scores$variable,] %>%
      rename(from = circRNA, to = variable, Score = value) %>%
      filter(Score >= input$circRnaCorrThresh) %>%
      add_column(MRE = NA, MRE.per.kb = NA, .before = "Score")
    edges <- rbind(edges, corr_scores) %>%
      mutate(title = paste0("<p>", from, "-", to,"<br>","weight = ",Score), color = "grey")
    nodes <- tibble(id = unique(c(edges$from, edges$to))) %>%
      mutate(label = c(id[1], rep("", n()-1))) %>%
      left_join(tissue_specific_summary, by = c("id" = "circRNA")) %>%
      select(id, label, Tissue) %>%
      mutate(color = tissue_colors[match(Tissue, names(tissue_colors))]) %>%
      mutate(color = ifelse(label != "", "#f0e0ce", ifelse(is.na(Tissue), "#B4C3DE", color))) %>%
      mutate(size = ifelse(label != "", 70, 30)) %>%
      mutate(shape = ifelse(color == "#f0e0ce", "circle", "dot")) %>%
      mutate(font.size = ifelse(length(unique(id)) > 30, 30, length(unique(id)))) %>%
      rename(group = Tissue) %>%
      mutate(group = ifelse(color == "#f0e0ce", "", ifelse(is.na(group), "unspecific", as.character(group)))) %>%
      mutate(title = paste0("<p>",id,"<br>","<strong>circBase ID:</strong> ", annots$circBase_ID[match(id, annots$circRNA)], "<br>","<strong>circAtlas ID:</strong> ",
                            annots$circAtlas_ID[match(id, annots$circRNA)], "<br>", "<strong>ciri2 ID:</strong> ", annots$ciri2_ID[match(id, annots$circRNA)],
                            "<br>","<strong>Gene Symbol:</strong> ", annots$gene_symbol[match(id, annots$circRNA)]))
    nodes$title[1] <- paste0("<p>",nodes$id[1],"</p>")
    vis_obj <- list(nodes = nodes, edges = edges)
    tissue_colors <- c(tissue_colors, unspecific = "#B4C3DE")
    labels <- unique(vis_obj$nodes$group)[unique(vis_obj$nodes$group) != ""]
    legend_dat <- data.frame(label = labels,
                             color.background = c(tissue_colors[match(labels,names(tissue_colors), nomatch = 0)]),
                             color.border = c(tissue_colors[match(labels,names(tissue_colors), nomatch = 0)]))
    legend_dat <- data.frame(rbind(legend_dat[legend_dat$label != "unspecific",], 
                                   legend_dat[legend_dat$label == "unspecific",]))
    legend_dat$shape <- "dot"
    legend_dat$size <- 10
    return(list(network_obj = vis_obj, legend_obj = legend_dat))
  })

 
  output$mirnaNetworkPlot <- renderVisNetwork({
    req(mirnaNetworkData())
    set.seed(1456)
    visNetwork(mirnaNetworkData()$network_obj$nodes, mirnaNetworkData()$network_obj$edges) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = list(enabled = TRUE, values = mirnaNetworkData()$network_obj$nodes$id)) %>%
      visIgraphLayout(layout = "layout_with_fr", weights = mirnaNetworkData()$network_obj$edges$Score) %>%
      visLegend(position = "right", useGroups = FALSE, addNodes = mirnaNetworkData()$legend_obj, stepY = 75, width = .2)
  })

  
  mirnaNetworkTableData <- eventReactive(input$plotMiRnaNetwork,{
    req(mirnaNetworkData())
    if(!input$mirnaNetworkCircNormCounts){
      circ.norm.summary$Mean.Expression <- circ.norm.summary$Mean.SRPBM
    }
    circ.expression <- pivot_wider(circ.norm.summary, id_cols = c(Tissue, circRNA), values_from = Mean.Expression, names_from = Tissue)
    network_dat <- mirnaNetworkData()$network_obj$nodes
    dat <- scores_circ_mirna[scores_circ_mirna$circRNA %in% network_dat$id[network_dat$label == ""] &
                            scores_circ_mirna$miRNA %in% input$selectMirnaNetwork,] %>%
      left_join(circ.expression, by = "circRNA") %>%
      mutate_at(11:21, round, digits = 2) %>%
      rename("Parental Gene" = gene_symbol)
    dat$MRE.per.kb <- as.numeric(formatC(dat$MRE.per.kb, digits = 4))
    dat$Score <- as.numeric(formatC(dat$Score, digits = 4))
    dat$color <- network_dat$color[match(dat$circRNA, network_dat$id)]
    dat$Tissue <- network_dat$group[match(dat$circRNA, network_dat$id)]
    dat$circRNA <- paste0("<a href='#' onclick='detect_click(this)'>", dat$circRNA, "</a>")
    dat <- dat[,c("circRNA","ciri2_ID","circAtlas_ID","circBase_ID","miRNA","Parental Gene","Tissue","MRE","MRE.per.kb","Score", 
                  colnames(circ.expression)[-1],"color")]
    return(dat)
  })
  
  output$mirnaNetworkTable <- renderDT({
    withProgress(message = 'Printing table...', min = 0, max = 1, value = .4, {
      dat <- mirnaNetworkTableData()
      dat$miRNA <- NULL
      if(nrow(dat) > 100){
        lengthMenu <- c(10,25,50,100, nrow(dat))
      } else{
        lengthMenu <- c(10,25,50,100)
      }
      dat$Tissue <- paste0("<strong style='color:", dat$color, "'>",dat$Tissue,"</strong>")
      sketch = htmltools::withTags(table(
        class = 'display',
        thead(
          tr(
            th(colspan = 1, 'Position'),
            th(colspan = 1, 'ciri2'),
            th(colspan = 1, 'circAtlas (hg19)'),
            th(colspan = 1, 'circBase'),
            th(colspan = 1, 'Parental Gene'),
            th(colspan = 1, 'Tissue'),
            th(colspan = 1, 'MRE'),
            th(colspan = 1, 'MRE.per.kb'),
            th(colspan = 1, 'Pareto.Score'),
            th(colspan = 1, 'Adrenal gland', style="background-color: #CDAA7D; color: #ffffff"),
            th(colspan = 1, 'Gastro medialis', style="background-color: #0000FF; color: #ffffff"),
            th(colspan = 1, 'Heart', style="background-color: #EE3B3B; color: #ffffff"),
            th(colspan = 1, 'Liver', style="background-color: #66CD00; color: #ffffff"),
            th(colspan = 1, 'Lower leg skin', style="background-color: #EEAD0E; color: #ffffff"),
            th(colspan = 1, 'Ovary', style="background-color: #9932CC; color: #ffffff"),
            th(colspan = 1, 'Prostate', style="background-color: #CD6600; color: #ffffff"),
            th(colspan = 1, 'Testis', style="background-color: #FF69B4; color: #ffffff"),
            th(colspan = 1, 'Thyroid gland', style="background-color: #90EE90; color: #ffffff"),
            th(colspan = 1, 'Uterus', style="background-color: #7A67EE; color: #ffffff"),
            th(colspan = 1, 'Vagina', style="background-color: #737373; color: #ffffff")
          )
        )
      ))
      datatable(
        dat[,-ncol(dat)], 
        style = "bootstrap4",class = "table-condensed",
        rownames = FALSE, selection = "none", escape = FALSE,container = sketch,
        options = list(deferRender = TRUE,columnDefs = list(list(visible = FALSE, targets = c(1,2,3))), autowidth = TRUE, scrollX = TRUE, lengthMenu = lengthMenu,
                       rowCallback = JS("function(row, data) {",
                                        "var full_text = 'ciri2_ID: ' + data[1] + ', ' + 'circAtlas_ID: ' + data[2] + ', ' + 'circBase_ID: ' + data[3]",
                                        "$('td:eq(0)', row).attr('title', full_text)","}"))
      )
    })
  })
  
  output$downloadMirnaNetworkTable <- downloadHandler(
    filename = function() {
      paste(input$selectMirnaNetwork, "_network_table", '.csv', sep='')
    },
    content = function(con) {
      dat <- mirnaNetworkTableData()
      dat$circRNA <- cleanFun(dat$circRNA)
      dat$miRNA <- NULL
      dat <- dat[,-ncol(dat)]
      colnames(dat)[1:9] <- c("Position","ciri2","circAtlas (hg19)","circBase","Parental Gene","Tissue","MRE","MRE.per.kb","Pareto.Score")
      write.csv(dat, con, row.names = FALSE)
    }
  )
  
  output$circRnaExpressionPlot <- renderUI({
    bsModal("modalBarPlot",HTML(paste("<strong>",input$clicked_text,"</strong>")),"clicked_text", size = "large",
            checkboxInput("circRnaLog2Scale", "log2 scale", FALSE),
            plotOutput("selectedCircRnaExpression")  
    )
  })
  
  observeEvent(input$clicked_text,{
    toggleModal(session, "modalBarPlot", toggle = "toggle")
  })
  
  output$selectedCircRnaExpression <- renderPlot({
    dat <- circ.norm.summary
    if(!input$mirnaNetworkCircNormCounts){
      dat$Expression <- dat$Mean.SRPBM
      if(input$circRnaLog2Scale){
        dat$Expression <- dat$Mean.Log2.SRPBM
      }
    } else{
      dat$Expression <- dat$Mean.Expression
      if(input$circRnaLog2Scale){
        dat$Expression <- dat$Mean.Log2.Expression
      }
    }
    p <- ggplot(filter(dat, circRNA == input$clicked_text), aes(x = Tissue, y = Expression)) + 
      stat_summary(geom = "bar", aes(fill = Tissue), width = .4) + 
      scale_fill_manual(values = tissue_colors) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
            legend.position = "none")
    return(p)
  }, height = 400)
  
  heatmapData <- reactive({
    if(!input$mirnaNetworkCircNormCounts){
      circ_norm_counts <- circ_srpbm_counts
    }
    circRnaAnnots <- mirnaNetworkData()$network_obj$nodes[-1,] %>% 
      select(id, color, group)
    circRnaData <- data.frame(t(log2(circ_norm_counts+.1)[match(circRnaAnnots$id, rownames(circ_norm_counts), nomatch = 0),]), check.names = FALSE) %>%
      mutate(Tissue = design$Tissue) %>%
      gather("circRNA", "Expression",-Tissue) %>%
      group_by(Tissue, circRNA) %>%
      summarise(Expression = mean(Expression)) %>%
      spread("circRNA","Expression")
    circRnaData <- data.frame(circRnaData, check.names = FALSE)
    rownames(circRnaData) <- circRnaData$Tissue
    circRnaData$Tissue <- NULL
    circRnaData <- data.frame(t(circRnaData), check.names = FALSE)
    circRnaData <- circRnaData[match(circRnaAnnots$id, rownames(circRnaData), nomatch = 0),]
    row_side_colors <- data.frame(circRnaAnnots[,"group",drop = F])
    row_colors <- unlist(circRnaAnnots[,"color"][-which(duplicated(circRnaAnnots[,"color"]) == TRUE),])
    names(row_colors) <- unlist(circRnaAnnots[,"group"][-which(duplicated(circRnaAnnots[,"color"]) == TRUE),])
    row_dist <- dist(circRnaData)
    row_clust <- fastcluster::hclust(row_dist, method = "ward.D2")
    col_dist <- dist(t(circRnaData))
    col_clust <- fastcluster::hclust(col_dist, method = "ward.D2")
    circRnaData[circRnaData < 0] <- 0
    col_side_colors <- data.frame(Tissue = names(tissue_colors), stringsAsFactors = FALSE)
    col_side_colors <- col_side_colors[match(colnames(circRnaData), col_side_colors$Tissue, nomatch = 0),,drop = F]
    col_side_colors <- data.frame(rbind(col_side_colors, "unspecific"))
    col_side_colors$Tissue <- factor(col_side_colors$Tissue, levels = c(col_side_colors$Tissue[nrow(col_side_colors):1]))
    col_colors <- tissue_colors[match(colnames(circRnaData), names(tissue_colors), nomatch = 0)]
    colnames(circRnaData)[grep("Gastro", colnames(circRnaData))] <- "Gastro medialis"
    valData <- round(as.matrix(circRnaData), 4)
    genes <- annots[match(rownames(circRnaData), annots$circRNA, nomatch = 0),"gene_symbol"]
    rows <- rownames(circRnaData)
    cols <- colnames(circRnaData)
    colData <- t(valData)
    colData[] <- paste0("column: ", cols)
    colData <- t(colData)
    textData <- valData 
    textData[] <- paste0("row: ", rows, "\n","Gene Symbol: ", genes)
    textData <- matrix(paste0("value: ",valData, "\n", colData,"\n",textData), nrow = nrow(valData), ncol = ncol(valData))
    return(list(data = circRnaData, col_side_colors = col_side_colors, row_side_colors = row_side_colors, text_data = textData, 
                col_colors = col_colors, row_colors = row_colors, col_clust = col_clust, row_clust = row_clust))
  })
  
  output$plotCircRnaHeatmap <- renderPlotly({
    p <- heatmaply(heatmapData()$data, colors = colorRampPalette(c("navy","yellow","firebrick3"))(100), 
                   row_side_colors = heatmapData()$row_side_colors,
                   row_side_palette = heatmapData()$row_colors,
                   col_side_colors = heatmapData()$col_side_colors[-nrow(heatmapData()$col_side_colors),,drop=F],
                   col_side_palette = heatmapData()$col_colors,
                   showticklabels = c(TRUE,FALSE),
                   subplot_widths = c(.80,.03,.12), 
                   subplot_heights = c(.12, .03, .85), 
                   colorbar_xanchor = "right",
                   colorbar_xpos = 1.08,
                   colorbar_ypos = .56,
                   row_dend_left = TRUE,
                   custom_hovertext = heatmapData()$text_data,
                   plot_method = "plotly") 
    for(i in 1:length(p$x$data)){
      if(plotly:::has_attr(p$x$data[[i]]$type, "showscale")){
        if(length(grep("title", names(p$x$data[[i]]$colorbar))) > 0){
          if(p$x$data[[i]]$colorbar$title == "group"){
            p$x$data[[i]]$showscale <- plotly:::default(FALSE)
          } 
          if(p$x$data[[i]]$colorbar$title == "Tissue"){
            p$x$data[[i]]$colorbar$yanchor = "bottom"
            p$x$data[[i]]$colorbar$y = .2
          }
        }
      }
    }
    return(p)
  }) 
  
  circForFuncAnalysis <- reactive({
    req(input$circForAnalysis)
    ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
    if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
      circ.id <- annotationIdentifier(id = input$selectCircForAnalysis, id.map = ids, default.ids = colnames(circ_mirna_mat), id.group = "circ")
    } else{
      circ.id <- annotationIdentifier(id = input$circForAnalysis, id.map = ids, default.ids = colnames(circ_mirna_mat), id.group = "circ") 
    }
    return(circ.id)
  })
  
  output$circMultipleMapping <- renderUI({
    req(input$circForAnalysis)
    dat <- annots[annots$gene_symbol %in% input$circForAnalysis,]
    if(length(unique(dat$circRNA)) >= 1 & !is.null(input$circForAnalysis)){
      return(
        selectizeInput("selectCircForAnalysis", "circRNA by Position", choices = unique(dat$circRNA), multiple = TRUE,
                       selected = unique(dat$circRNA)[1],options = list(maxItems = 1, maxOptions = 500))
      )
    } else{
      return(NULL)
    }
  })
  
  circNetworkData <- reactive({
    req(circForFuncAnalysis())
    edges <- scores_circ_mirna[scores_circ_mirna$circRNA %in% circForFuncAnalysis(),] %>%
      select(circRNA, miRNA, MRE, MRE.per.kb, Score) %>%
      rename(from = circRNA, to = miRNA) %>%
      filter(Score > 0) %>%
      mutate(color = "grey")
    nodes <- tibble(id = unique(c(edges$from, edges$to))) %>%
      mutate(label = id) %>%
      left_join(tissue_specific_summary, by = c("id" = "circRNA")) %>%
      select(id, label, Tissue) %>%
      left_join(edges, by = c("id" = "to")) %>%
      select(id, label, Tissue, MRE, MRE.per.kb, Score) %>%
      mutate(type = c("circRNA", rep("miRNA", n()-1))) %>%
      mutate(color = tissue_colors[match(Tissue, names(tissue_colors))]) %>%
      mutate(color = ifelse(type == "miRNA", "turquoise", ifelse(is.na(Tissue), "#B4C3DE", color))) %>%
      mutate(size = ifelse(type == "circRNA", 40, 30)) %>%
      mutate(shape = ifelse(type == "circRNA", "triangle", "dot")) %>%
      mutate(font.size = ifelse(length(unique(id)) > 30, 30, length(unique(id)))) %>%
      rename(group = Tissue) %>%
      mutate(group = ifelse(color == "turquoise", "", ifelse(is.na(group), "unspecific", as.character(group)))) %>%
      mutate(title = ifelse(type == "circRNA", paste0("<p>",id,"<br>","<strong>circBase:</strong> ", annots$circBase_ID[match(id, annots$circRNA)], "<br>","<strong>circAtlas (hg19):</strong> ",
                                                      annots$circAtlas_ID[match(id, annots$circRNA)], "<br>", "<strong>ciri2:</strong> ", annots$ciri2_ID[match(id, annots$circRNA)],
                                                      "<br>","<strong>Parental Gene:</strong> ", annots$gene_symbol[match(id, annots$circRNA)],
                                                      "<br>", "<strong>Tissue:</strong> ", group),
                            paste0("<p>",label,"<br>","Pareto Score: ", round(Score, 4),"<br>", "MRE: ", MRE, "<br>", "MRE per KB: ", round(MRE.per.kb, 4))))
    nodes$label[nodes$type == "circRNA"] <- input$circForAnalysis
    labels <- c("circRNA","miRNA") 
    color.background <- c(nodes$color[nodes$type == "circRNA"], "turquoise")
    color.border <- c(nodes$color[nodes$type == "circRNA"], "turquoise")
    shape <- c("triangle","dot")
    legend_dat <- data.frame(label = labels,
                             color.background = color.background,
                             color.border = color.border)
    legend_dat$shape <- shape
    legend_dat$size <- 10
    vis_obj <- list(nodes = nodes, edges = edges, legend_dat = legend_dat)
    return(vis_obj)
  })
  
  output$circNetworkPlot <- renderVisNetwork({
    req(circNetworkData())
    if(!is.null(input$circMirnaNetworkTable_rows_selected)){
      edges <- circNetworkData()$edges[circNetworkData()$edges$to %in% mirnasSelected(),]
      nodes <- circNetworkData()$nodes[circNetworkData()$nodes$id %in% c(edges$from, edges$to),]
      nodes$font.size <- ifelse(nrow(nodes) > 30, 30, nrow(nodes))
    } else{
      edges <- circNetworkData()$edges
      nodes <- circNetworkData()$nodes
    }
    set.seed(1456)
    visNetwork(nodes, edges) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = list(enabled = TRUE, values = nodes$id)) %>%
      visIgraphLayout(layout = "layout_with_fr", weights = edges$Score) %>%
      visLegend(position = "right", useGroups = FALSE, addNodes = circNetworkData()$legend_dat, stepY = 75, width = .2)
  })
  
  circMirnaNetworkData <- reactive({
    req(circNetworkData())
    dat <- circNetworkData()$edges[,2:5] %>% 
      rename(miRNA = to)
    dat$MRE.per.kb <- as.numeric(formatC(dat$MRE.per.kb, digits = 4))
    dat$Score <- as.numeric(formatC(dat$Score, digits = 4))
    dat <- dat[order(dat$Score, decreasing = TRUE),]
    if(input$circNetworkTableCpmCounts){
      mirna.norm.summary$Mean.Expression <- mirna.norm.summary$Mean.CPM
    }
    mirna.expression <- pivot_wider(mirna.norm.summary, id_cols = c(Tissue, miRNA), values_from = Mean.Expression, names_from = Tissue)
    dat <- left_join(dat, mirna.expression, by = "miRNA") %>%
      mutate_at(5:15, round, digits = 2)
    dat$miRNA <- paste0("<a href='#' onmousedown='event.preventDefault(); event.stopPropagation()' onclick='detect_click_mirna(this)'>", dat$miRNA, "</a>")
    return(dat)
  })
  
  output$mirnaExpressionPlot <- renderUI({
    bsModal("modalMirnaPlot",HTML(paste("<strong>",input$clicked_mirna_text,"</strong>")),"clicked_mirna_text", size = "large",
            checkboxInput("mirnaLog2Scale1", "log2 scale", FALSE),
            plotOutput("selectedMirnaExpression1")  
    )
  })
  
  observeEvent(input$clicked_mirna_text,{
    toggleModal(session, "modalMirnaPlot", toggle = "toggle")
  })
  
  output$selectedMirnaExpression1 <- renderPlot({
    dat <- mirna.norm.summary
    if(input$circNetworkTableCpmCounts){
      dat$Expression <- dat$Mean.CPM
      ylabel <- "CPM" 
      if(input$mirnaLog2Scale1){
        dat$Expression <- dat$Mean.Log2.CPM
      }
    } else{
      dat$Expression <- dat$Mean.Expression
      ylabel <- "Normalized Expression"
      if(input$mirnaLog2Scale1){
        dat$Expression <- dat$Mean.Log2.Expression
      }
    }
    p <- ggplot(filter(dat, miRNA == input$clicked_mirna_text), aes(x = Tissue, y = Expression)) + 
      stat_summary(geom = "bar", aes(fill = Tissue), width = .4) + 
      scale_fill_manual(values = tissue_colors) + ylab(ylabel) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
            legend.position = "none")
    return(p)
  }, height = 400)
  
  output$circMirnaNetworkTable <- renderDT({
    sketch = htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(colspan = 1, 'miRNA'),
          th(colspan = 1, 'MRE'),
          th(colspan = 1, 'MRE.per.kb'),
          th(colspan = 1, 'Pareto.Score'),
          th(colspan = 1, 'Adrenal gland', style="background-color: #CDAA7D; color: #ffffff"),
          th(colspan = 1, 'Gastro medialis', style="background-color: #0000FF; color: #ffffff"),
          th(colspan = 1, 'Heart', style="background-color: #EE3B3B; color: #ffffff"),
          th(colspan = 1, 'Liver', style="background-color: #66CD00; color: #ffffff"),
          th(colspan = 1, 'Lower leg skin', style="background-color: #EEAD0E; color: #ffffff"),
          th(colspan = 1, 'Ovary', style="background-color: #9932CC; color: #ffffff"),
          th(colspan = 1, 'Prostate', style="background-color: #CD6600; color: #ffffff"),
          th(colspan = 1, 'Testis', style="background-color: #FF69B4; color: #ffffff"),
          th(colspan = 1, 'Thyroid gland', style="background-color: #90EE90; color: #ffffff"),
          th(colspan = 1, 'Uterus', style="background-color: #7A67EE; color: #ffffff"),
          th(colspan = 1, 'Vagina', style="background-color: #737373; color: #ffffff")
        )
      )
    ))
    datatable(
      circMirnaNetworkData(), style = "bootstrap4", container = sketch, escape = FALSE,
      selection = "multiple", class = "table-condensed", rownames = FALSE,
      options = list(autowidth = TRUE, scrollX = TRUE)
    )
  })
  
  output$downloadCircMirnaNetworkTable <- downloadHandler(
    filename = function() {
      ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
      if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
        paste(circForFuncAnalysis(), "_network_table", '.csv', sep='')
      } else{
        paste(input$circForAnalysis, "_network_table", '.csv', sep='')
      }
    },
    content = function(con) {
      dat <- circMirnaNetworkData()
      colnames(dat)[1:4] <- c("miRNA","MRE","MRE.per.kb","Pareto.Score")
      write.csv(dat, con, row.names = FALSE)
    }
  )
  
  mirnasSelected <- reactive({
    if(is.null(input$circMirnaNetworkTable_rows_selected)){return(NULL)}
    mirnas <- cleanFun(circMirnaNetworkData()$miRNA[input$circMirnaNetworkTable_rows_selected])
    return(mirnas)
  })
  
  output$circSelected <- renderUI({
    req(circForFuncAnalysis())
    ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
    if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
      HTML(paste("<h5><strong>Position</strong><br>",input$selectCircForAnalysis,"</h5>"))
    } else{
      circSelected <- annots$circRNA[annots$circRNA == circForFuncAnalysis()]
      HTML(paste("<h5><strong>Position</strong><br>",circSelected,"</h5>"))
    }
  })
  
  output$circBaseID <- renderUI({
    req(circForFuncAnalysis())
    ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
    if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
      circBaseID <- annots$circBase_ID[annots$circRNA == input$selectCircForAnalysis]
      HTML(paste("<h5><strong>circBase</strong><br>",circBaseID,"</h5>"))
    } else{
      circBaseID <- annots$circBase_ID[annots$circRNA == circForFuncAnalysis()]
      HTML(paste("<h5><strong>circBase</strong><br>",circBaseID,"</h5>"))
    }
  })
  
  output$circAtlasID <- renderUI({
    req(circForFuncAnalysis())
    ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
    if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
      circAtlasID <- annots$circAtlas_ID[annots$circRNA == input$selectCircForAnalysis]
      HTML(paste("<h5><strong>circAtlas (hg19)</strong><br>",circAtlasID,"</h5>"))
    } else{
      circAtlasID <- annots$circAtlas_ID[annots$circRNA == circForFuncAnalysis()]
      HTML(paste("<h5><strong>circAtlas (hg19)</strong><br>",circAtlasID,"</h5>"))
    }
  })
  
  output$ciri2ID <- renderUI({
    req(circForFuncAnalysis())
    ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
    if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
      ciri2ID <- annots$ciri2_ID[annots$circRNA == input$selectCircForAnalysis]
      HTML(paste("<h5><strong>ciri2</strong><br>",ciri2ID,"</h5>"))
    } else{
      ciri2ID <- annots$ciri2_ID[annots$circRNA == circForFuncAnalysis()]
      HTML(paste("<h5><strong>ciri2</strong><br>",ciri2ID,"</h5>"))
    }
  })
  
  output$parentalGeneSymbol <- renderUI({
    req(circForFuncAnalysis())
    ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
    if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
      parentalGeneSymbol <- annots$gene_symbol[annots$circRNA == input$selectCircForAnalysis]
      HTML(paste("<h5><strong>Parental Gene</strong><br>",parentalGeneSymbol,"</h5>"))
    } else{
      parentalGeneSymbol <- annots$gene_symbol[annots$circRNA == circForFuncAnalysis()]
      HTML(paste("<h5><strong>Parental Gene</strong><br>",parentalGeneSymbol,"</h5>"))
    }
  })
  
  output$numberMirnasSelected <- renderUI({
    req(circForFuncAnalysis())
    number.mirnas <- nrow(circNetworkData()$edges)
    if(!is.null(mirnasSelected())){
      number.mirnas <- length(mirnasSelected())
    }
    HTML(paste("<h5><strong>Number of miRNAs</strong><br>",number.mirnas,"</h5>")) 
  })
  
  circSigGenes <- eventReactive(input$runGeneAnalysis,{
    if(!is.null(mirnasSelected())){
      res <- circRnaFuncAnalysis(circ_mirna_mat, mirna_linrna_mat, circForFuncAnalysis(), mirnas = mirnasSelected())
    } else{
      res <- circRnaFuncAnalysis(circ_mirna_mat, mirna_linrna_mat, circForFuncAnalysis())
    }
    return(list(pvals = res$p, mirnas = res$mirnas, scores = res$scores, circ.id = circForFuncAnalysis()))
  })
  
  geneMirnaTable <- reactive({
    req(circSigGenes())
    res <- circSigGenes()
    gene_mirna_tab <- scores_mirna_linrna[scores_mirna_linrna$miRNA %in% res$mirnas, ] %>%
      rename(Gene = linRNA, mirTarBase = MirTarBase, "miR-gene Score" = Score, "miR-gene MRE" = MRE) %>%
      select(Gene, miRNA, "miR-gene MRE", mirTarBase, "miR-gene Score")
    gene_mirna_tab$pvalue <- res$pvals[match(gene_mirna_tab$Gene, names(res$pvals))]
    gene_mirna_tab$fdr <- p.adjust(res$pvals, "fdr")[match(gene_mirna_tab$Gene, names(res$pvals))]
    gene_mirna_tab$`miR-gene Score` <- as.numeric(formatC(gene_mirna_tab$`miR-gene Score`, digits = 4))
    gene_mirna_tab$pvalue <- as.numeric(formatC(gene_mirna_tab$pvalue, digits = 4))
    gene_mirna_tab$fdr <- as.numeric(formatC(gene_mirna_tab$fdr, digits = 4))
    circMirnaNetworkData <- circMirnaNetworkData()
    circMirnaNetworkData$miRNA <- cleanFun(circMirnaNetworkData$miRNA)
    gene_mirna_tab <- left_join(gene_mirna_tab, circMirnaNetworkData, by = "miRNA") %>%
      select(Gene, miRNA, pvalue, fdr, "miR-gene MRE", mirTarBase, "miR-gene Score", MRE, MRE.per.kb, Score)
    return(gene_mirna_tab)
  })
  
  numberSigGenes <- reactive({
    pvals <- data.frame(pvalue = circSigGenes()$pvals, fdr = p.adjust(circSigGenes()$pvals, "fdr"))
    summary.table <- data.frame(pvalue = seq(.01,1,.01))
    summary.table$Unadjusted <- sapply(summary.table$pvalue, function(x) sum(pvals$pvalue <= x))
    summary.table$FDR <- sapply(summary.table$pvalue, function(x) sum(pvals$fdr <= x))
    summary.plot <- pvals
    summary.plot$raw.less.than <- rank(summary.plot$pvalue, ties.method = "max", na.last = "keep")
    summary.plot$fdr.less.than <- rank(summary.plot$fdr, ties.method = "max", na.last = "keep")
    rownames(summary.plot) <- names(circSigGenes()$pvals)
    return(list(table = summary.table, plot = summary.plot))
  })
  
  output$numberSigGenePlot <- renderPlot({
    req(numberSigGenes())
    withProgress(message = 'Generating figure...', min = 0, max = 1, value = .4, {
      p <- ggplot(numberSigGenes()$plot, aes(x = pvalue)) + geom_line(aes(y = raw.less.than, colour = "Unadjusted"), size = 2) + 
        geom_line(aes(x = fdr, y = fdr.less.than, colour = "FDR"), size = 2) + 
        scale_colour_manual(name = "", values = c(Unadjusted = "#67A9CF", FDR = "#e66065")) + 
        theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14)) +
        ylab("Number significant") + xlab("P-Value") + theme_classic()
      print(p)
    })
  })
  
  output$numberSigGeneTable <- renderDT({
    req(numberSigGenes())
    withProgress(message = 'Generating table...', min = 0, max = 1, value = .4, {
      #numberSigGenes()$table
      datatable(
        numberSigGenes()$table,
        style = "bootstrap4",
        class = "table-condensed",
        rownames = FALSE,
        selection = "none",
        options = list(autowidth = TRUE, scrollX = TRUE) 
      ) 
    })
  })
  
  output$geneMirnaTable <- renderDT({
    #req(numberSigGenes())
    withProgress(message = 'Generating table...', min = 0, max = 1, value = .4, {
      gene_mirna_tab <- geneMirnaTable()
      if(input$geneLevelTestFdrAdjusted){
        gene_mirna_tab <- gene_mirna_tab[gene_mirna_tab$fdr <= input$geneLevelTestThresh, ] 
      } else{
        gene_mirna_tab <- gene_mirna_tab[gene_mirna_tab$pvalue <= input$geneLevelTestThresh, ]
      }
      sketch = htmltools::withTags(table(
        class = 'display',
        thead(
          tr(
            th(colspan = 2, ''),
            th(colspan = 2, 'Gene', style='text-align:center; background-color: #ADD8E6'),
            th(colspan = 3, 'miRNA-Gene', style='text-align:center; background-color: #fed8b1'),
            th(colspan = 3, 'circRNA-miRNA', style='text-align:center; background-color: #90ee90')
          ),
          tr(
            th(colspan = 1, 'Gene'),
            th(colspan = 1, 'miRNA'),
            th(colspan = 1, 'pvalue'),
            th(colspan = 1, 'FDR'),
            th(colspan = 1, 'MRE'),
            th(colspan = 1, 'mirTarBase'),
            th(colspan = 1, 'Pareto.Score'),
            th(colspan = 1, 'MRE'),
            th(colspan = 1, 'MRE.per.KB'),
            th(colspan = 1, 'Pareto.Score')
          )
        )
      ))
      searchable <- FALSE
      disable.searching <- c(2:(ncol(gene_mirna_tab)-1))
      datatable(
        gene_mirna_tab,
        style = "bootstrap4",
        class = "table-condensed",
        container = sketch,
        rownames = FALSE,
        selection = "none",
        filter = list(position = "top", clear = FALSE),
        options = list(autowidth = TRUE, scrollX = TRUE,deferRender = TRUE,
                       columnDefs = list(list(searchable = searchable, targets = disable.searching)),
                       searchHighlight = TRUE, search = list(regex = TRUE))
      ) 
    })
  })
  
  output$downloadNumberSigGeneTable <- downloadHandler(
    filename = function() {
      ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
      if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
        paste(circForFuncAnalysis(), "_number_significant_genes_summary", '.csv', sep='')
      } else{
        paste(input$circForAnalysis, "_number_significant_genes_summary", '.csv', sep='') 
      }
    },
    content = function(con) {
      write.csv(numberSigGenes()$table, con, row.names = FALSE)
    }
  )
  
  output$downloadSigGeneList <- downloadHandler(
    filename = function() {
      ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
      if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
        paste(circForFuncAnalysis(), "_gene_list", '.csv', sep='')
      } else{
        paste(input$circForAnalysis, "_gene_list", '.csv', sep='')
      }
    },
    content = function(con) {
      dat <- data.frame(pvalue = circSigGenes()$pvals, row.names = names(circSigGenes()$pvals))
      dat <- dat[order(dat$pvalue, decreasing = FALSE),,drop=F]
      write.csv(dat, con, row.names = TRUE)
    }
  )
  
  output$downloadGeneMirnaTable <- downloadHandler(
    filename = function() {
      gene_mirna_tab <- geneMirnaTable()
      if(input$geneLevelTestFdrAdjusted){
        gene_mirna_tab <- gene_mirna_tab[gene_mirna_tab$fdr <= input$geneLevelTestThresh, ] 
      } else{
        gene_mirna_tab <- gene_mirna_tab[gene_mirna_tab$pvalue <= input$geneLevelTestThresh, ]
      }
      ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
      if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
        paste(circForFuncAnalysis(), "_gene_mirna_table", '.csv', sep='')
      } else{
        paste(input$circForAnalysis, "_gene_mirna_table", '.csv', sep='') 
      }
    },
    content = function(con) {
      write.csv(gene_mirna_tab, con, row.names = FALSE)
    }
  )
  
  pathways <- reactive({
    if(input$termType == "KEGG"){
      pathways <- kegg_genes
      pathway_mat <- kegg_mat
    } else if(input$termType == "GO BP"){
      pathways <- go_bp_genes
      pathway_mat <- go_bp_mat
    } else if(input$termType == "GO CC"){
      pathways <- go_cc_genes
      pathway_mat <- go_cc_mat
    } else{
      pathways <- go_mf_genes
      pathway_mat <- go_mf_mat
    }
    return(list(pathways = pathways, pathway_mat = pathway_mat))
  })
  
  
  categorize <- reactiveVal(NULL)
  runFuncCounter <- reactiveVal(0)
  
  observeEvent(input$runFuncAnalysis,{
    counter <- runFuncCounter() + 1
    runFuncCounter(counter)
  })
  
  observeEvent(if(runFuncCounter() == 0){
    return(NULL)
  } else{
    input$categorizeBy
  }, {
    categorize(input$categorizeBy)
  })
  
  observeEvent(input$runGeneAnalysis, {
    req(geneMirnaTable())
    runFuncCounter(0)
    show(selector = '#enrichmentWorkflow li a[data-value="func_enrich_val"]')
  })
  
  circFuncAnalysis <- eventReactive(c(runFuncCounter()),{
    if(runFuncCounter() == 0){return(NULL)}
    withProgress(message = 'Running functional analysis...', min = 0, max = 1, value = .4, {
      if(input$fdrAdjustedGenes){
        genes <- names(circSigGenes()$pvals)[p.adjust(circSigGenes()$pvals, "fdr") <= input$genePvalThresh]
      } else{
        genes <- names(circSigGenes()$pvals)[circSigGenes()$pvals <= input$genePvalThresh]
      }
      res <- pathwayAnalysis(genes, pathways()$pathways, fdr.adjusted = FALSE, pval.thresh = 1)
      return(res)
    })
  })
  
  circFuncResTable <- reactive({
    req(circFuncAnalysis())
    mre <- scores_circ_mirna
    mre$circRNA <- mre$circRNA
    res <- circFuncAnalysis()
    res$pvalue <- as.numeric(formatC(res$pvalue, digits = 4))
    res$FDR <- as.numeric(formatC(res$FDR, digits = 4))
    res$Genes <- gsub("c(","",res$Genes, fixed = TRUE)
    res$Genes <- gsub(")","",res$Genes, fixed = TRUE)
    res$Genes <- gsub('"',"",res$Genes, fixed = TRUE)
    res <- separate_rows(res, "Genes")
    mirna.dat <- data.frame(mirna_linrna_mat[,circSigGenes()$mirnas], check.names = FALSE)
    mirna.dat$linRNA <- rownames(mirna.dat)
    mirna.linrna.scores <- pivot_longer(mirna.dat, -linRNA, names_to = "miRNA", values_to = "Score") %>%
      select(miRNA, linRNA, Score) %>% rename(Genes = linRNA)
    if(input$categorizeBy == "Function"){
      if(!input$fdrAdjustedTerms){
        res <- left_join(res, mirna.linrna.scores, by = "Genes") %>% filter((Score > 0 | Genes == "NA") & pvalue <= input$termPvalThresh) %>%
          select(Function, miRNA, Genes, Count, Function.Size, pvalue, FDR) %>% group_by(Function) %>% 
          mutate(miRNA.Count = length(unique(miRNA))) %>%
          summarise(miRNA = paste(unique(miRNA), collapse = ", "), Genes = paste(unique(Genes), collapse = ", "), miRNA.Count = unique(miRNA.Count),Count = unique(Count), Function.Size = unique(Function.Size),
                    pvalue = unique(pvalue), FDR = unique(FDR)) %>%
          arrange(pvalue) %>%
          rename("#miRNAs" = miRNA.Count, "#Genes" = Count)
      } else{
        res <- left_join(res, mirna.linrna.scores, by = "Genes") %>% filter((Score > 0 | Genes == "NA") & FDR <= input$termPvalThresh) %>%
          select(Function, miRNA, Genes, Count, Function.Size, pvalue, FDR) %>% group_by(Function) %>% 
          mutate(miRNA.Count = n()) %>%
          summarise(miRNA = paste(unique(miRNA), collapse = ", "), Genes = paste(unique(Genes), collapse = ", "), miRNA.Count = unique(miRNA.Count), Count = unique(Count), Function.Size = unique(Function.Size),
                    pvalue = unique(pvalue), FDR = unique(FDR)) %>%
          arrange(pvalue) %>%
          rename("#miRNAs" = miRNA.Count, "#Genes" = Count)
      }
    } else if(input$categorizeBy == "miRNA"){
      if(!input$fdrAdjustedTerms){
        res <- left_join(res, mirna.linrna.scores, by = "Genes") %>% filter(Score > 0 & pvalue <= input$termPvalThresh) %>%
          select(miRNA, Function, Genes) %>% group_by(miRNA) %>% 
          summarise(Function = paste(unique(Function), collapse = "; "), Genes = paste(unique(Genes), collapse = ", "))
      } else{
        res <- left_join(res, mirna.linrna.scores, by = "Genes") %>% filter(Score > 0 & FDR <= input$termPvalThresh) %>%
          select(miRNA, Function, Genes) %>% group_by(miRNA) %>% 
          summarise(Function = paste(unique(Function), collapse = "; "), Genes = paste(unique(Genes), collapse = ", "))
      }
      circ.norm.summary$circRNA <- circ.norm.summary$circRNA
      circ.tissues <- circ.norm.summary[circ.norm.summary$circRNA %in% circSigGenes()$circ.id & circ.norm.summary$Above.Median == 1,]
      if(nrow(circ.tissues) > 0){
        mirna.tissues <- mirna.norm.summary[mirna.norm.summary$miRNA %in% res$miRNA & mirna.norm.summary$Above.Median == 1 & 
                                              mirna.norm.summary$Tissue %in% circ.tissues$Tissue,]
        if(nrow(mirna.tissues) > 0){
          mirna.tissues <- left_join(mirna.tissues, res, by = "miRNA") %>% 
            select(miRNA, Function, Genes, Tissue) %>% group_by(miRNA) 
          mirna.tissues$color <- tissue_colors[match(mirna.tissues$Tissue, names(tissue_colors))]
          mirna.tissues$Tissue <- paste0("<strong style='color:", mirna.tissues$color, "'>",mirna.tissues$Tissue,"</strong>")
          mirna.tissues <- mirna.tissues %>% 
            summarise(Function = unique(Function), Genes = unique(Genes), Tissue = paste(unique(Tissue), collapse = ", ")) %>% 
            select(miRNA, Tissue)
          res <- left_join(res, mirna.tissues, by = "miRNA")
        }
      }
      mre <- mre[mre$circRNA %in% circSigGenes()$circ.id, ]
      mre <- mre[match(res$miRNA, mre$miRNA, nomatch = 0),]
      res$MRE <- mre$MRE
      res$MRE.per.kb <- as.numeric(formatC(mre$MRE.per.kb, digits = 4))
      res$Score <- as.numeric(formatC(mre$Score, digits = 4))
      res <- res[order(res$Score, decreasing = TRUE),]
    } else{
      if(!input$fdrAdjustedTerms){
        res <- left_join(res, mirna.linrna.scores, by = "Genes") %>% filter(Score > 0 & pvalue <= input$termPvalThresh) %>%
          select(Genes, miRNA, Function) %>% group_by(Genes) %>% 
          summarise(Function = paste(unique(Function), collapse = "; "), miRNA = paste(unique(miRNA), collapse = ", "))
      } else{
        res <- left_join(res, mirna.linrna.scores, by = "Genes") %>% filter(Score > 0 & FDR <= input$termPvalThresh) %>%
          select(Genes, miRNA, Function) %>% group_by(Genes) %>% 
          summarise(Function = paste(unique(Function), collapse = "; "), miRNA = paste(unique(miRNA), collapse = ", "))
      }
      res$pvalue <- as.numeric(formatC(geneMirnaTable()$pvalue[match(res$Genes, geneMirnaTable()$Gene, nomatch = 0)], digits = 4))
      res$FDR <- as.numeric(formatC(geneMirnaTable()$fdr[match(res$Genes, geneMirnaTable()$Gene, nomatch = 0)], digits = 4))
    }
    return(res)
  })
  
  circFunResTableLong <- reactive({
    req(circFuncAnalysis())
    mre <- scores_circ_mirna
    mre$circRNA <- mre$circRNA
    res <- circFuncAnalysis()
    res$pvalue <- as.numeric(formatC(res$pvalue, digits = 4))
    res$FDR <- as.numeric(formatC(res$FDR, digits = 4))
    res$Genes <- gsub("c(","",res$Genes, fixed = TRUE)
    res$Genes <- gsub(")","",res$Genes, fixed = TRUE)
    res$Genes <- gsub('"',"",res$Genes, fixed = TRUE)
    res <- separate_rows(res, "Genes")
    mirna.dat <- data.frame(mirna_linrna_mat[,circSigGenes()$mirnas], check.names = FALSE)
    mirna.dat$linRNA <- rownames(mirna.dat)
    mirna.linrna.scores <- pivot_longer(mirna.dat, -linRNA, names_to = "miRNA", values_to = "Score") %>%
      select(miRNA, linRNA, Score) %>% rename(Genes = linRNA)
    if(!input$fdrAdjustedTerms){
      res <- left_join(res, mirna.linrna.scores, by = "Genes") %>% filter((Score > 0 | Genes == "NA") & pvalue <= input$termPvalThresh) %>%
        select(Function, miRNA, Genes, Count, Function.Size, pvalue, FDR) %>% group_by(Function) %>% 
        mutate(miRNA.Count = length(unique(miRNA))) %>%
        summarise(miRNA = paste(unique(miRNA), collapse = ", "), Genes = paste(unique(Genes), collapse = ", "), miRNA.Count = unique(miRNA.Count),Count = unique(Count), Function.Size = unique(Function.Size),
                  pvalue = unique(pvalue), FDR = unique(FDR)) %>%
        arrange(pvalue) %>%
        rename("#miRNAs" = miRNA.Count, "#Genes" = Count) %>%
        separate_rows(miRNA, sep = ", ") %>% 
        separate_rows(Genes, sep = ", ") %>%
        left_join(geneMirnaTable(), by = c("miRNA", "Genes" = "Gene")) %>%
        rename(Function.pvalue = pvalue.x, Function.FDR = FDR, Gene.pvalue = pvalue.y, Gene.FDR = fdr) %>%
        filter(!is.na(Gene.pvalue))
      if(input$categorizeBy == "miRNA"){
        res <- res %>%
          select(miRNA, everything())
      } else if(input$categorizeBy == "Gene"){
        res <- res %>%
          select(Genes, everything())
      }
    } else{
      res <- left_join(res, mirna.linrna.scores, by = "Genes") %>% filter((Score > 0 | Genes == "NA") & FDR <= input$termPvalThresh) %>%
        select(Function, miRNA, Genes, Count, Function.Size, pvalue, FDR) %>% group_by(Function) %>% 
        mutate(miRNA.Count = n()) %>%
        summarise(miRNA = paste(unique(miRNA), collapse = ", "), Genes = paste(unique(Genes), collapse = ", "), miRNA.Count = unique(miRNA.Count), Count = unique(Count), Function.Size = unique(Function.Size),
                  pvalue = unique(pvalue), FDR = unique(FDR)) %>%
        arrange(pvalue) %>%
        rename("#miRNAs" = miRNA.Count, "#Genes" = Count) %>%
        separate_rows(miRNA, sep = ", ") %>% 
        separate_rows(Genes, sep = ", ") %>%
        left_join(geneMirnaTable(), by = c("miRNA", "Genes" = "Gene")) %>%
        rename(Function.pvalue = pvalue.x, Function.FDR = FDR, Gene.pvalue = pvalue.y, Gene.FDR = fdr) %>%
        filter(!is.na(Gene.pvalue))
      if(input$categorizeBy == "miRNA"){
        res <- res %>%
          select(miRNA, everything())
      } else if(input$categorizeBy == "Gene"){
        res <- res %>%
          select(Genes, everything())
      }
    }
    return(res)
  })
 
  output$functionalResults <- renderDT({
    req(circFuncResTable())
    sketch = NULL
    dat <- circFuncResTable()
    if(input$longForm){
      dat <- circFunResTableLong()
      sketch = htmltools::withTags(table(
        class = 'display',
        thead(
          tr(
            th(colspan = 3, ''),
            th(colspan = 5, 'Function', style='text-align:center; background-color: #AB82FF'),
            th(colspan = 2, 'Gene', style='text-align:center; background-color: #ADD8E6'),
            th(colspan = 3, 'miRNA-Gene', style='text-align:center; background-color: #fed8b1'),
            th(colspan = 3, 'circRNA-miRNA', style='text-align:center; background-color: #90ee90')
          ),
          tr(
            th(colspan = 1, colnames(dat)[1]),
            th(colspan = 1, colnames(dat)[2]),
            th(colspan = 1, colnames(dat)[3]),
            th(colspan = 1, '#miRNAs'),
            th(colspan = 1, '#Genes'),
            th(colspan = 1, 'Size'),
            th(colspan = 1, 'pvalue'),
            th(colspan = 1, 'FDR'),
            th(colspan = 1, 'pvalue'),
            th(colspan = 1, 'FDR'),
            th(colspan = 1, 'MRE'),
            th(colspan = 1, 'mirTarBase'),
            th(colspan = 1, 'Pareto.Score'),
            th(colspan = 1, 'MRE'),
            th(colspan = 1, 'MRE.per.KB'),
            th(colspan = 1, 'Pareto.Score')
          )
        )
      ))
    }
    if(nrow(dat) > 100){
      lengthMenu <- c(10,25,50,100, nrow(dat))
    } else{
      lengthMenu <- c(10,25,50,100)
    }
    if(colnames(dat)[1] == "Function"){
      dat$Function <- paste0("<a href='#' onclick='detect_click_function(this)'>", dat$Function, "</a>")
    } else if(colnames(dat)[1] == "miRNA"){
      dat$miRNA <- paste0("<a href='#' onclick='detect_click_function(this)'>", dat$miRNA, "</a>")
      colnames(dat)[colnames(dat) == "Score"] <- "Pareto.Score"
    } else{
      dat$Genes <- paste0("<a href='#' onclick='detect_click_function(this)'>", dat$Genes, "</a>")
    }
    disable.searching <- c(4:ncol(dat))-1 
    searchable <- FALSE
    if(is.null(sketch)){
      datatable(
        dat, style = "bootstrap4", class = "table-condensed", rownames = FALSE, selection = "none", escape = FALSE,
        filter = list(position = "top", clear = FALSE),
        options = list(autowidth = TRUE, deferRender = TRUE, scrollX = TRUE, lengthMenu = lengthMenu,
                       columnDefs = list(list(searchable = searchable, targets = disable.searching)),
                       searchHighlight = TRUE, search = list(regex = TRUE))
      ) 
    } else{
      datatable(
        dat, style = "bootstrap4", class = "table-condensed", container = sketch, rownames = FALSE, selection = "none", escape = FALSE,
        filter = list(position = "top", clear = FALSE),
        options = list(autowidth = TRUE, deferRender = TRUE, scrollX = TRUE, lengthMenu = lengthMenu,
                       columnDefs = list(list(searchable = searchable, targets = disable.searching)),
                       searchHighlight = TRUE, search = list(regex = TRUE))
      )
    }
  })
  
  output$downloadFunctionalResults <- downloadHandler(
    filename = function() {
      ids <- annots[annots$circRNA %in% colnames(circ_mirna_mat),]
      if(isTruthy(input$selectCircForAnalysis) & length(unique(ids[ids$gene_symbol %in% input$circForAnalysis,]$circRNA)) > 1){
        paste(circForFuncAnalysis(), "_", input$termType, "_analysis", '.csv', sep='')
      } else{
        paste(input$circForAnalysis, "_", input$termType, "_analysis", '.csv', sep='')
      }
    },
    content = function(con) {
      dat <- circFuncResTable()
      if(input$longForm){
        dat <- circFunResTableLong()
      }
      write.csv(dat, con, row.names = FALSE)
    }
  )
  
  output$selectFunctionOrMirna <- renderUI({
    req(circFunResTableLong())
    if(input$categorizeBy == "miRNA"){
      terms <- unique(circFunResTableLong()$Function[circFunResTableLong()$miRNA == input$clicked_text_function])
      genes <- unique(circFunResTableLong()$Genes[circFunResTableLong()$miRNA == input$clicked_text_function])
      return(
        tagList(
          div(
            div(class = "row",
                div(class = "col-xs-4",
                    selectizeInput("selectMirnaForPlotting", "Select miRNA:", choices = unique(circFunResTableLong()$miRNA), 
                                   selected = input$clicked_text_function)
                ),
                div(class = "col-xs-4",
                    selectizeInput("filterOnFunction", "Select Functions:", choices = terms, selected = NULL, multiple = TRUE)
                ),
                div(class = "col-xs-4",
                    selectizeInput("filterOnGenes", "Select Genes:", choices = genes, selected = NULL, multiple = TRUE)
                )
            ),
            withSpinner(visNetworkOutput("mirnaGeneNetworkPlot", height = "50em"))  
          ),
          br(),
          fluidRow(
            column(width = 6,
                   htmlOutput("circRnaExpressionTitle"),
                   div(class = "row",
                       div(class = "col-xs-2",
                           checkboxInput("funcCircSrpbmCounts", "SRPBM", TRUE)
                       ),
                       div(class = "col-xs-6",
                           checkboxInput("funcCircNormCounts","DESeq2 normalized",FALSE)
                       )
                   ),
                   checkboxInput("funcCircLog2Scale", "log2 scale", FALSE),
                   plotOutput("funcCircRnaExpression")
            ),
            column(width = 6,
                   htmlOutput("mirnaExpressionTitle"),
                   div(class = "row",
                       div(class = "col-xs-2",
                           checkboxInput("funcMirnaCpmCounts", "CPM", TRUE)
                       ),
                       div(class = "col-xs-6",
                           checkboxInput("funcMirnaNormCounts","DESeq2 normalized",FALSE)
                       )
                   ),
                   checkboxInput("funcMirnaLog2Scale", "log2 scale", FALSE),
                   plotOutput("funcMiRnaExpression")
            ) 
          )
        )
      )
    } else if(input$categorizeBy == "Function"){
      genes <- unique(circFunResTableLong()$Genes[circFunResTableLong()$Function == input$clicked_text_function])
      mirnas <- unique(circFunResTableLong()$miRNA[circFunResTableLong()$Function == input$clicked_text_function])
      return(
        tagList(
          div(class = "row",
              div(class = "col-xs-4",
                  selectizeInput("selectFunctionForPlotting", "Select Function:", choices = unique(circFunResTableLong()$Function), 
                                 selected = input$clicked_text_function)
              ),
              div(class = "col-xs-4",
                  selectizeInput("filterOnMirna", "Select miRNAs:", choices = mirnas, selected = NULL, multiple = TRUE)
              ),
              div(class = "col-xs-4",
                  selectizeInput("filterOnGene", "Select Genes:", choices = genes, selected = NULL, multiple = TRUE)
              )
          ),
          withSpinner(visNetworkOutput("mirnaGeneNetworkPlot", height = "50em"))    
        )
      )
    } else{
      terms <- unique(circFunResTableLong()$Function[circFunResTableLong()$Genes == input$clicked_text_function])
      mirnas <- unique(circFunResTableLong()$miRNA[circFunResTableLong()$Genes == input$clicked_text_function])
      return(
        tagList(
          div(class = "row",
              div(class = "col-xs-4",
                  selectizeInput("selectGeneForPlotting", "Select Gene:", choices = unique(circFunResTableLong()$Genes), 
                                 selected = input$clicked_text_function)
              ),
              div(class = "col-xs-4",
                  selectizeInput("filterOnMirnas", "Select miRNAs:", choices = mirnas, selected = NULL, multiple = TRUE)
              ),
              div(class = "col-xs-4",
                  selectizeInput("filterOnFunctions", "Select Functions:", choices = terms, selected = NULL, multiple = TRUE)
              )
          ),
          withSpinner(visNetworkOutput("mirnaGeneNetworkPlot", height = "50em")),
          br(),
          fluidRow(
            column(width = 6,
                   htmlOutput("circRnaExpressionTitle"),
                   div(class = "row",
                       div(class = "col-xs-2",
                           checkboxInput("funcCircSrpbmCounts", "SRPBM", TRUE)
                       ),
                       div(class = "col-xs-6",
                           checkboxInput("funcCircNormCounts","DESeq2 normalized",FALSE)
                       )
                   ),
                   checkboxInput("funcCircLog2Scale", "log2 scale", FALSE),
                   plotOutput("funcCircRnaExpression")
            ),
            column(width = 6,
                   htmlOutput("geneExpressionTitle"),
                   div(class = "row",
                       div(class = "col-xs-2",
                           checkboxInput("funcGeneCpmCounts", "CPM", TRUE)
                       ),
                       div(class = "col-xs-6",
                           checkboxInput("funcGeneNormCounts","DESeq2 normalized",FALSE)
                       )
                   ),
                   checkboxInput("funcGeneLog2Scale", "log2 scale", FALSE),
                   plotOutput("funcGeneExpression")
            ) 
          )
        )
      )
    }
  })
  
  observe({
    req(circFunResTableLong())
    if(input$categorizeBy == "Function"){
      genes <- unique(circFunResTableLong()$Genes[circFunResTableLong()$Function == input$selectFunctionForPlotting])
      mirnas <- unique(circFunResTableLong()$miRNA[circFunResTableLong()$Function == input$selectFunctionForPlotting])
      updateSelectizeInput(session, "filterOnMirna", "Select miRNAs:", choices = mirnas, selected = NULL)
      updateSelectizeInput(session, "filterOnGene", "Select Genes:", choices = genes, selected = NULL)
    } else if(input$categorizeBy == "miRNA"){
      terms <- unique(circFunResTableLong()$Function[circFunResTableLong()$miRNA == input$selectMirnaForPlotting])
      genes <- unique(circFunResTableLong()$Genes[circFunResTableLong()$miRNA == input$selectMirnaForPlotting])
      updateSelectizeInput(session, "filterOnFunction", "Select Functions:", choices = terms, selected = NULL)
      updateSelectizeInput(session, "filterOnGenes", "Select Genes:", choices = genes, selected = NULL) 
    } else{
      terms <- unique(circFunResTableLong()$Function[circFunResTableLong()$Genes == input$selectGeneForPlotting])
      mirnas <- unique(circFunResTableLong()$miRNA[circFunResTableLong()$Genes == input$selectGeneForPlotting])
      updateSelectizeInput(session, "filterOnFunctions", "Select Functions:", choices = terms, selected = NULL)
      updateSelectizeInput(session, "filterOnMirnas", "Select miRNAs:", choices = mirnas, selected = NULL)
    }
  })
  
  
  observeEvent(input$filterOnMirna,{
    updateSelectizeInput(session, "filterOnGene", "Select Genes:", choices = filterNodes()$genes_for_select, selected = input$filterOnGene) 
  })
  
  observeEvent(input$filterOnGene,{
    updateSelectizeInput(session, "filterOnMirna", "Select miRNAs:", choices = filterNodes()$mirnas_for_select, selected = input$filterOnMirna) 
  })
  
  observeEvent({if(is.null(input$filterOnGene) & is.null(input$filterOnMirna)){
    return(TRUE)
  } else{
    return(NULL)
  }},{
    updateSelectizeInput(session, "filterOnGene", "Select Genes:", choices = filterNodes()$genes_for_select, selected = input$filterOnGene)
    updateSelectizeInput(session, "filterOnMirna", "Select miRNAs:", choices = filterNodes()$mirnas_for_select, selected = input$filterOnMirna)
  })
  
  observeEvent(input$filterOnFunction,{
    updateSelectizeInput(session, "filterOnGenes", "Select Genes:", choices = filterNodes()$genes_for_select, selected = input$filterOnGenes) 
  })
  
  observeEvent(input$filterOnGenes,{
    updateSelectizeInput(session, "filterOnFunction", "Select Functions:", choices = filterNodes()$terms_for_select, selected = input$filterOnFunction) 
  })
  
  observeEvent({if(is.null(input$filterOnGenes) & is.null(input$filterOnFunction)){
    return(TRUE)
  } else{
    return(NULL)
  }},{
    updateSelectizeInput(session, "filterOnGenes", "Select Genes:", choices = filterNodes()$genes_for_select, selected = input$filterOnGenes)
    updateSelectizeInput(session, "filterOnFunction", "Select Functions:", choices = filterNodes()$terms_for_select, selected = input$filterOnFunction)
  })
  
  observeEvent(input$filterOnFunctions,{
    updateSelectizeInput(session, "filterOnMirnas", "Select miRNAs:", choices = filterNodes()$mirnas_for_select, selected = input$filterOnMirnas) 
  })
  
  observeEvent(input$filterOnMirnas,{
    updateSelectizeInput(session, "filterOnFunctions", "Select Functions:", choices = filterNodes()$terms_for_select, selected = input$filterOnFunctions) 
  })
  
  observeEvent({if(is.null(input$filterOnFunctions) & is.null(input$filterOnMirnas)){
    return(TRUE)
  } else{
    return(NULL)
  }},{
    updateSelectizeInput(session, "filterOnFunctions", "Select Functions:", choices = filterNodes()$terms_for_select, selected = input$filterOnFunctions)
    updateSelectizeInput(session, "filterOnMirnas", "Select miRNAs:", choices = filterNodes()$mirnas_for_select, selected = input$filterOnMirnas)
  })
  
  filterNodes <- reactive({
    req(circFunResTableLong())
    if(input$categorizeBy == "Function"){
      genes <- unique(circFunResTableLong()$Genes[circFunResTableLong()$Function == input$selectFunctionForPlotting])
      mirnas <- unique(circFunResTableLong()$miRNA[circFunResTableLong()$Function == input$selectFunctionForPlotting])
      if(!is.null(input$filterOnMirna) | !is.null(input$filterOnGene)){
        if(!is.null(input$filterOnMirna) & !is.null(input$filterOnGene)){
          mirnas_to_remove <- c(setdiff(mirnas, input$filterOnMirna),
                                setdiff(mirnas, scores_mirna_linrna$miRNA[scores_mirna_linrna$linRNA %in% input$filterOnGene]))
          genes_to_remove <- c(setdiff(genes, input$filterOnGene),
                               setdiff(genes, scores_mirna_linrna$linRNA[scores_mirna_linrna$miRNA %in% input$filterOnMirna]))
          mirnas_for_select <- intersect(mirnas,scores_mirna_linrna$miRNA[scores_mirna_linrna$linRNA %in% input$filterOnGene])
          genes_for_select <- intersect(genes,scores_mirna_linrna$linRNA[scores_mirna_linrna$miRNA %in% input$filterOnMirna])
          if(length(genes_for_select) == 0){
            genes_for_select <- genes
          } 
        }
        if(!is.null(input$filterOnMirna) & is.null(input$filterOnGene)){
          mirnas_to_remove <- setdiff(mirnas, input$filterOnMirna)
          genes_to_remove <- setdiff(genes, scores_mirna_linrna$linRNA[scores_mirna_linrna$miRNA %in% input$filterOnMirna])
          mirnas_for_select <- mirnas
          genes_for_select <- intersect(genes,scores_mirna_linrna$linRNA[scores_mirna_linrna$miRNA %in% input$filterOnMirna])
          if(length(genes_to_remove) == 0){
            genes_to_remove <- NULL
          } 
        }
        if(!is.null(input$filterOnGene) & is.null(input$filterOnMirna)){
          genes_to_remove <- setdiff(genes, input$filterOnGene)
          mirnas_to_remove <- setdiff(mirnas, scores_mirna_linrna$miRNA[scores_mirna_linrna$linRNA %in% input$filterOnGene])
          genes_for_select <- genes
          mirnas_for_select <- intersect(mirnas,scores_mirna_linrna$miRNA[scores_mirna_linrna$linRNA %in% input$filterOnGene])
          if(length(mirnas_to_remove) == 0){
            mirnas_to_remove <- NULL
          }
        }
      } else{
        mirnas_to_remove <- genes_to_remove <- NULL
        mirnas_for_select <- mirnas
        genes_for_select <- genes
      }
      return(list(mirnas = mirnas_to_remove, genes = genes_to_remove, mirnas_for_select = mirnas_for_select, genes_for_select = genes_for_select)) 
    } else if(input$categorizeBy == "miRNA"){
      terms <- unique(circFunResTableLong()$Function[circFunResTableLong()$miRNA == input$selectMirnaForPlotting])
      genes <- unique(circFunResTableLong()$Genes[circFunResTableLong()$miRNA == input$selectMirnaForPlotting])
      mat <- reshape2::melt(pathways()$pathway_mat) %>%
        rename(Gene = Var1, Function = Var2) %>%
        filter(value > 0)
      if(!is.null(input$filterOnFunction) | !is.null(input$filterOnGenes)){
        if(!is.null(input$filterOnFunction) & !is.null(input$filterOnGenes)){
          terms_to_remove <- c(setdiff(terms, input$filterOnFunction),
                                setdiff(terms, mat$Function[mat$Gene %in% input$filterOnGenes]))
          genes_to_remove <- c(setdiff(genes, input$filterOnGenes),
                               setdiff(genes, mat$Gene[mat$Function %in% input$filterOnFunction]))
          terms_for_select <- intersect(terms, mat$Function[mat$Gene %in% input$filterOnGenes])
          genes_for_select <- intersect(genes, mat$Gene[mat$Function %in% input$filterOnFunction])
          if(length(genes_for_select) == 0){
            genes_for_select <- genes
          } 
        }
        if(!is.null(input$filterOnFunction) & is.null(input$filterOnGenes)){
          terms_to_remove <- setdiff(terms, input$filterOnFunction)
          genes_to_remove <- setdiff(genes, mat$Gene[mat$Function %in% input$filterOnFunction])
          terms_for_select <- terms
          genes_for_select <- intersect(genes,mat$Gene[mat$Function %in% input$filterOnFunction])
          if(length(genes_to_remove) == 0){
            genes_to_remove <- NULL
          } 
        }
        if(!is.null(input$filterOnGenes) & is.null(input$filterOnFunction)){
          genes_to_remove <- setdiff(genes, input$filterOnGenes)
          terms_to_remove <- setdiff(terms, mat$Function[mat$Gene %in% input$filterOnGenes])
          genes_for_select <- genes
          terms_for_select <- intersect(terms, mat$Function[mat$Gene %in% input$filterOnGenes])
          if(length(terms_to_remove) == 0){
            terms_to_remove <- NULL
          }
        }
      } else{
        terms_to_remove <- genes_to_remove <- NULL
        terms_for_select <- terms
        genes_for_select <- genes
      }
      return(list(terms = terms_to_remove, genes = genes_to_remove, terms_for_select = terms_for_select, genes_for_select = genes_for_select))
    } else{
      terms <- unique(circFunResTableLong()$Function[circFunResTableLong()$Genes == input$selectGeneForPlotting])
      mirnas <- unique(circFunResTableLong()$miRNA[circFunResTableLong()$Genes == input$selectGeneForPlotting])
      mat <- circFunResTableLong()[circFunResTableLong()$Genes %in% input$selectGeneForPlotting,c("miRNA","Function")]
      if(!is.null(input$filterOnFunctions) | !is.null(input$filterOnMirnas)){
        if(!is.null(input$filterOnFunctions) & !is.null(input$filterOnMirnas)){
          terms_to_remove <- c(setdiff(terms, input$filterOnFunctions),
                               setdiff(terms, mat$Function[mat$miRNA %in% input$filterOnMirnas]))
          mirnas_to_remove <- c(setdiff(mirnas, input$filterOnMirnas),
                               setdiff(mirnas, mat$miRNA[mat$Function %in% input$filterOnFunctions]))
          terms_for_select <- intersect(terms, mat$Function[mat$miRNA %in% input$filterOnMirnas])
          mirnas_for_select <- intersect(mirnas, mat$miRNA[mat$Function %in% input$filterOnFunctions])
          if(length(mirnas_for_select) == 0){
            mirnas_for_select <- mirnas
          } 
        }
        if(!is.null(input$filterOnFunctions) & is.null(input$filterOnMirnas)){
          terms_to_remove <- setdiff(terms, input$filterOnFunctions)
          mirnas_to_remove <- setdiff(mirnas, mat$miRNA[mat$Function %in% input$filterOnFunctions])
          terms_for_select <- terms
          mirnas_for_select <- intersect(mirnas,mat$miRNA[mat$Function %in% input$filterOnFunctions])
          if(length(mirnas_to_remove) == 0){
            mirnas_to_remove <- NULL
          } 
          if(length(terms_to_remove) == 0){
            terms_to_remove <- NULL
          }
        }
        if(!is.null(input$filterOnMirnas) & is.null(input$filterOnFunctions)){
          mirnas_to_remove <- setdiff(mirnas, input$filterOnMirnas)
          terms_to_remove <- setdiff(terms, mat$Function[mat$miRNA %in% input$filterOnMirnas])
          mirnas_for_select <- mirnas
          terms_to_remove <- NULL
          terms_for_select <- intersect(terms, mat$Function[mat$miRNA %in% input$filterOnMirnas])
          if(length(terms_to_remove) == 0){
            terms_to_remove <- NULL
          }
          if(length(mirnas_to_remove) == 0){
            mirnas_to_remove <- NULL
          }
        }
      } else{
        terms_to_remove <- mirnas_to_remove <- NULL
        terms_for_select <- terms
        mirnas_for_select <- mirnas
      }
      return(list(terms = terms_to_remove, mirnas = mirnas_to_remove, terms_for_select = terms_for_select, mirnas_for_select = mirnas_for_select))
    }
  })
  
  observe({
    req(circForFuncAnalysis())
    req(circFunResTableLong())
    if(input$categorizeBy == "miRNA"){
      output$mirnaExpressionTitle <- renderUI({
        HTML(paste("<h3>",input$selectMirnaForPlotting,"</h3>"))
      })
      output$circRnaExpressionTitle <- renderUI({
        HTML(paste("<h3>",circForFuncAnalysis(),"</h3>"))
      })
      output$funcCircRnaExpression <- renderPlot({
        withProgress(message = 'Plotting figures...', min = 0, max = 1, value = .4, {
          circ.id <- circSigGenes()$circ.id
          dat <- circ.norm.summary
          if(input$funcCircSrpbmCounts){
            dat$Expression <- dat$Mean.SRPBM
            if(input$funcCircLog2Scale){
              dat$Expression <- dat$Mean.Log2.SRPBM
            }
          } else{
            dat$Expression <- dat$Mean.Expression
            if(input$funcCircLog2Scale){
              dat$Expression <- dat$Mean.Log2.Expression
            }
          }
          label.colors <- tissue_colors
          label.colors[!names(label.colors) %in% tissues()] <- "black"
          label.face <- rep("plain", length(tissue_colors))
          names(label.face) <- names(tissue_colors)
          label.face[names(label.face) %in% tissues()] <- "bold"
          p <- ggplot(filter(dat, circRNA == circ.id), aes(x = Tissue, y = Expression)) + 
            stat_summary(geom = "bar", aes(fill = Tissue), width = .4) + 
            scale_fill_manual(values = tissue_colors) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, color = label.colors, face = label.face), axis.text = element_text(size = 12), 
                  axis.title = element_text(size = 14),
                  legend.position = "none")
          return(p)
        })
      })
      output$funcMiRnaExpression <- renderPlot({
        if(input$categorizeBy == "Function"){return(NULL)}
        mirna.id <- input$selectMirnaForPlotting
        dat <- mirna.norm.summary
        if(input$funcMirnaCpmCounts){
          dat$Expression <- dat$Mean.CPM
          if(input$funcMirnaLog2Scale){
            dat$Expression <- dat$Mean.Log2.CPM
          }
        } else{
          dat$Expression <- dat$Mean.Expression
          if(input$funcMirnaLog2Scale){
            dat$Expression <- dat$Mean.Log2.Expression
          }
        }
        label.colors <- tissue_colors
        label.colors[!names(label.colors) %in% tissues()] <- "black"
        label.face <- rep("plain", length(tissue_colors))
        names(label.face) <- names(tissue_colors)
        label.face[names(label.face) %in% tissues()] <- "bold"
        p <- ggplot(filter(dat, miRNA == mirna.id), aes(x = Tissue, y = Expression)) + 
          stat_summary(geom = "bar", aes(fill = Tissue), width = .4) + 
          scale_fill_manual(values = tissue_colors) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, color = label.colors, face = label.face), axis.text = element_text(size = 12), 
                axis.title = element_text(size = 14),
                legend.position = "none")
        return(p)
      })
      tissues <- reactive({
        if(input$categorizeBy == "Function"){return(NULL)}
        tissues <- gsub(", ",'","',circFuncResTable()$Tissue[circFuncResTable()$miRNA == input$selectMirnaForPlotting], fixed = TRUE)
        tissues <- eval(parse(text=paste0('c("', tissues,'")')))
        return(tissues)
      })
    } else if(input$categorizeBy == "Gene"){
      output$geneExpressionTitle <- renderUI({
        HTML(paste("<h3>",input$selectGeneForPlotting,"</h3>"))
      })
      output$circRnaExpressionTitle <- renderUI({
        HTML(paste("<h3>",circForFuncAnalysis(),"</h3>"))
      })
      output$funcCircRnaExpression <- renderPlot({
        withProgress(message = 'Plotting figures...', min = 0, max = 1, value = .4, {
          circ.id <- circSigGenes()$circ.id
          dat <- circ.norm.summary
          if(input$funcCircSrpbmCounts){
            dat$Expression <- dat$Mean.SRPBM
            if(input$funcCircLog2Scale){
              dat$Expression <- dat$Mean.Log2.SRPBM
            }
          } else{
            dat$Expression <- dat$Mean.Expression
            if(input$funcCircLog2Scale){
              dat$Expression <- dat$Mean.Log2.Expression
            } 
          }
          label.colors <- tissue_colors
          label.colors[!names(label.colors) %in% tissues()] <- "black"
          label.face <- rep("plain", length(tissue_colors))
          names(label.face) <- names(tissue_colors)
          label.face[names(label.face) %in% tissues()] <- "bold"
          p <- ggplot(filter(dat, circRNA == circ.id), aes(x = Tissue, y = Expression)) + 
            stat_summary(geom = "bar", aes(fill = Tissue), width = .4) + 
            scale_fill_manual(values = tissue_colors) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, color = label.colors, face = label.face), axis.text = element_text(size = 12), 
                  axis.title = element_text(size = 14),
                  legend.position = "none")
          return(p)
        })
      })
      output$funcGeneExpression <- renderPlot({
        if(input$categorizeBy == "Function"){return(NULL)}
        gene.id <- input$selectGeneForPlotting
        dat <- rna.norm.summary
        if(input$funcGeneCpmCounts){
          dat$Expression <- dat$Mean.CPM
          if(input$funcGeneLog2Scale){
            dat$Expression <- dat$Mean.Log2.CPM
          }
        } else{
          dat$Expression <- dat$Mean.Expression
          if(input$funcGeneLog2Scale){
            dat$Expression <- dat$Mean.Log2.Expression
          }
        }
        label.colors <- tissue_colors
        label.colors[!names(label.colors) %in% tissues()] <- "black"
        label.face <- rep("plain", length(tissue_colors))
        names(label.face) <- names(tissue_colors)
        label.face[names(label.face) %in% tissues()] <- "bold"
        p <- ggplot(filter(dat, gene_symbol == gene.id), aes(x = Tissue, y = Expression)) + 
          stat_summary(geom = "bar", aes(fill = Tissue), width = .4) + 
          scale_fill_manual(values = tissue_colors) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, color = label.colors, face = label.face), axis.text = element_text(size = 12), 
                axis.title = element_text(size = 14),
                legend.position = "none")
        return(p)
      })
      tissues <- reactive({
        if(input$categorizeBy == "Function"){return(NULL)}
        tissues <- gsub(", ",'","',circFuncResTable()$Tissue[circFuncResTable()$Genes == input$selectGeneForPlotting], fixed = TRUE)
        tissues <- eval(parse(text=paste0('c("', tissues,'")')))
        return(tissues)
      })
    } else if(input$categorizeBy == "Function"){
      output$selectedMiRnaExpression <- renderUI({
        return(NULL)
      })
    }
  })
  
  output$mirnaGeneNetworkPlot <- renderVisNetwork({
    req(mirnaGeneNetworkData())
    set.seed(1456)
    visNetwork(mirnaGeneNetworkData()$nodes, mirnaGeneNetworkData()$edges) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = list(enabled = TRUE, values = mirnaGeneNetworkData()$nodes$id)) %>%
      visIgraphLayout(layout = "layout_with_fr", weights = mirnaGeneNetworkData()$edges$Score) %>%
      visLegend(position = "right", useGroups = FALSE, addNodes = mirnaGeneNetworkData()$legend_dat, stepY = 75, width = .15)
  })
  
  observeEvent(input$clicked_text_function,{
    updateTabItems(session, "functionalEnrichmentTabs", "Figure")
  })
  
  observe({
    if(input$rnaNormCounts){
      updateCheckboxInput(session, "rnaCpmCounts", "CPM", FALSE)
    } else{
      updateCheckboxInput(session, "rnaCpmCounts", "CPM", TRUE)
    }
  })
  
  observe({
    if(input$rnaCpmCounts){
      updateCheckboxInput(session, "rnaNormCounts", "DESeq2 normalized", FALSE)
    } else{
      updateCheckboxInput(session, "rnaNormCounts", "DESeq2 normalized", TRUE)
    }
  })
  
  observe({
    if(input$mirnaNormCounts){
      updateCheckboxInput(session, "mirnaCpmCounts", "CPM", FALSE)
    } else{
      updateCheckboxInput(session, "mirnaCpmCounts", "CPM", TRUE)
    }
  })
  
  observe({
    if(input$mirnaCpmCounts){
      updateCheckboxInput(session, "mirnaNormCounts", "DESeq2 normalized", FALSE)
    } else{
      updateCheckboxInput(session, "mirnaNormCounts", "DESeq2 normalized", TRUE)
    }
  })
  
  observe({
    if(input$circNormCounts){
      updateCheckboxInput(session, "circSrpbmCounts", "SRPBM", FALSE)
    } else{
      updateCheckboxInput(session, "circSrpbmCounts", "SRPBM", TRUE)
    }
  })
  
  observe({
    if(input$circSrpbmCounts){
      updateCheckboxInput(session, "circNormCounts", "DESeq2 normalized", FALSE)
    } else{
      updateCheckboxInput(session, "circNormCounts", "DESeq2 normalized", TRUE)
    }
  })
  
  observe({
    if(input$cpmCounts){
      updateCheckboxInput(session, "normCounts", "DESeq2 normalized", FALSE)
    } else{
      updateCheckboxInput(session, "normCounts", "DESeq2 normalized", TRUE)
    }
  })
  
  observe({
    if(input$normCounts){
      updateCheckboxInput(session, "cpmCounts", "CPM/SRPBM", FALSE)
    } else{
      updateCheckboxInput(session, "cpmCounts", "CPM/SRPBM", TRUE)
    }
  })
  
  
  observe({
    if(input$mirnaNetworkCircSrpbmCounts){
      updateCheckboxInput(session, "mirnaNetworkCircNormCounts", "DESeq2 normalized", FALSE)
    } else{
      updateCheckboxInput(session, "mirnaNetworkCircNormCounts", "DESeq2 normalized", TRUE)
    }
  })
  
  observe({
    if(input$mirnaNetworkCircNormCounts){
      updateCheckboxInput(session, "mirnaNetworkCircSrpbmCounts", "SRPBM", FALSE)
    } else{
      updateCheckboxInput(session, "mirnaNetworkCircSrpbmCounts", "SRPBM", TRUE)
    }
  })
  
  observe({
    if(input$circNetworkTableNormCounts){
      updateCheckboxInput(session, "circNetworkTableCpmCounts", "CPM", FALSE)
    } else{
      updateCheckboxInput(session, "circNetworkTableCpmCounts", "CPM", TRUE)
    }
  })
  
  observe({
    if(input$circNetworkTableCpmCounts){
      updateCheckboxInput(session, "circNetworkTableNormCounts",  "DESeq2 Normalized", FALSE)
    } else{
      updateCheckboxInput(session, "circNetworkTableNormCounts", "DESeq2 Normalized", TRUE)
    }
  })
  
  
  observe({
    req(input$funcCircNormCounts)
    if(input$funcCircNormCounts){
      updateCheckboxInput(session, "funcCircSrpbmCounts", "SRPBM", FALSE)
    } else{
      updateCheckboxInput(session, "funcCircSrpbmCounts", "SRPBM", TRUE)
    }
  })
  
  observe({
    req(input$funcCircSrpbmCounts)
    if(input$funcCircSrpbmCounts){
      updateCheckboxInput(session, "funcCircNormCounts", "DESeq2 Normalized", FALSE)
    } else{
      updateCheckboxInput(session, "funcCircNormCounts", "DESeq2 Normalized", TRUE)
    }
  })
  
  
  observe({
    req(input$funcMirnaNormCounts)
    if(input$funcMirnaNormCounts){
      updateCheckboxInput(session, "funcMirnaCpmCounts", "CPM", FALSE)
    } else{
      updateCheckboxInput(session, "funcMirnaCpmCounts", "CPM", TRUE)
    }
  })
  
  observe({
    req(input$funcMirnaCpmCounts)
    if(input$funcMirnaCpmCounts){
      updateCheckboxInput(session, "funcMirnaNormCounts", "DESeq2 Normalized", FALSE)
    } else{
      updateCheckboxInput(session, "funcMirnaNormCounts", "DESeq2 Normalized", TRUE)
    }
  })
  
  
  observe({
    req(input$funcGeneNormCounts)
    if(input$funcGeneNormCounts){
      updateCheckboxInput(session, "funcGeneCpmCounts", "CPM", FALSE)
    } else{
      updateCheckboxInput(session, "funcGeneCpmCounts", "CPM", TRUE)
    }
  })
  
  observe({
    req(input$funcGeneCpmCounts)
    if(input$funcGeneCpmCounts){
      updateCheckboxInput(session, "funcGeneNormCounts", "DESeq2 Normalized", FALSE)
    } else{
      updateCheckboxInput(session, "funcGeneNormCounts", "DESeq2 Normalized", TRUE)
    }
  })
  
  
  mirnaGeneNetworkData <- reactive({
    req(circFuncResTable())
    if(input$categorizeBy == "Function"){
      genes <- gsub(", ",'","',circFuncResTable()$Genes[circFuncResTable()$Function == input$selectFunctionForPlotting], fixed = TRUE)
      genes <- eval(parse(text=paste0('c("', genes,'")')))
      mirnas <- gsub(", ",'","',circFuncResTable()$miRNA[circFuncResTable()$Function == input$selectFunctionForPlotting], fixed = TRUE)
      mirnas <- eval(parse(text=paste0('c("', mirnas,'")')))
      circ.gene.scores <- circSigGenes()$scores
      circ.gene.scores$from <- circSigGenes()$circ.id
      circ.gene.scores$to <- rownames(circ.gene.scores)
      circ.gene.scores$Score1 <- NA
      circ.gene.scores$Score2 <- NA
      circ.gene.scores$color <- "grey"
      circ.gene.scores <- circ.gene.scores[,c("from","to","Score1","Score2","Score","color")]
      edges_circ_mir <- circFunResTableLong() %>% 
        filter(Function == input$selectFunctionForPlotting) %>%
        mutate(circRNA = circSigGenes()$circ.id, color = 'grey') %>%
        select(circRNA, miRNA, MRE, MRE.per.kb, Score, color) %>%
        rename(from = circRNA, to = miRNA, Score1 = MRE, Score2 = MRE.per.kb)
      edges_mir_gene <- circFunResTableLong() %>%
        filter(Function == input$selectFunctionForPlotting) %>%
        mutate(color = 'grey') %>%
        select(miRNA, Genes, 'miR-gene MRE', mirTarBase,'miR-gene Score', color) %>%
        rename(from = miRNA, to = Genes, Score1 = 'miR-gene MRE', Score2 = mirTarBase, Score = 'miR-gene Score')
      circ.gene.scores <- circ.gene.scores %>%
        filter(to %in% edges_mir_gene$to)
      edges <- rbind(edges_circ_mir, edges_mir_gene, circ.gene.scores)
      nodes <- tibble(id = c(circSigGenes()$circ.id, edges_mir_gene$from, edges_mir_gene$to), 
                      shape = c("triangle",rep("square",nrow(edges_mir_gene)), rep("dot", nrow(edges_mir_gene))),
                      type = c("circRNA", rep("miRNA", nrow(edges_mir_gene)), rep("Gene", nrow(edges_mir_gene)))) %>%
        distinct() %>%
        mutate(label = id, color = ifelse(shape == "square", "orange", 
                                          ifelse(shape == "triangle", circNetworkData()$nodes$color[circNetworkData()$nodes$type == "circRNA"], "turquoise")),
               size = 20, font.size = ifelse(length(unique(id)) > 30, 30, length(unique(id))), Score = NA)
      nodes$Score[nodes$type == "Gene"] <- circ.gene.scores$Score[match(nodes$id[nodes$type == "Gene"], circ.gene.scores$to, nomatch = 0)]
      nodes <- nodes %>%
        left_join(circNetworkData()$nodes, by = "id") %>%
        select(id, shape.x, type.x, label.x, color.x, size.x, font.size.x, Score.x, group, MRE, MRE.per.kb, Score.y) %>%
        mutate(Score = ifelse(type.x == "Gene", Score.x, Score.y)) %>%
        select(id, shape.x, type.x, label.x, color.x, size.x, font.size.x, group, MRE, MRE.per.kb, Score) %>%
        rename(shape = shape.x, size = size.x, type = type.x, label = label.x, color = color.x, font.size = font.size.x) %>%
        mutate(title = ifelse(type == "circRNA", paste0("<p>",id,"<br>","<strong>circBase ID:</strong> ", annots$circBase_ID[match(id, annots$circRNA)], "<br>","<strong>circAtlas ID:</strong> ",
                                                        annots$circAtlas_ID[match(id, annots$circRNA)], "<br>", "<strong>ciri2 ID:</strong> ", annots$ciri2_ID[match(id, annots$circRNA)],
                                                        "<br>","<strong>Gene Symbol:</strong> ", annots$gene_symbol[match(id, annots$circRNA)],
                                                        "<br>", "<strong>Tissue:</strong> ", group),
                              ifelse(type == "miRNA", paste0("<p>",label,"<br>","Circ-miR Score: ", round(Score, 4),"<br>", "MRE: ", MRE, "<br>", "MRE per KB: ", round(MRE.per.kb, 4)),
                                     paste0("<p>",label,"<br>","Circ-gene Score: ", round(Score, 4),"</p>"))), size = ifelse(type != "circRNA", Score*size, size))
      if(!is.null(filterNodes()$mirnas) | !is.null(filterNodes()$genes)){
        if(!is.null(filterNodes()$mirnas) & !is.null(filterNodes()$genes)){
          nodes <- nodes[-which(nodes$id %in% filterNodes()$mirnas | nodes$id %in% filterNodes()$genes),]
          edges <- edges[-which(edges$from %in% filterNodes()$mirnas | edges$to %in% filterNodes()$mirnas | edges$from %in% filterNodes()$genes | edges$to %in% filterNodes()$genes),] 
        } else{
          if(!is.null(filterNodes()$mirnas)){
            nodes <- nodes[-which(nodes$id %in% filterNodes()$mirnas),]
            edges <- edges[-which(edges$from %in% filterNodes()$mirnas | edges$to %in% filterNodes()$mirnas),] 
          }
          if(!is.null(filterNodes()$genes)){
            nodes <- nodes[-which(nodes$id %in% filterNodes()$genes),]
            edges <- edges[-which(edges$from %in% filterNodes()$genes | edges$to %in% filterNodes()$genes),] 
          }
        }
      }
      nodes$font.size <- ifelse(nrow(nodes) > 30, 30, nrow(nodes))
      labels <- c("circRNA","miRNA","Gene") 
      color.background <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise")
      color.border <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise")
      shape <- c("triangle","square","dot")
      legend_dat <- data.frame(label = labels,
                               color.background = color.background,
                               color.border = color.border)
      legend_dat$shape <- shape
      legend_dat$size <- 10
      vis_obj <- list(nodes = nodes, edges = edges, legend_dat = legend_dat)
    } else if(input$categorizeBy == "miRNA"){
      terms <- gsub("; ",'","',circFuncResTable()$Function[circFuncResTable()$miRNA == input$selectMirnaForPlotting], fixed = TRUE)
      terms <- eval(parse(text=paste0('c("', terms,'")')))
      genes <- gsub(", ",'","',circFuncResTable()$Genes[circFuncResTable()$miRNA == input$selectMirnaForPlotting], fixed = TRUE)
      genes <- eval(parse(text=paste0('c("', genes,'")')))
      circ.gene.scores <- circSigGenes()$scores
      circ.gene.scores$from <- circSigGenes()$circ.id
      circ.gene.scores$to <- rownames(circ.gene.scores)
      circ.gene.scores$Score1 <- NA
      circ.gene.scores$Score2 <- NA
      circ.gene.scores <- circ.gene.scores[,c("from","to","Score1","Score2","Score")]
      edges_circ_mir <- circFunResTableLong() %>%
        filter(miRNA %in% input$selectMirnaForPlotting) %>%
        mutate(circRNA = circSigGenes()$circ.id, color = 'grey') %>%
        select(circRNA, miRNA, MRE, MRE.per.kb, Score, color) %>%
        rename(from = circRNA, to = miRNA, Score1 = MRE, Score2 = MRE.per.kb)
      edges_circ_term <- circFunResTableLong() %>%
        filter(miRNA == input$selectMirnaForPlotting) %>%
        mutate(circRNA = circSigGenes()$circ.id, color = 'grey', Score1 = NA, Score2 = NA, Score = 1) %>%
        select(circRNA, Function, Score1, Score2, Score, color) %>%
        rename(from = circRNA, to = Function)
      edges_mir_gene <- circFunResTableLong() %>%
        filter(miRNA %in% input$selectMirnaForPlotting) %>%
        mutate(color = 'grey') %>%
        select(miRNA, Genes, 'miR-gene MRE', mirTarBase,'miR-gene Score', color) %>%
        rename(from = miRNA, to = Genes, Score1 = 'miR-gene MRE', Score2 = mirTarBase, Score = 'miR-gene Score')
      edges_term_gene <- circFunResTableLong() %>%
        filter(miRNA %in% input$selectMirnaForPlotting) %>%
        mutate(color = 'grey', Score1 = NA, Score2 = NA, Score = 1) %>%
        select(Function, Genes, Score1, Score2, Score, color) %>%
        rename(from = Function, to = Genes)
      edges_mir_term <- circFunResTableLong() %>%
        filter(miRNA == input$selectMirnaForPlotting) %>%
        mutate(color = 'grey', Score1 = NA, Score2 = NA, Score = 1) %>%
        select(miRNA, Function, Score1, Score2, Score, color) %>%
        rename(from = miRNA, to = Function)
      edges_circ_gene <- circ.gene.scores %>%
        filter(to %in% edges_mir_gene$to) %>%
        mutate(color = 'grey')
      edges <- rbind(edges_circ_mir, edges_circ_gene, edges_circ_term, edges_mir_gene, edges_term_gene, edges_mir_term)
      nodes <- tibble(id = c(circSigGenes()$circ.id, input$selectMirnaForPlotting, edges_mir_gene$to, edges_mir_term$to), 
                      shape = c("star","triangle",rep("dot",nrow(edges_mir_gene)), rep("square", nrow(edges_mir_term))),
                      type = c("circRNA","miRNA", rep("Gene", nrow(edges_mir_gene)), rep("Function", nrow(edges_mir_term)))) %>%
        distinct() %>%
        mutate(label = id, color = ifelse(shape == "triangle", "orange", 
                                          ifelse(shape == "square","navy", 
                                                 ifelse(shape == "dot", "turquoise", 
                                                        circNetworkData()$nodes$color[circNetworkData()$nodes$id %in% circSigGenes()$circ.id]))),
               size = 20, font.size = ifelse(length(unique(id)) > 30, 30, length(unique(id))))
      nodes$MRE <- nodes$MRE.per.kb <- nodes$mirTarBase <- nodes$Score <- nodes$pvalue <- NA
      nodes$MRE[nodes$shape %in% "triangle"] <- circNetworkData()$nodes$MRE[circNetworkData()$nodes$id %in% nodes$id[nodes$shape %in% "triangle"]]
      nodes$MRE[nodes$type %in% "Gene"] <- edges_mir_gene$Score1[match(nodes$id[nodes$type %in% "Gene"], edges_mir_gene$to, nomatch = 0)]
      nodes$mirTarBase[nodes$type %in% "Gene"] <- edges_mir_gene$Score2[match(nodes$id[nodes$type %in% "Gene"], edges_mir_gene$to, nomatch = 0)]
      nodes$Score[nodes$shape %in% "triangle"] <- circNetworkData()$nodes$Score[circNetworkData()$nodes$id %in% nodes$id[nodes$shape %in% "triangle"]]
      nodes$Score[nodes$type %in% "Gene"] <- edges_mir_gene$Score[match(nodes$id[nodes$type %in% "Gene"], edges_mir_gene$to, nomatch = 0)]
      nodes$pvalue[nodes$type %in% "Gene"] <- circFunResTableLong()$Gene.pvalue[match(nodes$id[nodes$type %in% "Gene"], circFunResTableLong()$Genes, nomatch = 0)]
      nodes <- nodes %>%
        mutate(title = ifelse(type == "miRNA", paste0("<p>",label,"<br>","Circ-miR Score: ", round(Score, 4),"<br>", "MRE: ", MRE, "<br>", "MRE per KB: ", round(MRE.per.kb, 4)),
                              ifelse(shape == "Gene", paste0("<p>",label,"<br>","miR-gene Score: ", round(Score, 4),"<br>", "MRE: ", MRE,
                                                            "<br>", "mirTarBase: ", mirTarBase, "<br>", "p-value: ", pvalue, "<br>"),
                                     paste0("<p>",label))))
      
      if(!is.null(filterNodes()$terms) | !is.null(filterNodes()$genes)){
        if(!is.null(filterNodes()$terms) & !is.null(filterNodes()$genes)){
          nodes <- nodes[-which(nodes$id %in% filterNodes()$terms | nodes$id %in% filterNodes()$genes),]
          edges <- edges[-which(edges$from %in% filterNodes()$terms | edges$to %in% filterNodes()$terms | edges$from %in% filterNodes()$genes | edges$to %in% filterNodes()$genes),] 
        } else{
          if(!is.null(filterNodes()$terms)){
            nodes <- nodes[-which(nodes$id %in% filterNodes()$terms),]
            edges <- edges[-which(edges$from %in% filterNodes()$terms | edges$to %in% filterNodes()$terms),] 
          }
          if(!is.null(filterNodes()$genes)){
            nodes <- nodes[-which(nodes$id %in% filterNodes()$genes),]
            edges <- edges[-which(edges$from %in% filterNodes()$genes | edges$to %in% filterNodes()$genes),] 
          }
        }
      }
      nodes$font.size <- ifelse(nrow(nodes) > 30, 30, nrow(nodes))
      if(length(which(nodes$shape %in% "square")) <= 10){
        nodes$color[which(nodes$shape %in% "square")] <- pathway_colors[1:length(which(nodes$shape %in% "square"))] 
      }
      if(length(unique(nodes$id[nodes$shape == "square"])) <= 10){
        labels <- c("circRNA","miRNA","Gene",nodes$id[nodes$shape == "square"])
        color.background <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise", nodes$color[nodes$shape == "square"])
        color.border <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise", nodes$color[nodes$shape == "square"])
        shape <- c("star","triangle","dot",rep("square", length(unique(nodes$id[nodes$shape == "square"]))))
      } else{
        labels <- c("circRNA","miRNA","Gene","Function") 
        color.background <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise", nodes$color[nodes$shape == "square"][1])
        color.border <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise", nodes$color[nodes$shape == "square"][1])
        shape <- c("star","triangle","dot","square")
      }
      legend_dat <- data.frame(label = labels,
                               color.background = color.background,
                               color.border = color.border)
      legend_dat$shape <- shape
      legend_dat$size <- 10
      vis_obj <- list(nodes = nodes, edges = edges, legend_dat = legend_dat) 
    } else{
      terms <- gsub("; ",'","',circFuncResTable()$Function[circFuncResTable()$Genes == input$selectGeneForPlotting], fixed = TRUE)
      terms <- eval(parse(text=paste0('c("', terms,'")')))
      mirnas <- gsub(", ",'","',circFuncResTable()$miRNA[circFuncResTable()$Genes == input$selectGeneForPlotting], fixed = TRUE)
      mirnas <- eval(parse(text=paste0('c("', mirnas,'")')))
      circ.gene.scores <- circSigGenes()$scores
      circ.gene.scores$from <- circSigGenes()$circ.id
      circ.gene.scores$to <- rownames(circ.gene.scores)
      circ.gene.scores$Score1 <- NA
      circ.gene.scores$Score2 <- NA
      circ.gene.scores <- circ.gene.scores[,c("from","to","Score1","Score2","Score")]
      edges_circ_mir <- circFunResTableLong() %>%
        filter(Genes %in% input$selectGeneForPlotting) %>%
        mutate(circRNA = circSigGenes()$circ.id, color = 'grey') %>%
        select(circRNA, miRNA, MRE, MRE.per.kb, Score, color) %>%
        rename(from = circRNA, to = miRNA, Score1 = MRE, Score2 = MRE.per.kb)
      edges_circ_term <- circFunResTableLong() %>%
        filter(Genes %in% input$selectGeneForPlotting) %>%
        mutate(circRNA = circSigGenes()$circ.id, color = 'grey', Score1 = NA, Score2 = NA, Score = 1) %>%
        select(circRNA, Function, Score1, Score2, Score, color) %>%
        rename(from = circRNA, to = Function)
      edges_gene_mir <- circFunResTableLong() %>%
        filter(Genes %in% input$selectGeneForPlotting) %>%
        mutate(color = 'grey') %>%
        select(Genes, miRNA, 'miR-gene MRE', mirTarBase,'miR-gene Score', color) %>%
        rename(from = Genes, to = miRNA, Score1 = 'miR-gene MRE', Score2 = mirTarBase, Score = 'miR-gene Score')
      edges_gene_term <- circFunResTableLong() %>%
        filter(Genes %in% input$selectGeneForPlotting) %>%
        mutate(color = 'grey', Score1 = NA, Score2 = NA, Score = 1) %>%
        select(Genes, Function, Score1, Score2, Score, color) %>%
        rename(from = Genes, to = Function)
      edges_mir_term <- circFunResTableLong() %>%
        filter(Genes %in% input$selectGeneForPlotting) %>%
        mutate(color = 'grey', Score1 = NA, Score2 = NA, Score = 1) %>%
        select(miRNA, Function, Score1, Score2, Score, color) %>%
        rename(from = miRNA, to = Function)
      edges_circ_gene <- circ.gene.scores %>%
        filter(to %in% edges_gene_mir$from) %>%
        mutate(color = 'grey')
      edges <- rbind(edges_circ_mir,edges_circ_gene,edges_circ_term,edges_gene_mir, edges_gene_term, edges_mir_term)
      nodes <- tibble(id = c(circSigGenes()$circ.id, input$selectGeneForPlotting, edges_gene_mir$to, edges_mir_term$to), 
                      shape = c("star","triangle",rep("dot",nrow(edges_gene_mir)), rep("square", nrow(edges_mir_term))),
                      type = c("circRNA","Gene", rep("miRNA", nrow(edges_gene_mir)), rep("Function", nrow(edges_mir_term)))) %>%
        distinct() %>%
        mutate(label = id, color = ifelse(shape == "triangle", "orange", 
                                          ifelse(shape == "square","navy", 
                                                 ifelse(shape == "dot", "turquoise", 
                                                        circNetworkData()$nodes$color[circNetworkData()$nodes$id %in% circSigGenes()$circ.id]))),
               size = 20, font.size = ifelse(length(unique(id)) > 30, 30, length(unique(id))))
      nodes$MRE.gene <- nodes$MRE.circ <- nodes$MRE.per.kb <- nodes$mirTarBase <- nodes$Score.gene <- nodes$Score.circ <- nodes$Score <- nodes$pvalue <- NA
      nodes$MRE.gene[nodes$shape %in% "dot"] <- edges_gene_mir$Score1[match(nodes$id[nodes$shape %in% "dot"], edges_gene_mir$to, nomatch = 0)]
      nodes$mirTarBase[nodes$shape %in% "dot"] <- edges_gene_mir$Score2[match(nodes$id[nodes$shape %in% "dot"], edges_gene_mir$to, nomatch = 0)]
      nodes$Score.gene[nodes$shape %in% "dot"] <- edges_gene_mir$Score[match(nodes$id[nodes$shape %in% "dot"], edges_gene_mir$to, nomatch = 0)]
      nodes$MRE.circ[nodes$shape %in% "dot"] <- circNetworkData()$nodes$MRE[match(nodes$id[nodes$shape %in% "dot"],circNetworkData()$nodes$id, nomatch = 0)]
      nodes$MRE.per.kb[nodes$shape %in% "dot"] <- circNetworkData()$nodes$MRE.per.kb[match(nodes$id[nodes$shape %in% "dot"],circNetworkData()$nodes$id, nomatch = 0)]
      nodes$Score.circ[nodes$shape %in% "dot"] <- circNetworkData()$nodes$Score[match(nodes$id[nodes$shape %in% "dot"],circNetworkData()$nodes$id, nomatch = 0)]
      nodes$pvalue[nodes$shape %in% "triangle"] <- circFunResTableLong()$Gene.pvalue[circFunResTableLong()$Genes %in% input$selectGeneForPlotting]
      nodes <- nodes %>%
        mutate(title = ifelse(shape == "dot", paste0("<p>",label,"<br>","Circ-miR Score: ", round(Score.circ, 4),"<br>", "Circ-miR MRE: ", MRE.circ, 
                                                     "<br>", "MRE per KB: ", round(MRE.per.kb, 4), "<br>", "miR-gene Score: ", round(Score.gene, 4), "<br>",
                                                     "miR-gene MRE: ", MRE.gene, "<br>", "mirTarBase: ", mirTarBase),
                              ifelse(shape == "triangle", paste0("<p>",label,"<br>","p-value: ", pvalue, "<br>"),
                                     paste0("<p>",label))))
      
      if(!is.null(filterNodes()$terms) | !is.null(filterNodes()$mirnas)){
        if(!is.null(filterNodes()$terms) & !is.null(filterNodes()$mirnas)){
          nodes <- nodes[-which(nodes$id %in% filterNodes()$terms | nodes$id %in% filterNodes()$mirnas),]
          edges <- edges[-which(edges$to %in% filterNodes()$terms | edges$from %in% filterNodes()$mirnas | edges$to %in% filterNodes()$mirnas),] 
        } else{
          if(!is.null(filterNodes()$terms)){
            nodes <- nodes[-which(nodes$id %in% filterNodes()$terms),]
            edges <- edges[-which(edges$from %in% filterNodes()$terms | edges$to %in% filterNodes()$terms),] 
          }
          if(!is.null(filterNodes()$mirnas)){
            nodes <- nodes[-which(nodes$id %in% filterNodes()$mirnas),]
            edges <- edges[-which(edges$from %in% filterNodes()$mirnas | edges$to %in% filterNodes()$mirnas),] 
          }
        }
      }
      nodes$font.size <- ifelse(nrow(nodes) > 30, 30, nrow(nodes))
      if(length(which(nodes$shape %in% "square")) <= 10){
        nodes$color[which(nodes$shape %in% "square")] <- pathway_colors[1:length(which(nodes$shape %in% "square"))] 
      }
      if(length(unique(nodes$id[nodes$shape == "square"])) <= 10){
        labels <- c("circRNA","Gene","miRNA",nodes$id[nodes$shape == "square"])
        color.background <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise", nodes$color[nodes$shape == "square"])
        color.border <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise", nodes$color[nodes$shape == "square"])
        shape <- c("star","triangle","dot",rep("square", length(unique(nodes$id[nodes$shape == "square"]))))
      } else{
        labels <- c("circRNA","Gene","miRNA","Function") 
        color.background <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise", nodes$color[nodes$shape == "square"][1])
        color.border <- c(nodes$color[nodes$type == "circRNA"], "orange", "turquoise", nodes$color[nodes$shape == "square"][1])
        shape <- c("star","triangle","dot","square")
      }
      legend_dat <- data.frame(label = labels,
                               color.background = color.background,
                               color.border = color.border)
      legend_dat$shape <- shape
      legend_dat$size <- 10
      vis_obj <- list(nodes = nodes, edges = edges, legend_dat = legend_dat)
    }
    return(vis_obj)
  })
}
