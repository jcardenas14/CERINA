ui <- div(id = "main_content",
          dashboardPage(skin = "blue",
                        header <- dashboardHeader(title = strong("CERINA"),
                                                  dropdownMenuOutput("messageMenu")),
                        
                        sidebar <- dashboardSidebar(
                          sidebarMenu(
                            menuItem(strong("Documentation"), icon = icon("book"), tabName = "documentation"),
                            menuItem(strong("Data Exploration"), icon = icon("line-chart"), tabName = "dataExploration"),
                            menuItem(strong("miRNA-circRNA Network"), icon = icon("project-diagram"), tabName = "networkAnalysis"),
                            menuItem(strong("Functional Enrichment"), icon = icon("chart-bar"), tabName = "functionalEnrichment")
                          )
                        ),
                        
                        body <-  dashboardBody(
                          includeScript("helper_functions.js"),
                          useShinyjs(),
                          tags$style(HTML(".content-wrapper, .right-side {background: #ffff; overflow-y: auto}
                    .skin-blue .main-header .navbar {background-color: #241c49}
                    .skin-blue .left-side, .skin-blue .main-sidebar, .skin-blue .wrapper {background-color: #2c343a}
                    .skin-blue .main-header .logo {background-color: #241c49; color: #ffff}
                    .skin-blue .main-header .logo:hover {background-color: #ffff; color: #241c49}
                    .skin-blue .main-header .navbar .sidebar-toggle:hover {background-color: #ffff; color: #241c49}
                    .skin-blue .sidebar-menu {color: #ffff}
                    .skin-blue .sidebar-menu>li.active>a {background: #767381; border-left-color: #1a162d; color: #ffff}
                    .skin-blue .sidebar-menu>li:hover>a {background: #ffff; border-left-color: #1a162d; color: #241c49}
                    .box.box-solid.box-primary {border: 1px solid #0000000d}
                    .box.box-primary {border-top-color: #241c49}
                    .box.box-solid.box-primary>.box-header {background: #241c49}
                    .nav-tabs-custom>.nav-tabs>li.active {border-top-color: #241c49}
                    #unsupervisedTabs.nav.nav-tabs.shiny-tab-input.shiny-bound-input>.nav-tabs-custom>.nav-tabs>li>a, 
                    .nav-tabs-custom>.nav-tabs>li>a:hover {background-color: #7d7799; color: #ffff}
                    #rnaDendrogram.shiny-plot-output.shiny-bound-output {height: 100% !important}
                    #mirnaDendrogram.shiny-plot-output.shiny-bound-output {height: 100% !important}
                    #circDendrogram.shiny-plot-output.shiny-bound-output {height: 100% !important}
                    #rnaExpressionProfilePlot.shiny-plot-output.shiny-bound-output {height: auto !important}
                    #mirnaExpressionProfilePlot.shiny-plot-output.shiny-bound-output {height: auto !important}
                    #circExpressionProfilePlot.shiny-plot-output.shiny-bound-output {height: auto !important}
                    #correlationScatterPlot.shiny-plot-output.shiny-bound-output {height: auto !important}
                    #networkBox {height: auto !important}
                    .box-body {width: 100%}
                    .shiny-input-container:not(.shiny-input-container-inline) {width: 100%}
                    .btn-default {background-color: #4f4c4c; opacity: .8; color: #ffff}
                    .btn-default:hover {background-color: #4f4c4c; opacity: .8; color: #ffff}
                    .btn-default:active {background-color: #4f4c4c; opacity: .8; color: #ffff}
                    .btn-default:focus {background-color: #4f4c4c; opacity: .8; color: #ffff}
                    .btn-default:active:focus {background-color: #4f4c4c; opacity: .8; color: #ffff}
                    .btn-default:active:hover {background-color: #4f4c4c; opacity: .8; color: #ffff}
                    .btn {border-radius: .5rem; border-color: #4f4c4c}
                    .btn.btn-mini {border-color: white;}
                    .fa-info-circle {background-color: white; color: black;}
                    .panel panel-default {aria-expanded = 'true'}
                    .modal-dialog{margin: 0px; padding-left: 15%; padding-right: 15%; padding-top: 5%; width: 100%}
                    .shiny-notification{position: fixed; top: 90%; left: 80%; right: 2%; background: #241c49; color: white}")),
                          tabItems(
                            tabItem("documentation",
                                    tabBox(title = "", id = "documentationTabs", width = 12,
                                           tabPanel(strong("Welcome to CERINA"), icon = icon("info"),
                                                    includeHTML("welcome.html")
                                           ),
                                           tabPanel(strong("Data Exploration"), icon = icon("line-chart"),
                                                    includeHTML("dataExploration.html")
                                           ),
                                           tabPanel(strong("miRNA-circRNA Network"), icon = icon("project-diagram"),
                                                    includeHTML("mirnaCircNetwork.html")
                                           ),
                                           tabPanel(strong("Functional Enrichment"), icon = icon("chart-bar"),
                                                    includeHTML("functionalEnrichment.html")
                                           )
                                    )
                            ),
                            tabItem(tabName = "dataExploration",
                                    fluidRow(
                                      box(title = strong("Hierarchical Clustering"), width = 6, status = "primary", solidHeader = TRUE,
                                          fluidRow(
                                            tabBox(title = "", id = "clusteringTabs", width = 12,
                                                   tabPanel("RNA", 
                                                            plotOutput("rnaDendrogram")
                                                   ),
                                                   tabPanel("miRNA", 
                                                            plotOutput("mirnaDendrogram")
                                                   ),
                                                   tabPanel("circRNA", 
                                                            plotOutput("circDendrogram")
                                                   )
                                            )
                                          )
                                      ),
                                      box(title = strong("Expression Profile"), width = 6, status = "primary", solidHeader = TRUE,
                                          fluidRow(
                                            tabBox(title = "", id = "expressionProfileTabs", width = 12,
                                                   tabPanel("RNA",
                                                            tipify(selectizeInput("selectGene", "Enter Gene", 
                                                                           selected = NULL, multiple = TRUE, choices = NULL,
                                                                           options = list(maxItems = 1, maxOptions = 500), width = "50%"),
                                                                   title = "Ensembl ID or gene symbol", placement = "top"),
                                                            div(class = "row",
                                                                div(class = "col-xs-2",
                                                                    checkboxInput("rnaCpmCounts", "CPM", TRUE)
                                                                ),
                                                                div(class = "col-xs-6",
                                                                    checkboxInput("rnaNormCounts","DESeq2 normalized",FALSE)
                                                                )
                                                            ),
                                                            checkboxInput("rnaLog2Scale", "log2 scale", FALSE, width = '50%'),
                                                            plotOutput("rnaExpressionProfilePlot")
                                                   ),
                                                   tabPanel("miRNA",
                                                            selectizeInput("selectMiRna", "Enter miRNA", 
                                                                           selected = NULL, multiple = TRUE, choices = NULL,
                                                                           options = list(maxItems = 1, maxOptions = 500), width="50%"),
                                                            div(class = "row",
                                                                div(class = "col-xs-2",
                                                                    checkboxInput("mirnaCpmCounts", "CPM", TRUE)
                                                                ),
                                                                div(class = "col-xs-6",
                                                                    checkboxInput("mirnaNormCounts","DESeq2 normalized",FALSE)
                                                                )
                                                            ),
                                                            checkboxInput("mirnaLog2Scale", "log2 scale", FALSE, width = "50%"),
                                                            plotOutput("mirnaExpressionProfilePlot")
                                                   ),
                                                   tabPanel("circRNA",
                                                            div(class = "row",
                                                                div(class = "container-fluid col-xs-6",
                                                                    tipify(selectizeInput("selectCircRna", "Enter circRNA", 
                                                                                   selected = NULL, multiple = TRUE, choices = NULL,
                                                                                   options = list(maxItems = 1, maxOptions = 500)),
                                                                    title = "e.g. CDR1-AS; chrX:139865339-139866824; hsa_circ_0001946; hsa-intergenic_004850",
                                                                    placement = "top")
                                                                ),
                                                                div(class = "container-fluid col-xs-6",
                                                                    uiOutput("multipleMapping")      
                                                                )
                                                            ),
                                                            div(class = "row",
                                                                div(class = "col-xs-2",
                                                                    checkboxInput("circSrpbmCounts", "SRPBM", TRUE)
                                                                ),
                                                                div(class = "col-xs-6",
                                                                    checkboxInput("circNormCounts","DESeq2 normalized",FALSE)
                                                                )
                                                            ),
                                                            checkboxInput("circLog2Scale", "log2 scale", FALSE, width = "50%"),
                                                            plotOutput("circExpressionProfilePlot")
                                                   )
                                            )
                                          )
                                      )
                                    ),
                                    fluidRow(
                                      box(title = strong("Correlation Analysis"), width = 12, status = "primary", solidHeader = TRUE,
                                          fluidRow(
                                            div(class = "container-fluid col-xs-3",
                                                selectizeInput("xaxis", "x-axis variable", 
                                                               selected = NULL, multiple = TRUE, choices = NULL,
                                                               options = list(maxItems = 1, maxOptions = 500)),
                                                selectizeInput("yaxis", "y-axis variable", 
                                                               selected = NULL, multiple = TRUE, choices = NULL,
                                                               options = list(maxItems = 1, maxOptions = 500)),
                                                div(style="display:inline-block",checkboxInput("cpmCounts", "CPM/SRPBM", TRUE)),
                                                div(style="display:inline-block", 
                                                    infoPopup("","CPM plotted for miRNA and linear RNA. SRPBM used for circRNA",
                                                              placement = "right", trigger = "hover")),
                                                checkboxInput("normCounts", "DESeq2 normalized", FALSE)
                                            ),
                                            div(class = "container-fluid col-xs-9",
                                                plotOutput("correlationScatterPlot")  
                                            )
                                          )
                                      )
                                    )
                            ),
                            tabItem("networkAnalysis",
                                    fluidRow(
                                      div(class = "container-fluid",
                                          bsCollapse(multiple = FALSE, open = "networkOptions",
                                                     bsCollapsePanel("Network Filtering Options", value = "networkOptions",
                                                                     div(class = "row",
                                                                         div(class = "col-xs-3",
                                                                             selectizeInput("selectMirnaNetwork", "Enter miRNA", 
                                                                                            selected = NULL, multiple = TRUE, choices = NULL,
                                                                                            options = list(maxItems = 1, maxOptions = 500)),
                                                                             actionButton("plotMiRnaNetwork", "Create Network", icon = icon("marker"), 
                                                                                          style = "background: rgb(153, 50, 204)")
                                                                         ),
                                                                         div(class = "col-xs-3",
                                                                             tipify(numericInput("filterByScore", "Filter by Pareto score", 
                                                                                          min = 0, max = 1, value = .8, step = .1),
                                                                             title = "circRNAs with pareto interaction score >= threshold included", 
                                                                             placement = "top")
                                                                         ),
                                                                         div(class = "col-xs-3",
                                                                             tipify(numericInput("filterByExpression", "Filter by circRNA expression", 
                                                                                          min = 0, value = 1, step = .5),
                                                                                    title = "circRNAs with mean normalized expression >= threshold in at least one tissue included",
                                                                                    placement = "top"),
                                                                             div(class = "row",
                                                                                 div(class = "col-xs-3",
                                                                                     checkboxInput("mirnaNetworkCircSrpbmCounts", "SRPBM", TRUE)
                                                                                 ),
                                                                                 div(class = "col-xs-9",
                                                                                     checkboxInput("mirnaNetworkCircNormCounts","DESeq2 normalized",FALSE)
                                                                                 )
                                                                             )
                                                                         ),
                                                                         div(class = "col-xs-3",
                                                                             tipify(numericInput("circRnaCorrThresh", 
                                                                                          "Filter by circ-circ spearman rho",
                                                                                          min = 0, max = 1, value = .8, step = .05),
                                                                                    title = "Edges drawn between circRNAs with spearman rho >= threshold",
                                                                                    placement = 'top')
                                                                         )
                                                                     ),
                                                                     style = "default") 
                                          )
                                      )
                                    ),
                                    fluidRow(
                                      div(class = "container-fluid",
                                          tabBox(title = "", id = "networkAnalysisTabs", width = 12,
                                                 tabPanel("Network",
                                                          uiOutput("tableViewButton"),
                                                          br(),
                                                          br(),
                                                          withSpinner(visNetworkOutput("mirnaNetworkPlot", height = '50em'))
                                                 ),
                                                 tabPanel("Table",
                                                          fluidRow(
                                                            downloadButton('downloadMirnaNetworkTable', 'Download'),
                                                            br(),
                                                            br(),
                                                            DT::DTOutput("mirnaNetworkTable"),
                                                            uiOutput("circRnaExpressionPlot")
                                                          ) 
                                                 ),
                                                 tabPanel("Heatmap",
                                                          withSpinner(plotlyOutput("plotCircRnaHeatmap", height = "55em", width = "auto"))
                                                 )
                                          )
                                      )
                                    )
                            ),
                            tabItem("functionalEnrichment",
                                    tabBox(title = "", id = "enrichmentWorkflow", width = 12,
                                           tabPanel(strong("miRNA Exploration"), icon = icon("angle-double-right"),
                                                    fluidRow(
                                                      br(),
                                                      div(class = "container-fluid row",
                                                          div(class = "col-xs-3",
                                                              tipify(selectizeInput("circForAnalysis", "Enter circRNA", 
                                                                             selected = NULL, multiple = TRUE, choices = NULL,
                                                                             options = list(maxItems = 1, maxOptions = 500)), 
                                                                     title = "e.g. CDR1-AS; chrX:139865339-139866824; hsa_circ_0001946; hsa-intergenic_004850", 
                                                                     placement = "top")
                                                          ),
                                                          div(class = "col-xs-3",
                                                              uiOutput("circMultipleMapping")
                                                          )
                                                      ),
                                                      div(class = "container-fluid",
                                                          actionButton("viewNetwork", "View Network"),
                                                          br(),
                                                          br(),
                                                          downloadButton('downloadCircMirnaNetworkTable', 'Download'),
                                                          br(),
                                                          div(class = "row",
                                                              div(class = "col-xs-1",
                                                                  checkboxInput("circNetworkTableCpmCounts", strong("CPM"), TRUE)
                                                              ),
                                                              div(class = "col-xs-6",
                                                                  checkboxInput("circNetworkTableNormCounts",strong("DESeq2 normalized"),FALSE)
                                                              )
                                                          ),
                                                          DT::DTOutput("circMirnaNetworkTable"),
                                                          uiOutput("mirnaExpressionPlot"),
                                                          div(class = "col-xs-12",
                                                              bsModal("modalNetwork","Network","viewNetwork", size = "large",
                                                                      visNetworkOutput("circNetworkPlot", height = "50em") 
                                                              ) 
                                                          )
                                                      )
                                                    )
                                           ),
                                           tabPanel(strong("Gene Level Test"), icon = icon("angle-double-right"),
                                                    fluidRow(
                                                      div(class = "col-xs-2 container-fluid", 
                                                          div(class = "container-fluid", style = "border: 1px solid #0000000d",
                                                              HTML('<h4><strong>Gene level test</strong></h4>'),
                                                              hr(),
                                                              htmlOutput("circSelected"),
                                                              htmlOutput("circBaseID"),
                                                              htmlOutput("circAtlasID"),
                                                              htmlOutput("ciri2ID"),
                                                              htmlOutput("parentalGeneSymbol"),
                                                              htmlOutput("numberMirnasSelected"),
                                                              hr(),
                                                              actionButton("runGeneAnalysis", "Run Analysis", icon = icon("play"), style = "background: rgb(153, 50, 204)"),
                                                              br(),
                                                              br()
                                                          )
                                                      ),
                                                      column(width = 10,
                                                             fluidRow(
                                                               box(title = strong("Significant Genes Summary"), width = 12, status = "primary", solidHeader = TRUE,
                                                                   fluidRow(
                                                                     column(width = 5,
                                                                            div(class = "row",
                                                                                div(class = "container-fluid col-xs-6",
                                                                                    downloadButton('downloadNumberSigGeneTable', 'Download Table')
                                                                                ),
                                                                                div(class = "container-fluid col-xs-6",
                                                                                    downloadButton('downloadSigGeneList', 'Download Gene List')
                                                                                )
                                                                            ),
                                                                            br(),
                                                                            DT::DTOutput("numberSigGeneTable")
                                                                     ),
                                                                     column(width = 7,
                                                                            div(style="margin-left:1%",
                                                                                plotOutput("numberSigGenePlot", width = "auto")   
                                                                            )
                                                                     )
                                                                   )
                                                               )
                                                             )
                                                      ),
                                                      box(title = strong("Gene-miRNA Table"), width = 12, status = "primary", solidHeader = TRUE,
                                                          numericInput("geneLevelTestThresh", "Gene p-value threshold", value = .05, 
                                                                       min = 0, max = 1, step = .01, width = '25%'),
                                                          checkboxInput("geneLevelTestFdrAdjusted", "FDR adjusted", FALSE),
                                                          downloadButton('downloadGeneMirnaTable', 'Download Table'),
                                                          br(),
                                                          br(),
                                                          DT::DTOutput("geneMirnaTable")
                                                      ) 
                                                    )
                                           ),
                                           tabPanel(strong('Functional Enrichment'), icon = icon("angle-double-right"), value = 'func_enrich_val',
                                                    fluidRow(
                                                      div(class = "container-fluid",
                                                          bsCollapse(multiple = FALSE, open = "options",
                                                                     bsCollapsePanel("Test Options", value = "options",
                                                                                     div(class = "row",
                                                                                         div(class = "col-xs-3",
                                                                                             selectizeInput("termType", "Select Function type:", choices = c("KEGG","GO BP","GO CC","GO MF"), selected = "KEGG", multiple = FALSE),
                                                                                             actionButton("runFuncAnalysis", "Run Analysis", icon = icon("play"), style = "background: rgb(153, 50, 204)")
                                                                                         ),
                                                                                         div(class = "col-xs-3",
                                                                                             numericInput("genePvalThresh", "Gene p-value threshold", value = .05, min = 0, max = 1, step = .01),
                                                                                             checkboxInput("fdrAdjustedGenes", "FDR adjusted", FALSE)
                                                                                         ),
                                                                                         div(class = "col-xs-3",
                                                                                             numericInput("termPvalThresh", "Function p-value threshold", value = .05, min = 0, max = 1, step = .01),
                                                                                             checkboxInput("fdrAdjustedTerms", "FDR adjusted", FALSE)
                                                                                         )
                                                                                     ),
                                                                                     style = "default") 
                                                          )
                                                      )
                                                    ),
                                                    fluidRow(
                                                      box(title = strong("Results"), width = 12, status = "primary", solidHeader = TRUE,
                                                          tabBox(title = "", id = "functionalEnrichmentTabs", width = 12,
                                                                 tabPanel("Table",
                                                                          selectizeInput("categorizeBy", "Categorize by:", choices = c("Function", "miRNA", "Gene"), 
                                                                                         selected = "Function", multiple = FALSE, width = "25%"),
                                                                          checkboxInput("longForm", "Display table in long form?", FALSE),
                                                                          downloadButton('downloadFunctionalResults', 'Download'),
                                                                          br(),
                                                                          br(),
                                                                          DT::DTOutput("functionalResults")
                                                                 ),
                                                                 tabPanel("Figure",
                                                                          uiOutput("selectFunctionOrMirna"),
                                                                          uiOutput("selectedMiRnaExpression")
                                                                 ) 
                                                          )
                                                      )
                                                    )
                                           )
                                    )
                            )
                          )
                        )
          )
)

