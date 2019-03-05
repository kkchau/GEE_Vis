# TODO: Reorganize panels and controls

library(tidyverse)
library(shiny)
source("build_graph.R")

example_data <- read_csv("switch_enrichments_unordered_moderate.csv")

ui <- fluidPage(
    
    # Application title
    titlePanel("Gene EnrichmEnt VISualizer"),
   
    # Sidebar with a slider input for number of bins 
    fluidRow(
        column(
            3,
            tabsetPanel(
                type = "tabs",
                tabPanel(
                    "Upload gProfileR File",
                    downloadButton(
                        outputId = "downloadExample",
                        label = "Download an example file"
                    ),
                    fileInput(
                        inputId = "enrfile", 
                        label = "Choose Enrichment File",
                        multiple = FALSE,
                        accept = c(
                            "text/csv",
                            "text/comma-separated-values",
                            "text/tab-separated-values",
                            ".csv",
                            ".tsv",
                            ".txt",
                            ".tab"
                        ),
                        placeholder = ".csv, .tsv, .txt"
                    ),
                    uiOutput("PValueColumn"),
                    uiOutput("TermNameColumn"),
                    uiOutput("GeneOverlapColumn"),
                    checkboxInput(
                        inputId = "showNegativeFlag",
                        label = "+/- flags?", 
                        value = FALSE
                    ),
                    conditionalPanel(
                        condition = "input.showNegativeFlag == true",
                        uiOutput("NegativeFlagColumnUI"),
                        uiOutput("PositiveFlagsInput"),
                        uiOutput("NegativeFlagsInput")
                    )
                ),
                tabPanel(
                    "Display Parameters",
                    uiOutput("numDisplay"),
                    numericInput(
                        inputId = "PValueThresh",
                        label = "P-Value Threshold Marker",
                        value = 0.05
                    ),
                    sliderInput(
                        inputId = "fontSizeLab",
                        label = "Label Font Size",
                        min = 5, max = 30, value = 12
                    ),
                    sliderInput(
                        inputId = "fontSizeAx",
                        label = "Axis Font Size",
                        min = 5, max = 30, value = 12
                    ),
                    checkboxGroupInput(
                        inputId = "domains",
                        label = NULL
                    ),
                    uiOutput("termSelect"),
                    uiOutput("downloadFigure")
                ),
                tabPanel(
                    "Download Figure",
                    numericInput(
                        inputId = "figwidth",
                        label = "Figure width (in.)",
                        value = 16
                    ),
                    numericInput(
                        inputId = "figheight",
                        label = "Figure height (in.)",
                        value = 9
                    ),
                    downloadButton(
                        outputId = "figureDownload",
                        label = "Download your figure"
                    )
                )
            )
        ),
        column(
            9,
            tabsetPanel(
                tabPanel(
                    "Enrichment Plot",
                    plotOutput("enrPlot"),
                    textOutput("intersection")
                ),
                tabPanel(
                    "Enrichment Table",
                    DT::dataTableOutput("enrContents")
                )
            )
        )
    )
)

server <- function(input, output, session) {
    
    output$downloadExample <- downloadHandler(
        filename = "switch_enrichments_unordered_moderate.csv",
        content = function(con) {
            write_csv(example_data, con)
        }
    )
    
    output$figureDownload <- downloadHandler(
        filename = paste0("enrichment_figure_", Sys.Date(), ".pdf"),
        content = function(con) {
            ggsave(
                plot = gg_plot(),
                filename = con,
                width = input$figwidth,
                height = input$figheight,
                device = "pdf"
            )
        }
    )
    
    enr_data_in <- reactive({
        switch(
            EXPR = input$enrfile$type,
            "text/plain" = read_tsv(input$enrfile$datapath),
            "text/tab-separated-values" = read_tsv(input$enrfile$datapath),
            "text/comma-separated-values" = read_csv(input$enrfile$datapath),
            "text/csv" = read_csv(input$enrfile$datapath)
        )
    })
    
    observe({
        req(input$enrfile)
        
        updateCheckboxGroupInput(
            session,
            "domains", 
            label = "Select domains to display",
            choices = sort(unique(enr_data_in()$domain)),
            selected = sort(unique(enr_data_in()$domain))
        )
    })
    
    output$intersection <- renderText({
        req(input$enrfile)
        if (is.null(input$termIntersect)) {
            return("")
        } else {
            return(
                paste(str_split(pull(filter(enr_data_in(), term.name == input$termIntersect), intersection), ",", simplify = TRUE), collapse = "\n")
            )
        }
    })

    enr_data <- reactive({
        if (is.null(input$numDisplay)) {
            x <- arrange(filter(enr_data_in(), domain %in% input$domains), p.value, desc(term.name))[
                1:nrow(enr_data_in()), 
            ]
        } else {
            x <- arrange(filter(enr_data_in(), domain %in% input$domains), p.value, desc(term.name))[
                input$numDisplay[[1]]:input$numDisplay[[2]], 
            ]
        }
        x[["p.value"]] <- x[[input$PValColLoc]]
        x[["term.name"]] <- x[[input$TermNameColLoc]]
        x[["intersection"]] <- x[[input$GnOverlapColLoc]]
        return(x)
    })
    
    gg_plot <- reactive({
        return(build_graph_gg(
            enr_data(),
            fontsize_labels = input$fontSizeLab,
            fontsize_axes = input$fontSizeAx,
            colorpalette = setNames(
                RColorBrewer::brewer.pal(length(unique(enr_data_in()$domain)), name = "Set2"),
                unique(enr_data_in()$domain)
            ),
            pvaluecolumn = input$PValColLoc,
            termnamecolumn = input$TermNameColLoc,
            geneoverlapcolumn = input$GnOverlapColLoc,
            pncolumn = input$NegFlagColLoc,
            pflag = input$PosFlag,
            nflag = input$NegFlag,
            vcut = -log10(input$PValueThresh)
        ))
    })
    
    output$numDisplay <- renderUI({
        req(input$enrfile)
        sliderInput(
            inputId = "numDisplay",
            label = "Terms to display",
            min = 1, max = nrow(enr_data_in()), value = c(1, 5), step = 1
        )
    })
   
    output$termSelect <- renderUI({
        req(input$enrfile)
        term_choices <- enr_data_in()$term.name
        selectInput(
            inputId = "termIntersect",
            label = "Term to show intersection",
            choices = term_choices 
        )
    })
    
    output$PValueColumn <- renderUI({
        req(input$enrfile)
        col_choices <- colnames(enr_data_in())
        selectInput(
            inputId = "PValColLoc",
            label = "PValue Column",
            choices = col_choices 
        )
    })
    
    output$TermNameColumn <- renderUI({
        req(input$enrfile)
        col_choices <- colnames(enr_data_in())
        selectInput(
            inputId = "TermNameColLoc",
            label = "Term Name Column",
            choices = col_choices 
        )
    })
    
    output$GeneOverlapColumn <- renderUI({
        req(input$enrfile)
        col_choices <- colnames(enr_data_in())
        selectInput(
            inputId = "GnOverlapColLoc",
            label = "Gene Overlap Column",
            choices = col_choices 
        )
    })
    
    output$NegativeFlagColumnUI <- renderUI({
        req(input$enrfile)
        col_choices <- colnames(enr_data_in())
        selectInput(
            inputId = "NegFlagColLoc",
            label = "+/- Flag Column",
            choices = col_choices
        )
    })
    
    output$PositiveFlagsInput <- renderUI({
        req(input$enrfile)
        textInput(
            inputId = "PosFlag",
            label = "+ Flag",
            value = ""
        )
    })
    
    output$NegativeFlagsInput <- renderUI({
        req(input$enrfile)
        textInput(
            inputId = "NegFlag",
            label = "- Flag",
            value = ""
        )
    })
    
    output$enrPlot <- renderPlot({
        req(input$enrfile)
        return(gg_plot())
    })
    
    output$enrContents <- DT::renderDataTable({
        req(input$enrfile)
        return(DT::datatable(enr_data()))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
