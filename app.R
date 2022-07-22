library(shiny)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(dplyr)
library(DT)
library(Biostrings)
library(msa)
library(ggseqlogo)

synovium <- read.csv("synovium_TRB_curated.csv", header = T, stringsAsFactors = F)
epitopes <- read.csv("annotated_epitopes_paratopes.csv", header = T, stringsAsFactors = F)



ui <- fluidPage(
  titlePanel("Rheumatoid Arthritis Synovial TCR Repertoire Database"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("cdr3AAlen_Input", "TRB CDR3 Length (AA)", 1, 30, c(10, 25)),
      checkboxGroupInput("dataset_Input", "Dataset Origin",
                   choices = c("GSE89408", "PEAC", "McPAS", "VDJdb"),
                   selected = c("GSE89408", "PEAC", "McPAS", "VDJdb")),
      checkboxGroupInput("MHC_Input", "MHC Class",
                   choices = c("MHCI", "MHCII", "Unknown"),
                   selected = c("MHCI", "Unknown")),
      checkboxGroupInput("Pathogen_Input", "Conditions to Consider (Please do not select all)", 
                         choices = c("RA", "CMV", "SARS-CoV-2", "EBV", "HIV", 
                                     "HCV", "Influenza", "DENV", "YFV"),
                         selected =c("RA", "SARS-CoV-2")
                         )
      
          ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Size Distribution", plotOutput("SizeDistPlot")),
                  tabPanel("Repertoire Table", dataTableOutput("SynvTable")),
                  tabPanel("Annotated Epitopes", dataTableOutput("AnnoTable")),
                  tabPanel("Motif Logos", plotOutput("MotifLogoPlot", height = '1500px'))
      
      )
)))

              

server <- function(input, output) {
  output$SizeDistPlot <- renderPlot({
    filtered_synovium <-
      synovium %>%
      filter(TRB_Length >= input$cdr3AAlen_Input[1],
             TRB_Length <= input$cdr3AAlen_Input[2],
             DatasetOrigin == input$dataset_Input
      )
    p1 <- 
      ggplot(filtered_synovium, aes(ReadCount)) +
      geom_histogram(color= 'midnightblue', fill = 'midnightblue')+theme_classic()+
      scale_x_continuous(trans = 'log2')+scale_y_continuous(trans = 'log2')+
      labs(x = 'Read Counts of Clonotype', y = 'Number of Clonotypes')+
      theme(text = element_text(size = 14))
    
    p2 <- 
      ggplot(filtered_synovium, aes(ClonalFreq)) +
      geom_histogram(color= 'midnightblue', fill = 'midnightblue')+theme_classic()+
      scale_y_continuous(trans = 'log2')+
      labs(x = 'In-Sample Frequency', y = 'Number of Clonotypes')+
      theme(text = element_text(size = 14))
    
    p1+p2
    
  })
  
  output$SynvTable <- renderDataTable({
    filtered_synovium <-
      synovium %>%
      filter(TRB_Length >= input$cdr3AAlen_Input[1],
             TRB_Length <= input$cdr3AAlen_Input[2],
             DatasetOrigin %in% input$dataset_Input
      )
    filtered_synovium
  }, rownames = FALSE) 
  
  output$AnnoTable <- renderDataTable({
  filtered_epitopes <-
    epitopes %>%
    filter(TRB_Length >= input$cdr3AAlen_Input[1],
           TRB_Length <= input$cdr3AAlen_Input[2],
           Condition %in% input$Pathogen_Input,
           MHC_Class %in% input$MHC_Input, 
           Database %in% input$dataset_Input)
  
  list_pathogens <- as.list(input$Pathogen_Input)
  partial_epitopes <- lapply(list_pathogens, function(x)  
    filtered_epitopes[filtered_epitopes$Condition %in% x, ])
  
  partial_epitopes <- lapply(partial_epitopes, "[", 2)
  partial_epitopes <- lapply(partial_epitopes, function(x) unlist(as.list(x)))
    joint_index <- Reduce(intersect, as.list(partial_epitopes))
    
    annotated_epitopes <-filtered_epitopes[filtered_epitopes$index %in% joint_index, ]
 
 
     annotated_epitopes
  }, rownames = FALSE) 
  
  output$MotifLogoPlot <- renderPlot({
    filtered_epitopes <-
      epitopes %>%
      filter(TRB_Length >= input$cdr3AAlen_Input[1],
             TRB_Length <= input$cdr3AAlen_Input[2],
             Condition %in% input$Pathogen_Input,
             MHC_Class %in% input$MHC_Input, 
		         Database %in% input$dataset_Input)
    
    list_pathogens <- as.list(input$Pathogen_Input)
    partial_epitopes <- lapply(list_pathogens, function(x)  
      filtered_epitopes[filtered_epitopes$Condition %in% x, ])
    
    partial_epitopes <- lapply(partial_epitopes, "[", 2)
    partial_epitopes <- lapply(partial_epitopes, function(x) unlist(as.list(x)))
    joint_index <- Reduce(intersect, as.list(partial_epitopes))
    
    annotated_epitopes <-filtered_epitopes[filtered_epitopes$index %in% joint_index, ]
    
    epitope_motifs <- levels(as.factor(annotated_epitopes$index))
    list_motifs <- lapply(epitope_motifs, function(x) 
      annotated_epitopes[annotated_epitopes$index %in%  x, ])
    list_motifs <- lapply(list_motifs, "[", 8)
    list_motifs <- lapply(list_motifs, function(x) unlist(as.list(x)))
    
    len_list_motifs <- lapply(list_motifs, function(x) length(x))
    list_motifs <- list_motifs[len_list_motifs > 4]
    len_list_motifs <- len_list_motifs[len_list_motifs > 4]
    list_motifs <- list_motifs[len_list_motifs < 11]
    
    list_motifs_align <- lapply(list_motifs, function(x) AAStringSet(x))
    list_motifs_align <- lapply(list_motifs_align, function(x) msa(x,method="ClustalOmega", type = 'protein'))
    out_aligns <- lapply(list_motifs_align, function(x) x@unmasked)
    out_aligns <- lapply(out_aligns, function(x) as.character(x))
      
  col_logos <- as.list(lapply(out_aligns, function(x) ggplot() + 
                        geom_logo(x, method = 'prob')+ theme_logo()+
                        theme(text = element_text(size = 14))+
                                guides(color=guide_legend(ncol =2))
                      )) 
  do.call("grid.arrange", c(col_logos, ncol = 3))
  
  })
    
     }


shinyApp(ui = ui, server = server)