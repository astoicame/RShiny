library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker)
library(readr)
library(tidyverse)
library(DESeq2)
library(dplyr)
library(plotly)
library(heatmaply)
library(gplots)
library(matrixStats)
library(scales)
library(stats)
library(sparseMatrixStats)
library(DelayedMatrixStats)
library(fgsea)
library(biomaRt)
library(DT)
library(pheatmap)
library(ggbeeswarm)

# setup options and vectors for later use
options(shiny.maxRequestSize=30*1024^2)
deseq_choices <-
  c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
viz_choices <- c("Post_Mortem_Interval", "Age_of_Death", "RNA_Integrity_Number", "mRNAseq_Reads")


# UI setup
ui <- fluidPage(
  theme = bs_theme(bootswatch = "minty"),
  titlePanel("BF591 Final Project - Andreea Stoica"),
      tabsetPanel(
        # First tab, input for sample summary
        tabPanel("Samples",
          sidebarLayout(
            sidebarPanel(
              fileInput("file", "Upload CSV file"),
              radioButtons("x_variable", "Select X Variable", choices = viz_choices)),
          mainPanel(tabsetPanel(
            tabPanel("Summary", dataTableOutput("summary_table")),
            tabPanel("Table", DTOutput("data_table")),
            tabPanel("Plots", plotOutput("plots"))
                )
              )
            )
          ),
        
        # Second tab, input for gene counts analysis
        tabPanel(
          "Counts",
          sidebarLayout(
            sidebarPanel(
              fileInput("counts_file", "Upload normalized counts matrix (CSV)"),
              sliderInput("variance", "Minimum Variance Percentile", min = 1, max = 100, value = 90),
              sliderInput("non_zero_samples", "Minimum Non-zero Samples", min = 1, max = 100, value = 5)
            ),
  
            mainPanel(tabsetPanel(
                tabPanel("Filter Results", 
                         tableOutput("summary")),
                tabPanel("Diagnostics",
                         plotOutput("counts_variance_median", width = "1000px"),
                         plotOutput('counts_nonzero_vs_median', width = "1000px")),
                tabPanel("Heatmap", 
                         plotOutput("heatmap")),
                tabPanel("PCA", 
                         sliderInput("n_components_slider", "Select (N) components ", min = 1, max = 10, value = 8),
                         plotOutput("counts_PCA", width = "600"))
                )
              )
            )
          ),
         
        # Third tab, input for DESeq analysis display 
        tabPanel(
          "DESeq",
          sidebarLayout(
            sidebarPanel(
              fileInput(
                "deseq_file",
                label = "Load differential expression results",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv"),
                placeholder = "deseq_res.csv"
              ),
              HTML(
                paste(
                  "A volcano plot can be generated with",
                  "<b>\"log<sub>2</sub> fold-change\"</b> on the x-axis and",
                  "<b>\"p-adjusted\"</b> on the y-axis.<br>"
                )
              ),
              br(),
              radioButtons(
                inputId = "x_axis",
                label = "Choose the column for the x-axis",
                choices = deseq_choices,
                selected = "log2FoldChange"
              ),
              radioButtons(
                inputId = "y_axis",
                label = "Choose the column for the y-axis",
                choices = deseq_choices,
                selected = "padj"
              ),
              colourInput(
                inputId = "base",
                label = "Base point color",
                value = "#22577A",
                closeOnClick = T
              ),
              colourInput(
                inputId = "highlight",
                label = "Highlight point color",
                value = "#FFCF56",
                closeOnClick = T
              ),
              sliderInput(
                "slider",
                "Select the magnitude of the p adjusted coloring:",
                min = -50,
                max = 0,
                value = -5,
                step = 1
              ),
            ),
            mainPanel(tabsetPanel(
              tabPanel("Plot", {
                plotOutput("volcano")
              }),
              tabPanel("Table",
                       tableOutput("table"))
              )
            )
          )
        ),
        
        #Fourth tab, input for GSEA analysis results
        tabPanel("GSEA",
                 sidebarLayout(
                   sidebarPanel(
                     fileInput("fgsea_results_file", "Upload fgsea results file", accept = c(".csv", ".tsv"))
                   ),
                   mainPanel(
                     tabsetPanel(
                       # Tab 1 - Top Pathways Plot
                       tabPanel("Top Pathways Plot",
                                sidebarLayout(
                                  sidebarPanel(
                                    sliderInput("number_paths", "Number of Pathways", min = 0, max =50, value = 5)
                                  ),
                                  mainPanel(
                                    plotOutput("top_pathways_plot") #  plot for the top pathways!
                            )
                          )
                       ),
                       
                       tabPanel("Table Results",
                                sidebarLayout(
                                  sidebarPanel(
                                    sliderInput("top_paths", "Number of Top Pathways to Display:",
                                                min = 1, max = 50, value = 10),
                                    radioButtons("nes_direction", "NES Direction:",
                                                 choices = c("All", "Positive", "Negative"), selected = "All"),
                                    downloadButton("download_table", "Download Table")
                                  ),
                                  mainPanel(
                                    DTOutput("fgsea_table")
                                  )
                                )
                              ),
                      
                       tabPanel("Scatter Plot",
                                sidebarLayout(
                                  sidebarPanel(
                                    sliderInput("pvalue_threshold", "Adjusted p-value Threshold:",
                                                min = 1.281914e-11, max = 1, value = 0.05),
                                    width = 5
                                  ),
                                  mainPanel(
                                    plotOutput("scatter_plot", height = "400px"),
                                    width = 9
                                  )
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )



server <- function(input, output, session) {
  #read in csv with metadata
  sample_data <- reactive({
    req(input$file)
    inFile <- input$file
    read.csv(inFile$datapath)
  })
  
  
  # summary table 
  output$summary_table <- renderDT({
    req(sample_data())
    col_types <- sapply(sample_data(), class)
    col_means <- sapply(sample_data(), function(x) ifelse(is.numeric(x), mean(x, na.rm = TRUE), NA))
    
    summary_df <- data.frame(
      Data_Type = col_types,
      Mean_Value = col_means
    )
    
    datatable(summary_df, options = list(pageLength = 11))
  })
  
  
  # show data table
  output$data_table <- renderDT({
    req(sample_data())
    datatable(sample_data(), options = list(pageLength = 10))
  })
  
  
  # create distribution plots
  output$plots <- renderPlot({
    df <- sample_data()
    selected_var <- input$x_variable
    p <- isolate(hist(as.numeric(df[[selected_var]]), main = paste("Distribution of", selected_var), xlab = selected_var))
    return(p)
  })
  
  
  #DESEQ read in data and create volcano plot as well as table to display results
  deseq_data <- reactive({
    df <- read.csv(input$deseq_file$datapath)
    colnames(df)[1] <- "gene"
    return(df)
  })
  
  
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
      p <- ggplot(dataf, aes(x = !!sym(x_name),
                             y = -log10(!!sym(y_name)))) +
        geom_point(aes(color = !!sym(y_name) < 1 * 10 ^ (as.numeric(slider)))) +
        theme_bw() +
        scale_color_manual(values = c(color1, color2)) +
        theme(legend.position = "bottom") +
        labs(color = paste0(y_name, " < 1 × 10^", slider))
      return(p)
    }
  draw_table <- function(dataf, slider) {
    df_out <- dataf[which(dataf$padj < 1 * 10 ^ (as.numeric(slider))),]
    df_out$pvalue <- formatC(df_out$pvalue, digits = -2)
    df_out$padj <- formatC(df_out$padj, digits = -2)
    return(df_out)
  }
  
  output$volcano <- renderPlot({
    req(input$deseq_file)
    df <- deseq_data()
    p <-volcano_plot(df,
                     input$x_axis,
                     input$y_axis,
                     input$slider,
                     input$base,
                     input$highlight)
    return(p)}, height = 700)
  
  # Same here, just return the table as you want to see it in the web page
  output$table <- renderTable({
    req(input$deseq_file)
    table <- deseq_data()
    colnames(table)[1] <- "gene"
    return(draw_table(dataf = table, slider = input$slider))
  }, striped = T)
  
  
  #read in counts data
  counts_data <- reactive({
    req(input$counts_file)
    inFile <- input$counts_file
    if (is.null(inFile))
      return(NULL)
    read_csv(inFile$datapath)
  })
  
  
  #function to filter genes based on parameters set by sliders in app
  filter_genes <- reactive({
    
    counts_matrix <- counts_data()[-1]
    
    # variance threshold
    variance_threshold <- quantile(apply(counts_matrix, 1, var), input$variance / 100)
    genes_passing_variance <- rownames(counts_matrix)[apply(counts_matrix,1, var )  >= variance_threshold]
    
    # filter genes based on non-zero samples
    nonzero_samples_threshold <- input$non_zero_samples / 100 * ncol(counts_matrix)
    genes_passing_nonzero_samples <- rownames(counts_matrix)[apply(counts_matrix > 0,1, sum) >= nonzero_samples_threshold]
    
    # selecting the genes by intersecting 
    selected_genes <- intersect(genes_passing_variance, genes_passing_nonzero_samples)
    
    # creating a filtered matrix 
    filtered_counts_matrix <- counts_matrix[selected_genes, ]
    
    
    # creating the metrics
    total_genes <- nrow(counts_matrix)
    passing_genes <- length(selected_genes)
    non_passing_genes <- total_genes - passing_genes
    percent_passing <- (passing_genes / total_genes) * 100
    percent_non_passing <- (non_passing_genes / total_genes) * 100
    
    # creating a summary tibble
    summary_data <- tibble::tibble(
      Metric = c("Number of Genes", "Genes Passing Filter", "% of Genes Passing Filter", "Genes Not Passing Filter", "% of Genes Not Passing Filter"),
      Value = c(total_genes, passing_genes, percent_passing, non_passing_genes, percent_non_passing)
    )
    
    print(summary_data)
  })
   
  
  #PCA beeswarm plot of # of PCs as input by slider
  filter_genes_pca_beeswarm <- reactive({
    counts_matrix <- counts_data()[-1]
    
    # variance threshold
    variance_threshold <- quantile(apply(counts_matrix, 1, var), input$variance/ 100)
    genes_passing_variance <- rownames(counts_matrix)[apply(counts_matrix,1, var )  >= variance_threshold]
    
    # filter genes based on non-zero samples
    nonzero_samples_threshold <- input$non_zero_samples / 100 * ncol(counts_matrix)
    genes_passing_nonzero_samples <- rownames(counts_matrix)[apply(counts_matrix > 0,1, sum) >= nonzero_samples_threshold]
    
    # selecting the genes by intersecting parameters
    selected_genes <- intersect(genes_passing_variance, genes_passing_nonzero_samples)
    
    # creating a filtered matrix 
    filtered_count_matrix <- counts_matrix[selected_genes, ]
    
    #creating the plot
    pca_results_filter <- prcomp(t(filtered_count_matrix))
    pc_variances <- pca_results_filter$sdev^2 / sum(pca_results_filter$sdev^2) * 100
    pca_data_filter <- as.data.frame(pca_results_filter$x)
    
    pca_data_subset_filter <- pca_data_filter[, 1:input$n_components_slider]
    
    plot_data_filter <- reshape2::melt(pca_data_subset_filter)
    
    print(head(plot_data_filter))
    pca_plot <- ggplot(plot_data_filter) +
      geom_beeswarm(aes(x = variable, y = value), color = 'red') +
      xlab("Principal Components") +
      ylab("Values") +
      ggtitle("PC Beeswarm Plot") 
      
    return(pca_plot)
  })
  
  
  #output for summary table of filtered genes and the plots based on slider parameters
  output$summary <- renderTable({
    filter_genes()
  })
  
  #output variance vs median plot
  output$counts_variance_median <- renderPlot({
    # function to plot the variance vs median counts
    counts_matrix <- counts_data()[-1]
    
    # finding median and variance of matrix
    median_count <- apply(counts_matrix, 1, median)
    variance_count <- apply(counts_matrix, 1, var)
    
    plotting_data <- data.frame(MedianCount = median_count, Variance = variance_count)
    
    # ranking by medians
    plotting_data$rank <- rank(median_count)
    
    variance_quantiles <- quantile(variance_count, probs = c(0, input$variance / 100, 1), na.rm = TRUE)
    plotting_data$color <- cut(variance_count, breaks = variance_quantiles, labels = c("Below Threshold", "Above Threshold"), include.lowest = TRUE)
    
    
    ggplot(plotting_data, aes(x = rank, y = Variance, color = color)) + 
      geom_point() +
      labs(x = "Rank(median)", y = "Variance", title = "Median vs Variance") +
      scale_y_log10() +
      scale_color_manual(values = c("Above Threshold" = "#FFCF56", "Below Threshold" = "#22577A"))  # Customize colors as needed
    
  })
  
  
  output$counts_nonzero_vs_median <- renderPlot({
    
    counts_matrix <- counts_data()[-1]
    
    median_count <- apply(counts_matrix, 1, median)
    # count of zeroes
    zero_count <- apply(counts_matrix == 0, 1, sum)
    
    plotting_data <- data.frame(MedianCount = median_count, ZeroCount = zero_count)
    
    # ranking by medians
    plotting_data$rank <- rank(median_count)
    
    # highlighting based on the zero slider
    plotting_data$color <- ifelse(zero_count >= input$non_zero_samples / 100 * ncol(counts_matrix), "Above Threshold", "Below Threshold")
    
    ggplot(plotting_data, aes(x = rank, y = ZeroCount, color = color)) + 
      geom_point() +
      labs(x = "Rank(median)", y = "Zero Count", title = "Median vs Zero Count") +
      #scale_y_reverse() +
      scale_color_manual(values = c("Below Threshold" = "#FFCF56", "Above Threshold" = "#22577A"))  # Customize colors as needed
  })
  
  
  # heatmap based on gene counts
  output$heatmap <- renderPlot({
    counts_matrix <- counts_data()[-1]
    
    
    # variance threshold
    variance_threshold <- quantile(apply(counts_matrix, 1, var), input$variance / 100)
    genes_passing_variance <- rownames(counts_matrix)[apply(counts_matrix,1, var )  >= variance_threshold]
    
    # filter genes based on non-zero samples
    nonzero_samples_threshold <- input$non_zero_samples / 100 * ncol(counts_matrix)
    genes_passing_nonzero_samples <- rownames(counts_matrix)[apply(counts_matrix > 0,1, sum) >= nonzero_samples_threshold]
    
    # selecting the genes by intersecting parameters
    selected_genes <- intersect(genes_passing_variance, genes_passing_nonzero_samples)
    
    # creating a filtered matrix 
    heat_matrix <- counts_matrix[selected_genes, ]
    
    
    heat_matrix <- as.matrix(heat_matrix)
    heat_matrix <- log2(heat_matrix + 1)
    
    # creating a heatmap based on the filtered matrix, this also needs to be reactive, similar to the summary count
    plot <- pheatmap(heat_matrix, 
                     clustering_distance_rows = "euclidean", 
                     clustering_distance_cols = "euclidean", 
                     scale = "row", 
                     main = "Clustered Heatmap",
                     width =15,
                     height = 8)
  })
  
  
  #display beeswarm plot
  output$counts_PCA <- renderPlot({
    
    filter_genes_pca_beeswarm()
    
  })
  
  
  #input fgsea file for GSEA analysis display
  fgsea_data <- reactive({
    req(input$fgsea_results_file)
    inFile <- input$fgsea_results_file
    read.csv(inFile$datapath)
  })
  
  
  #get top affected pathways
  top_pathways <- reactive({ 
    # filtering pathways based on log-transformed padj threshold
    threshold <- input$pvalue_threshold_slider
      top_pos <- fgsea_data() %>% slice_max(NES, n=input$number_paths) %>% pull(pathway)
      top_neg <- fgsea_data() %>% slice_min(NES, n=input$number_paths) %>% pull(pathway)
      
      subset <- fgsea_data() %>% 
        filter(pathway %in% c(top_pos, top_neg)) %>%
        mutate(pathway = factor(pathway)) %>%
        mutate(plot_name = str_replace_all(pathway, '_', ' '))
      
      plot <- subset %>% 
        mutate(plot_name = forcats::fct_reorder(factor(plot_name), NES)) %>%
        ggplot() +
        geom_bar(aes(x=plot_name, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
        scale_fill_manual(values = c('TRUE' = 'red', 'FALSE' = 'blue')) + 
        theme_minimal(base_size = 8) +
        ggtitle('fgsea results for Hallmark MSigDB gene sets') +
        ylab('Normalized Enrichment Score (NES)') +
        xlab('') +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 80)) +
        coord_flip()
      return(plot)
  })
  
  
  filter_results <- reactive({
    filtered_results <- switch(input$nes_direction,
                               "All" = fgsea_data(),
                               "Positive" = fgsea_data()[fgsea_data()$NES > 0, ],
                               "Negative" = fgsea_data()[fgsea_data()$NES < 0, ])
    return(filtered_results)
  })
  
  
  #output top pathways bar graph
  output$top_pathways_plot <- renderPlot({
    # stacked  bar chart of the top pathways identified
    top_pathways()
  })
  
  
  #output table displaying top pathways
  output$fgsea_table <- renderDT({
    # filtering based on NES direction and select top pathways by adjusted p-value
    filtered_data <- head(filter_results()[order(filter_results()$padj), ], input$top_paths)
    datatable(filtered_data, options = list(order = list(list(3, 'asc'))))
  })
  
  
  # creating download button
  output$download_table <- downloadHandler(
    filename = function() {
      paste("fgsea_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # write filtered data table to a CSV file
      write.csv(filter_results(), file, row.names = FALSE)
    }
  )
  
  #inputting threshold value for p adjusted and filtering by it
  filtered_fgsea <- reactive({
    fgsea_data() %>%
      filter(padj <= input$pvalue_threshold)
  })
  
  
  #creating scatter plot
  output$scatter_plot <- renderPlot({
    ggplot(fgsea_data, aes(x = NES, y = -log10(padj), color = padj <= input$pvalue_threshold)) +
      geom_point(aes(fill = padj <= input$pvalue_threshold), shape = 21, size = 3) +
      scale_fill_manual(values = c("grey", "#FFCF56"), guide = FALSE) +
      scale_color_manual(values = c("grey", "#FFCF56"), guide = FALSE) +
      labs(x = "NES", y = "-log10(Adjusted p-value)") +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
}



shinyApp(ui = ui, server = server)




