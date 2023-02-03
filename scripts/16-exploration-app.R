# load packages
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
if(!require(Seurat)){
  install.packages("Seurat")
  library(Seurat)
}
if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)
}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(patchwork)){
  install.packages("patchwork")
  library(patchwork)
}

# load data and scripts
load("results/0-objects/obj_annotated_diet.Rdata")
source("scripts/etc/colors.R")
deg <- read.csv("results/12-differential-gene-expression/dge_no_threshold.csv")
deg <- deg %>% mutate(significant = if_else(p_val_adj <= 0.01, "Yes", "No"))

# ui
ui <- fluidPage(
  verticalLayout(
    titlePanel("Bottermann et al. scRNAseq-DREADD-24h"),
    wellPanel(textInput("gene",
                        label = "Gene (e.g. Col1a1)"),
              splitLayout(cellWidths = c("50%", "50%"), 
                          plotOutput("Dim"),
                          plotOutput("Feature")),
              splitLayout(cellWidths = c("100%"), 
                          plotOutput("Violin_split")),
              splitLayout(cellWidths = c("100%"),
                          dataTableOutput("table"))),
    helpText("daniel.gorski@uni-duesseldorf.de")))

# server
server <- function(input, output) {
  output$Dim <- renderPlot({
    a <- DimPlot(obj,
                 reduction = "umap",
                 pt.size = .3,
                 raster = F,
                 label = F,
                 cols = colors,
                 group.by = "basic_annotation") +
      xlab("UMAP-1") +
      ylab("UMAP-2") +
      labs(color = "Identity") +
      theme(legend.text = element_text(size = 13),
            legend.title = element_text(size = 16),
            legend.justification = "top",
            legend.key.size = unit(3, "point"),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_text(size = 16),
            plot.title = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 4.25),
                                  nrow = 37))
    d <- LabelClusters(plot = a,
                       id = "basic_annotation",
                       repel = T,
                       force = 0.25,
                       box = T,
                       fill = alpha("white", 0.45),
                       size = 4,
                       label.r = unit(0.25, "lines"),
                       label.size = NA)
    print(d)
  })
  output$Feature <- renderPlot({
   FeaturePlot(obj, feature = input$gene) +
      theme(
        plot.title = element_text(face = "italic", size = 20)
      )
  })
  output$Violin_split <- renderPlot({
    v <- VlnPlot(obj,
                 features = input$gene,
                 split.by = "genotype",
                 split.plot = TRUE) +
      theme(
        plot.title = element_text(face = "italic", size = 20),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.position = "right"
      )
    print(v)
  })
  output$table <- renderDataTable(
    deg %>% filter(gene == input$gene) %>% select(cluster,
                                                  gene,
                                                  p_val_adj,
                                                  significant,
                                                  cluster,
                                                  regulation,
                                                  avg_log2FC,
                                                  pct.1,
                                                  pct.2),
    options = list(paging = FALSE, searching = FALSE)
  )
}

# app----
runApp(list(ui = ui, server = server), launch.browser = TRUE)
