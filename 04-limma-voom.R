# RNA-seq pipeline
# Ben Laufer

# Modifies and expands on these references:
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

# Load packages -----------------------------------------------------------

setwd("~/Box/PEBBLES/RNA")

packages <- c("edgeR", "tidyverse", "annotables", "RColorBrewer", "org.Mm.eg.db", "EnhancedVolcano",
              "enrichR", "openxlsx", "gt", "glue", "Glimma", "sva", "DMRichR")
enrichR:::.onAttach() # Needed or else "EnrichR website not responding"
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

#BiocManager::install("ben-laufer/DMRichR")
#BiocManager::install("stephenturner/annotables")

# To test and develop, assign the variable tissue and then just run the main sections

sink("RNA-seq_log.txt", type = "output", append = FALSE, split = TRUE)

purrr::pwalk(tidyr::crossing(tissue = c("placenta", "brain"),
                             sex = c("male", "female")),
             function(tissue, sex){
               
               dir.create(glue::glue("{tissue}_{sex}"))
  
  # Count Matrix ------------------------------------------------------------
  
  #name <- gsub( "(?:[^_]+_){4}([^_ ]+)*$","", files)
  
  # STAR quantMode geneCounts output:
  #column 1: gene ID
  #column 2: counts for unstranded RNA-seq
  #column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
  #column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
  
  # KAPA mRNA HyperPrep Kit reads are reverse stranded, so select column 4
  # Confirm by looking at the N_noFeature line for the 3rd and 4th column and pick the column with the lowest count.
  
  sampleNames <- list.files(path = glue::glue(getwd(), "/GeneCounts"),
                            pattern = "*.ReadsPerGene.out.tab") %>%
    stringr::str_split_fixed("_", n = 3) %>%
    tibble::as_tibble() %>%
    tidyr::unite(Name, c(V1:V2), sep = "-") %>%
    dplyr::select(Name) %>% 
    purrr::flatten_chr()
  
  # Could alternatively use edgeR::readDGE() but that calls to the slower read.delim()
  geneSymbols <- list.files(path = glue::glue(getwd(), "/GeneCounts"),
                            pattern = "*.ReadsPerGene.out.tab", full.names = TRUE)[1] %>% 
    data.table::fread(select = 1) %>%
    purrr::flatten_chr()
  
  countMatrix <- list.files(path = glue::glue(getwd(), "/GeneCounts"),
                            pattern = "*.ReadsPerGene.out.tab", full.names = TRUE) %>%
    purrr::map_dfc(data.table::fread, select = 4, data.table = FALSE) %>%
    magrittr::set_colnames(sampleNames) %>% 
    magrittr::set_rownames(geneSymbols)
  
  # Remove meta info
  countMatrix <- countMatrix[-c(1:4),]
  
  # Design Matrix -----------------------------------------------------------
  
  designMatrix <- readxl::read_xlsx("sample_info.xlsx") %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::mutate(Name = as.character(Name))
  
  # # Recode sex
  # designMatrix$Sex <- as.character(designMatrix$Sex)
  # designMatrix$Sex[designMatrix$Sex  == "F"] <- "0"
  # designMatrix$Sex[designMatrix$Sex  == "M"] <- "1"
  # designMatrix$Sex <- as.factor(designMatrix$Sex)
  
  samples.idx <- pmatch(designMatrix$Name, colnames(countMatrix))
  designMatrix <- designMatrix[order(samples.idx),]
  
  # Preprocessing -----------------------------------------------------------
  
  print(glue::glue("Preprocessing {sex} {tissue} samples"))
 
   # Select sample subset
  designMatrix <- designMatrix %>%
    dplyr::filter(Tissue == tissue & Sex == sex)
  
  countMatrix <- countMatrix %>%
    dplyr::select(contains(designMatrix$Name))
  
  # Create DGE list and calculate normalization factors
  countMatrix <- countMatrix %>%
    DGEList() %>%
    calcNormFactors()
  
  # Reorder design matrix 
  samples.idx <- pmatch(designMatrix$Name, rownames(countMatrix$samples))
  designMatrix <- designMatrix[order(samples.idx),]
  stopifnot(rownames(countMatrix$samples) == designMatrix$Name)
  
  # Add sample info from design matrix to DGE list
  countMatrix$samples$group <- designMatrix$Treatment
  countMatrix$samples$Sex <- designMatrix$Sex
  countMatrix$samples$Litter <- designMatrix$Litter
  countMatrix$samples$Tissue <- designMatrix$Tissue
  
  # Raw density of log-CPM values
  
  L <- mean(countMatrix$samples$lib.size) * 1e-6
  M <- median(countMatrix$samples$lib.size) * 1e-6
  
  logCPM <- cpm(countMatrix, log = TRUE)
  logCPM.cutoff <- log2(10/M + 2/L)
  nsamples <- ncol(countMatrix)
  col <- brewer.pal(nsamples, "Paired")
  
  pdf(glue::glue("{tissue}_{sex}/{tissue}_{sex}_density_plot.pdf"), height = 8.5, width = 11)
  par(mfrow = c(1,2))
  
  plot(density(logCPM[,1]), col = col[1], lwd = 2, las = 2, main = "", xlab = "")
  title(main = "A. Raw data", xlab = "Log-cpm")
  abline(v = logCPM.cutoff, lty = 3)
  for (i in 2:nsamples){
    den <- density(logCPM[,i])
    lines(den$x, den$y, col = col[i], lwd = 2)
  }
  legend("topright", designMatrix$Name, text.col = col, bty = "n", cex = 0.5)
  
  # Filter genes with low expression
  
  rawCount <- dim(countMatrix)
  
  keep.exprs <- filterByExpr(countMatrix,
                             group = countMatrix$samples$group,
                             lib.size = countMatrix$samples$lib.size)
  
  countMatrix <- countMatrix[keep.exprs,, keep.lib.sizes = FALSE] %>%
    calcNormFactors() 
  
  filterCount <- dim(countMatrix)
  
  print(glue::glue("{100 - round((filterCount[1]/rawCount[1])*100)}% of genes were filtered from {rawCount[2]} samples, \\
             where there were {rawCount[1]} genes before filtering and {filterCount[1]} genes after filtering for {tissue}"))
  
  # Filtered density plot of log-CPM values 
  logCPM <- cpm(countMatrix, log = TRUE)
  plot(density(logCPM[,1]), col = col[1], lwd = 2, las =2 , main = "", xlab = "")
  title(main = "B. Filtered data", xlab = "Log-cpm")
  abline(v = logCPM.cutoff, lty = 3)
  for (i in 2:nsamples){
    den <- density(logCPM[,i])
    lines(den$x, den$y, col = col[i], lwd = 2)
  }
  legend("topright", designMatrix$Name, text.col = col, bty = "n", cex = 0.5)
  dev.off()
  
  # Interactive MDS plot
  Glimma::glMDSPlot(countMatrix,
                    groups = designMatrix,
                    path = getwd(),
                    folder = "interactiveMDS",
                    html = glue::glue("{tissue}_{sex}_MDS-Plot"),
                    launch = FALSE)

  # Surrogate variables analysis --------------------------------------------
  
  # # Create model matrices, with null model for svaseq, and don't force a zero intercept
  # mm <- model.matrix(~Treatment + Litter,
  #                    data = designMatrix)
  # 
  # mm0 <- model.matrix(~1 + Litter,
  #                     data = designMatrix)
  # 
  # # svaseq requires normalized data that isn't log transformed
  # cpm <- cpm(countMatrix, log = FALSE)
  # 
  # # Calculate number of surrogate variables
  # nSv <- num.sv(cpm,
  #               mm,
  #               method = "leek")
  # 
  # # Estimate surrogate variables
  # svObj <- svaseq(cpm,
  #                 mm,
  #                 mm0,
  #                 n.sv = nSv)
  # 
  # # Update model to include surrogate variables
  # mm <- model.matrix(~Treatment + svObj$sv,
  #                    data = designMatrix)
  
  # Voom transformation and calculation of variance weights -----------------
  
  print(glue::glue("Normalizing {sex} {tissue} samples"))
  
  # Design
  mm <- model.matrix(~Treatment,
                     data = designMatrix)
  
  # Voom
  pdf(glue::glue("{tissue}_{sex}/{tissue}_{sex}_voom_mean-variance_trend.pdf"), height = 8.5, width = 11)
  voomLogCPM <- voom(countMatrix,
                     mm,
                     plot = TRUE)
  dev.off()
  
  # Make litter a random effect, since limma warns "coefficients not estimable" for some litters
  # Ref: https://support.bioconductor.org/p/11956/
  # Obstacle: Cannot do this properly with surrogtate variables, since there's an error when including litter in null model
  # Duplicate correlations alternative for other scenarios:
  # https://support.bioconductor.org/p/68916/
  # https://support.bioconductor.org/p/110987/
  
  correlations <- duplicateCorrelation(voomLogCPM,
                                       mm,
                                       block = designMatrix$Litter)

  # Extract intraclass correlation within litters
  correlations <- correlations$consensus.correlation
  
  # Boxplots of logCPM values before and after normalization
  pdf(glue::glue("{tissue}_{sex}/{tissue}_{sex}_normalization_boxplots.pdf"), height = 8.5, width = 11)
  par(mfrow=c(1,2))
  
  boxplot(logCPM, las = 2, col = col, main = "")
  title(main = "A. Unnormalised data", ylab = "Log-cpm")
  
  boxplot(voomLogCPM$E, las = 2, col = col, main = "")
  title(main = "B. Normalised data", ylab = "Log-cpm")
  
  dev.off()
  
  # Fitting linear models in limma ------------------------------------------
  
  print(glue::glue("Testing {sex} {tissue} samples for differential expression"))
  
  # Wieght standard errors of log fold changes by within litter correlation 
  fit <- lmFit(voomLogCPM,
               mm,
               correlation = correlations,
               block = designMatrix$Litter)
  
  head(coef(fit))
  
  # Save normalized expression values for WGCNA
  voomLogCPM$E %>%
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "Gene") %>%
    openxlsx::write.xlsx(glue::glue("{tissue}_{sex}/{tissue}_{sex}_voomLogCPMforWGCNA.xlsx"))
  
  # Create DEG tibble -------------------------------------------------------
  
  print(glue::glue("Creating DEG list of {sex} {tissue} samples"))
  
  DEGs <- fit %>%
    contrasts.fit(coef = 2) %>% # Change for different models
    eBayes() %>%
    topTable(sort.by = "P", n = Inf)
  
  DEGs <- DEGs %>%
    rownames_to_column() %>% 
    tibble::as_tibble() %>%
    dplyr::rename(symbol = rowname) %>% 
    dplyr::left_join(annotables::grcm38, by = "symbol") %>% 
    dplyr::select(symbol, logFC, P.Value, adj.P.Val, ensgene, description) %T>%
    openxlsx::write.xlsx(file = glue::glue("{tissue}_{sex}/{tissue}_{sex}_DEGs.xlsx"))
  
  # Volcano Plot ------------------------------------------------------------
  
  volcano <- DEGs %>% 
    EnhancedVolcano::EnhancedVolcano(title = "",
                                     labSize = 5,
                                     lab = .$symbol,
                                     x = 'logFC',
                                     y = 'P.Value', # P.Value 'adj.P.Val'
                                     col = c("grey30", "royalblue", "royalblue", "red2"),
                                     pCutoff = 0.05,
                                     FCcutoff = 0.0) +
    ggplot2::coord_cartesian(xlim = c(-3, 3),
                             ylim = c(0, 4))
  
  ggplot2::ggsave(glue::glue("{tissue}_{sex}/{tissue}_{sex}_volcano.pdf"),
                  plot = volcano,
                  device = NULL,
                  width = 11,
                  height = 8.5)
  
  # HTML report -------------------------------------------------------------
  
  print(glue::glue("Saving html report of {sex} {tissue} samples"))
  
  DEGs <- DEGs %>%
    dplyr::filter(P.Value < 0.05) %T>%
    openxlsx::write.xlsx(file = glue::glue("{tissue}_{sex}/{tissue}_{sex}_filtered_DEGs.xlsx"))
  
  DEGs %>%
    dplyr::rename(Gene = symbol,
                  "p-value" = P.Value,
                  "adjusted  p-value" = adj.P.Val,
                  Description = description) %>% 
    dplyr::select(-ensgene) %>%
    dplyr::mutate(Description = purrr::map_chr(strsplit(DEGs$description, split='[', fixed=TRUE),function(x) (x[1]))) %>% 
    gt() %>%
    tab_header(
      title = glue::glue("{nrow(DEGs)} Differentially Expressed Genes"),
      subtitle = glue::glue("{round(sum(DEGs$logFC > 0) / nrow(DEGs), digits = 2)*100}% up-regulated, \\
                              {round(sum(DEGs$logFC < 0) / nrow(DEGs), digits = 2)*100}% down-regulated")) %>% 
    fmt_number(
      columns = vars("logFC"),
      decimals = 2) %>% 
    fmt_scientific(
      columns = vars("p-value", "adjusted  p-value"),
      decimals = 2) %>%
    as_raw_html(inline_css = TRUE) %>%
    write(glue::glue("{tissue}_{sex}/{tissue}_{sex}_DEGs.html")) 
  
  # Heatmap -----------------------------------------------------------------
  
  print(glue::glue("Plotting heatmap of {sex} {tissue} samples"))
  
  voomLogCPM$E[which(rownames(voomLogCPM$E) %in% DEGs$symbol),] %>%
    as.matrix() %>% 
    pheatmap::pheatmap(.,
                       scale = "row",
                       annotation_col = designMatrix %>%
                         tibble::column_to_rownames(var = "Name") %>% 
                         dplyr::select(Treatment, Litter),
                       color = RColorBrewer::brewer.pal(11, name = "RdBu") %>%
                         rev(),
                       show_colnames = FALSE,
                       show_rownames = F,
                       #angle_col = 45,
                       border_color = "grey",
                       main = glue::glue("Z-Scores of {nrow(DEGs)} Differentially Expressed Genes"),
                       fontsize = 16,
                       filename = glue::glue("{tissue}_{sex}/{tissue}_{sex}_heatmap.pdf"),
                       width = 11,
                       height = 8.5,
                       annotation_colors = list(Treatment = c("PCB" = "#F8766D",
                                                             "Control" = "#619CFF")))
 
  # Ontologies and Pathways -------------------------------------------------
  
  print(glue::glue("Performing GO and pathway analysis of {sex} {tissue} samples"))

  tryCatch({
  DEGs %>% 
    dplyr::select(symbol) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2018",
                       "GO_Cellular_Component_2018",
                       "GO_Molecular_Function_2018",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2016",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %T>% # %>% 
    #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %>% 
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";"))) %>% 
    openxlsx::write.xlsx(file = glue::glue("{tissue}_{sex}/{tissue}_{sex}_enrichr.xlsx")) %>%
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = "org.Mm.eg.db",
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("{tissue}_{sex}/{tissue}_{sex}_rrvgo_enrichr.xlsx")) %>%
    DMRichR::GOplot() %>%
    ggplot2::ggsave(glue::glue("{tissue}_{sex}/{tissue}_{sex}_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) },
  error = function(error_condition) {
    print(glue::glue("ERROR: Gene Ontology pipe didn't finish for {sex} {tissue}"))
  })
  print(glue::glue("The pipeline has finished for {sex} {tissue} samples"))
})

sink()

