#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# messages
GenerateMessageText <- function(message.text){
  paste0("[ ", date(), " ]: ", message.text)
}

# parse the directory names
SplitPath <- function(x) {
  if (dirname(x)==x) {
    x
  } else {
    c(basename(x), SplitPath(dirname(x)))
  }
}
# load files
LoadMfuzzResults <- function(pattern) {
  results.files <- list.files(path = "output", pattern = pattern,
                              full.names = TRUE, recursive = TRUE)  
  names(results.files) <- sapply(results.files, function(x) SplitPath(x)[3])
  results.tables <- lapply(results.files, readRDS)
  rbindlist(results.tables, idcol = "Species")
}

lfc.table <- LoadMfuzzResults("lfc_table.Rds")
cluster.results <- LoadMfuzzResults("annotated_clusters.Rds")
var.genes.table <- LoadMfuzzResults("var_genes_table.Rds")

# tidy the data: rename species column, order species column, 
TidyPlotData <- function(results.table){
  my.dt <- copy(results.table)
  my.dt[, Species := factor(
    Species, levels = c("at", "sp", "pp", "sl", "sm", "cr", "os"))]
  my.dt[, Species := plyr::mapvalues(
    Species,
    from = c("at", "sl" ,"os", "sp", "sm", "pp", "cr"),
    to = c("A. thaliana", "S. lycopersicum", "O. sativa",
           "S. polyrhiza", "S. moellendorffi",
           "P. patens", "C. rheinhardtii"))]
  my.dt
}

lfc <- TidyPlotData(lfc.table)
clusters <- TidyPlotData(cluster.results)
var.genes <- TidyPlotData(var.genes.table)

setkey(lfc, "gene")
setkey(clusters, "gene")
setkey(var.genes, "gene")

# set up plots
excluded <- lfc[!gene %in% var.genes$gene]
unclustered <- var.genes[!gene %in% clusters$gene]
clustered <- lfc[clusters, .(Species, gene, bs, hw, cluster)]

# check we got everything
check <- dim(excluded)[1] + dim(unclustered)[1] +
  dim(clustered)[1] == dim(lfc.table)[1]
if (!check) {
  stop(paste0("Number of genes in lfc.table doesn't match number in plot")) 
  quit(save = "no", status = 1)
}

# draw a mondo plot
pal <- RColorBrewer::brewer.pal(12, "Set3")[-9]
g <- ggplot(mapping = aes(x = bs, y = hw)) +
  theme_minimal(base_size = 10) +
  theme(legend.position = c(5/6, 2/6),
        strip.text = element_text(face = "italic")) + 
  coord_fixed() +
  facet_wrap(~Species, dir = "v") +
  xlab(expression(Log[2]*'-fold change (BS)')) +
  ylab(expression(Log[2]*'-fold change (HW)')) +
  geom_point(data = excluded,
             colour = "#D9D9D9", alpha = 0.2, size = 0.5) +
  geom_point(data = unclustered,
             colour = "black", alpha = 0.5, size = 0.5) +
  geom_point(data = clustered,
             mapping = aes(colour = as.factor(cluster)),
             size = 1, alpha = 0.5) +
  scale_colour_manual(values = pal,
                      guide = guide_legend(title = "Cluster"))

# save output
message(GenerateMessageText("Saving output"))
out.dir <- "output/merged/mfuzz"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

plot.file <- paste0(out.dir, "/cluster_plot.pdf")
pdf(plot.file, width = 3.937 * 2, height = 3.937 * 2)
g
dev.off()
