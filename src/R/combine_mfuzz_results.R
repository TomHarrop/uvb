#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
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
GenerateMessage("Loading files")
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
           "S. polyrhiza", "S. moellendorffii",
           "P. patens", "C. reinhardtii"))]
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
GenerateMessage("Making sure we didn't lose any genes")
check <- dim(excluded)[1] + dim(unclustered)[1] +
  dim(clustered)[1] == dim(lfc.table)[1]
if (!check) {
  stop(paste0("Number of genes in lfc.table doesn't match number in plot")) 
  quit(save = "no", status = 1)
}

# draw a mondo plot
GenerateMessage("Producing plot")
pal <- RColorBrewer::brewer.pal(12, "Set3")[-9]
#pal <- clustered[, wesanderson::wes_palette(
#  "Zissou", n = max(cluster), type = "continuous")]

g <- ggplot(mapping = aes(x = bs, y = hw)) +
  theme_minimal(base_size = 10) +
  theme(legend.position = c(5/6, 2/6),
        legend.key.height = unit(3/4, "lines"),
        strip.text = element_text(face = "italic"),
        axis.ticks.length = unit(0, "mm")) + 
  coord_fixed() +
  facet_wrap(~Species, dir = "v") +
  xlab(expression("BS ("*L[2]*"FC)")) +
  ylab(expression("HW ("*L[2]*"FC)")) +
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
GenerateMessage("Saving output")
out.dir <- "output/merged/mfuzz"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

plot.file <- paste0(out.dir, "/cluster_plot.pdf")
pdf(plot.file, width = 7.874, height = 7.874)
g
dev.off()

# tsv for SI
setkey(clusters, "Species", "cluster", "gene")
write.table(clusters[, .(Species, gene, cluster, briefDescription)],
            file = paste0(out.dir, "/clusters.tab"), sep = "\t",
            quote = FALSE, na = "", row.names = FALSE)


# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.mfuzz.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")
quit(save = "no", status = 0)
