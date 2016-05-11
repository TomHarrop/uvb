#!/usr/bin/Rscript

library(data.table)
library(Mfuzz)
library(ggplot2)

set.seed(1)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

# catch species from cli
args <- commandArgs(trailingOnly = TRUE)
if ("-s" %in% args) {
  idx <- which(args == "-s")
  species <- args[idx + 1]
} else {
  stop("Couldn't catch species from input!\n")
  quit(save = "no", status = 1)
}
if (!dir.exists(paste0("output/", species))) {
  stop(paste0("No such species: ", species, "\n")) 
  quit(save = "no", status = 1)
}

GenerateMessage(paste("Running mfuzz for", species))

# load the deseq objects
GenerateMessage("Loading DESeq2 objects")
hw.file <- list.files(paste0("output/", species, "/deseq2"),
                      pattern = "hw", recursive = TRUE,
                      full.names = TRUE)
bs.file <- list.files(paste0("output/", species, "/deseq2"),
                      pattern = "bs", recursive = TRUE,
                      full.names = TRUE)
vst.file <- list.files(paste0("output/", species, "/deseq2"),
                       pattern = "vst", recursive = TRUE,
                       full.names = TRUE)
hw.ds <- readRDS(hw.file)
bs.ds <- readRDS(bs.file)
vst.ds <- readRDS(vst.file)

# make a data.table of hw and bs for each gene
hw <- data.table(data.frame(hw.ds), keep.rownames = TRUE, key = "rn")
bs <- data.table(data.frame(bs.ds), keep.rownames = TRUE, key = "rn")

lfc.table <- hw[bs, .(gene = rn, hw = log2FoldChange, bs = i.log2FoldChange,
                      baseMean)]
setkey(lfc.table, "gene")

# find the top n genes by fc (sum of fold changes)
GenerateMessage("Choosing genes")
n <- 2000
bm.cutoff <- lfc.table[, quantile(baseMean, 0.5)]
message(paste0("            n = ", n))
message(paste0("    bm.cutoff = ", bm.cutoff))

Hypot <- function(a, b) {sqrt((a^2)+(b^2))}

var.genes.table <- lfc.table[
  baseMean > bm.cutoff][
    rev(order(Hypot(hw, bs), na.last = FALSE))][
      1:n]
var.genes <- var.genes.table[, gene]

# take mean across samples
gm_mean <- function(x, na.rm = TRUE) {
  exp(sum(log(x[x > 0]), na.rm = na.rm)/length(x))
}
vst.means.all <- sapply(levels(vst.ds$uvb), function(x)
  apply(SummarizedExperiment::assay(vst.ds[, vst.ds$uvb == x]), 1, gm_mean))

# subset vst values
vst.means <- vst.means.all[var.genes,]

# set up the expression set for mfuzz
GenerateMessage("Setting up Mfuzz object")
p.data <- data.frame(sample = as.factor(colnames(vst.means)),
                     row.names = colnames(vst.means))
pheno.data <- new("AnnotatedDataFrame", data = p.data)
vg.e <- ExpressionSet(assayData = vst.means, data = pheno.data)

# standardise
vg.s <- standardise(vg.e)

# clustering parameters
GenerateMessage("Guessing clustering parameters")
maxClust <- 25
maxFuzz <- 2

# fuzzifier
m <- min(maxFuzz, mestimate(vg.s))
message(paste0("     m = ", m))

# try to guess the cluster number based on inflection of Dmin vs c.
centroids <- data.frame(
  c = 2:(maxClust),
  Dmin = Dmin(vg.s, m = m, crange = seq(2, maxClust, 1), repeats = 1,
              visu = FALSE)
)
points <- seq(2, maxClust, length.out = 10000)
pred <- predict(loess(centroids$Dmin ~ centroids$c), points)
c.est <- points[min(which(diff(pred, differences = 2) < 0))] 
c1 <- c.est %/% 1
message(paste0("    c1 = ", c1))

# diagnostic plot for cluster number
c1.estimate.plot <- ggplot(centroids, aes(x = c, y = Dmin)) +
  geom_smooth() +
  geom_point() +
  geom_point(aes(x = c.est, y = pred[points == c.est]), colour = "red")

# run clustering
GenerateMessage("Running clustering")
clusters <- mfuzz(vg.s, c = c1, m = m)

# assign genes to clusters
GenerateMessage("Parsing clustering results")
mem.cutoff <- 0.7
mem <- clusters$membership
clustered.genes <- mem[rowMax(mem) > mem.cutoff, ]
assigned.clusters <- apply(clustered.genes, 1, which.max)

# annotate clusters with phytozome script
GenerateMessage("Calling phytomine.py to annotate clustered genes")
cluster.table <- data.table(gene = names(assigned.clusters),
                            cluster = assigned.clusters, key = "gene")
tmp <- tempfile(fileext = ".txt")
write.table(cluster.table[, gene], tmp, quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")
phytomine.results.file <- system(paste("uvb/phytomine.py --input", tmp),
                                 intern = TRUE)
phytomine.results <- data.table(
  read.table(phytomine.results.file, sep = "\t", header = TRUE,
             stringsAsFactors = FALSE, na.strings = "None"), key = "name")
cluster.results <- phytomine.results[cluster.table, .(
  gene = name, cluster, organism.shortName, briefDescription
)]
setkey(cluster.results, "gene")

# save output
GenerateMessage("Saving output")
out.dir <- paste0("output/", species, "/mfuzz")
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
so <- function(obj, lab){
  saveRDS(object = obj, file = paste0(out.dir, "/", lab, ".Rds"))
}
so(clusters, "mfuzz_object")
so(cluster.results, "annotated_clusters")
so(lfc.table, "lfc_table")
so(var.genes.table, "var_genes_table")
ggsave(paste0(out.dir, "/c_plot.pdf"), c1.estimate.plot)

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.mfuzz.txt")
writeLines(sInf, logLocation)

GenerateMessage("Done")
quit(save = "no", status = 0)

# 
# # plot lfc with clusters
# 
# pal <- RColorBrewer::brewer.pal(12, "Set3")[-9]
# pd <- cluster.results[lfc.table,
#                       .(gene, cluster, hw, bs)]
# ggplot(mapping = aes(x = bs, y = hw)) +
#   theme_minimal() +
#   #scale_colour_manual(values = pal) +
#   geom_point(data = pd[!gene %in% var.genes.table$gene],
#              colour = "grey90", alpha = 0.2, size = 1) +
#   geom_point(data = var.genes.table[!gene %in% names(assigned.clusters)],
#              colour = "black", alpha = 0.5, size = 1) +
#   geom_point(data = pd[!is.na(cluster)],
#              mapping = aes(colour = as.factor(cluster)), alpha = 0.8)



