#!/usr/bin/Rscript

library(data.table)
library(ggplot2) # remove later? only for pca plot
library(scales) # remove later? only for density plot

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

# load the quant files
count.data.files <- list.files(paste0("output/", species, "/star"),
                               pattern = "ReadsPerGene", recursive = TRUE,
                               full.names = TRUE)
count.data.raw <- lapply(count.data.files, read.table, stringsAsFactors = FALSE)
names(count.data.raw) <- gsub(".ReadsPerGene.out.tab", "",
                              basename(count.data.files), fixed = TRUE)

# get the unstranded counts column
extractUnstrandedCounts <- function(x){
  x <- copy(x)
  data.frame(row.names = x$V1, x$V2)
}
unstranded.count.data <- lapply(count.data.raw, extractUnstrandedCounts)
count.data <- do.call(cbind, unstranded.count.data)
names(count.data) <- names(unstranded.count.data)

# tidy up gene names
TrimGeneId <- function(x, species = species) {
  if (species %in% c("at", "sp")) {
    gsub("([^\\.]+).*", "\\1", x)
  } else if (species %in% c("cr", "sl")) {
    gsub("([^\\.]+\\.[^\\.]+).*", "\\1", x)
  } else if (species == "os") {
    gsub("(.*)\\.MSU.*", "\\1", x)
  } else if (species == "sm") {
    gsub("^g", "", x)
  } else if (species == "pp") {
    x
  } else {
    stop("Species not matched in TrimGeneId")
  }
}
# at: AT1G01010.TAIR10
# TrimGeneId("AT1G01010.TAIR10", "at")
# cr: Cre01.g000017.v5.5
# TrimGeneId("Cre01.g000017.v5.5", "cr")
# os: LOC_Os01g01010.MSUv7.0
# TrimGeneId("LOC_Os01g01010.MSUv7.0", "os")
# TrimGeneId("ChrSy.fgenesh.gene.5.MSUv7.0", "os")
# TrimGeneId("ChrUn.fgenesh.gene.90.MSUv7.0", "os")
# pp: Pp3c1_20
# TrimGeneId("Pp3c1_20", "pp")
# sl: Solyc00g005000.2.iTAGv2.3
# TrimGeneId("Solyc00g005000.2.iTAGv2.3", "sl")
# sm: g401996
# TrimGeneId("g401996", "sm")
# sp: Spipo0G0000100.v2
# TrimGeneId("Spipo0G0000100.v2", "sp")

rownames(count.data) <- sapply(rownames(count.data), TrimGeneId,
                               species = species)

# set up colData
meta.data.table <- data.table(rn = colnames(count.data))
meta.data.table[, sample := gsub("^uv([[:alnum:]]+)\\.[[:digit:]]+$",
                                 "\\1", rn)]
meta.data.table[sample %in% c("PAR", "320"), uvb := "none"]
meta.data.table[sample == "295", uvb := "bs"]
meta.data.table[sample == "305", uvb := "hw"]

meta.data <- data.frame(meta.data.table, row.names = "rn")
meta.data$uvb <- factor(meta.data$uvb, levels = c("none", "bs", "hw"))

# set up DESeq2 object
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count.data[!grepl("^N_", rownames(count.data)), ],
  colData = meta.data, design = ~ uvb)  
dds <- DESeq2::DESeq(dds)

######
# QC #
######

rld <- DESeq2::rlogTransformation(dds, blind = TRUE)
transformed.counts <- SummarizedExperiment::assay(rld)

# choose some genes with high counts
qS <- quantile(rowSums(DESeq2::counts(dds)), 0.7)
qM <- quantile(rowMeans(DESeq2::counts(dds)), 0.7)
expressed.genes <- rownames(DESeq2::counts(dds)[
  rowSums(DESeq2::counts(dds)) > qS | rowMeans(DESeq2::counts(dds)) > qM, ])
expressed.tc <- transformed.counts[expressed.genes,]

# plot pca
pca <- prcomp(t(expressed.tc))
percentVar <- round(pca$sdev^2/sum(pca$sdev^2) * 100, 1)

pcaPlotData <- data.frame(
  label = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  UVB = SummarizedExperiment::colData(rld)$uvb
)
bump <- 0.025 * (range(pcaPlotData$PC2)[2] - range(pcaPlotData$PC2)[1])
pca.plot <- ggplot(pcaPlotData, aes(x = PC1, y = PC2, label = label, colour = UVB)) +
  guides(label=FALSE) +
  xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
  ylab(paste0("PC2 (", percentVar[2], "% variance)")) +
  geom_point() +
  geom_text(mapping = aes(y = PC2 + bump), show.legend = FALSE)

# plot density
expressed.counts <- DESeq2::counts(dds)[expressed.genes,]
density.plot.data <- reshape2::melt(expressed.counts)
density.plot.data$col <- substr(density.plot.data$Var2, 1, 5)

density.plot <- ggplot(density.plot.data, aes(x = value, fill = col))  +
  xlab("Transformed counts") + ylab(NULL) +
  guides(fill = FALSE, colour = FALSE) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  scale_fill_brewer(palette = "Set1") +
  geom_density(alpha = 0.8, colour = NA) +
  facet_wrap(~Var2, ncol = 2)

# downstream (blind = FALSE)
vst <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
rld <- DESeq2::rlogTransformation(dds, blind = FALSE)

# results
bs <- DESeq2::results(dds, contrast = c("uvb", "bs", "none"),
                      lfcThreshold = log2(1.5))
# subset(data.frame(bs), padj < 0.1)
# DESeq2::plotCounts(dds, "AT3G56290", intgroup = "uvb")

hw <- DESeq2::results(dds, contrast = c("uvb", "hw", "none"),
                      lfcThreshold = log2(1.5))
# subset(data.frame(hw), padj < 0.1)
# DESeq2::plotCounts(dds, "LOC_Os04g43800", intgroup = "uvb")

# save output

# make output folder
outDir <- paste0("output/", species, "/deseq2")
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

saveRDS(dds, paste0(outDir, "/dds.Rds"))
saveRDS(rld, paste0(outDir, "/rld.Rds"))
saveRDS(vst, paste0(outDir, "/vst.Rds"))
saveRDS(bs, paste0(outDir, "/bs.Rds"))
saveRDS(hw, paste0(outDir, "/hw.Rds"))
ggsave(paste0(outDir, "/pca_plot.pdf"), pca.plot)
ggsave(paste0(outDir, "/density_plot.pdf"), density.plot)

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
writeLines(sInf, paste0(outDir, "/SessionInfo.txt"))

# exit 0
quit(save = "no", status = 0)
