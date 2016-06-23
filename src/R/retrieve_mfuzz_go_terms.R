#!/usr/bin/env Rscript

library(data.table)

# can do this on one go. load all the cluster results, annotate them in
# phytomine, parse them with data.table

# find mfuzz objects
cluster.files <- list.files(pattern = "mfuzz_object.Rds",
                            recursive = TRUE, full.names = TRUE)
names(cluster.files) <- sapply(cluster.files, function(x)
  strsplit(x, "/", fixed = TRUE)[[1]][3]
)

# read them into R
ReadMfuzzObjects <- function(species, mem.cutoff = 0.7) {
  mfuzz.object <- readRDS(cluster.files[species])
  mfuzz.dt <- data.table(gene = unique(rownames(mfuzz.object$membership)),
                         species = species)
  mfuzz.mem <- mfuzz.object$membership
  clustered.genes <- mfuzz.mem[Biobase::rowMax(mfuzz.mem) > mem.cutoff, ]
  assigned.clusters <- apply(clustered.genes, 1, which.max)
  clusters <- data.table(gene = names(assigned.clusters),
                         cluster = assigned.clusters)
  merge(mfuzz.dt, clusters, by = "gene", all.x = TRUE)
}

# make a data.table of gene names w/ clusters
mfuzz.clusters <- rbindlist(
  lapply(names(cluster.files), ReadMfuzzObjects))
setkey(mfuzz.clusters, gene, species)
mfuzz.clusters <- unique(mfuzz.clusters)

# write to tempfile for python
tmp <- tempfile(fileext = ".txt")
write.table(mfuzz.clusters[, gene], tmp, quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")

# call python script
phytomine.results.file <- system(
  paste("uvb/go_annotations.py --input", tmp), intern = TRUE)

# read the phytomine results
phytomine.results <- data.table(
  read.table(phytomine.results.file, sep = "\t", header = TRUE,
             stringsAsFactors = FALSE, na.strings = "None"),
  key = "primaryIdentifier")
setnames(phytomine.results, "primaryIdentifier", "gene")

# merge with clusters
clusters.go <- merge(mfuzz.clusters, phytomine.results, by = "gene", all = TRUE)

# genes that didn't get annotated
no.go.term <- clusters.go[is.na(organism.shortName)]

# successfully annotated genes
clusters.go.annotated <- clusters.go[!is.na(organism.shortName)]
clusters.go.annotated[, length(unique(gene)), by = species]

# check that right species was returned, e.g. "at" == "A. thaliana"
clusters.go.annotated[, all(species == tolower(gsub(
  "^([[:alpha:]]).*[[:blank:]]+([[:alpha:]]).*", "\\1\\2",
  organism.shortName
)))]
