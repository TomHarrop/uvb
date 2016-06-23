#!/usr/bin/env Rscript

library(data.table)
library(goseq)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Perform GO enrichment analsysis on mfuzz clusters")

# find mfuzz objects
GenerateMessage("Loading mfuzz files")
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
GenerateMessage("Combining mfuzz results")
mfuzz.clusters <- rbindlist(
  lapply(names(cluster.files), ReadMfuzzObjects))
setkey(mfuzz.clusters, gene, species)
mfuzz.clusters <- unique(mfuzz.clusters)

# write to tempfile for python
GenerateMessage("Writing mfuzz results to tempfile for python")
tmp <- tempfile(fileext = ".txt")
write.table(mfuzz.clusters[, gene], tmp, quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")

# call python script
GenerateMessage("Finding GO annotation using Phytomine python3 interface")
phytomine.results.file <- system(
  paste("uvb/go_annotations.py --input", tmp), intern = TRUE)

# read the phytomine results
GenerateMessage("Reading phytomine results")
phytomine.results <- data.table(
  read.table(phytomine.results.file, sep = "\t", header = TRUE,
             stringsAsFactors = FALSE, na.strings = "None", quote = ""),
  key = "primaryIdentifier")
setnames(phytomine.results, "primaryIdentifier", "gene")

# merge with clusters
GenerateMessage("Merging phytomine results with mfuzz results")
clusters.go <- merge(mfuzz.clusters, phytomine.results, by = "gene",
                     all = TRUE)

# genes that didn't get annotated
no.phytozome.match <- clusters.go[is.na(organism.shortName)]

# unhelpful go annotations
no.go.namespace <- clusters.go[
  is.na(ontologyAnnotations.ontologyTerm.namespace)]

# successfully annotated genes
gene.ontologies <- c("molecular_function", "biological_process",
                     "cellular_component")
clusters.go.annotated <- clusters.go[
  !is.na(organism.shortName) &
    ontologyAnnotations.ontologyTerm.namespace %in% gene.ontologies]
clusters.go.annotated[, length(unique(gene)), by = species]

# check that right species was returned, e.g. "at" == "A. thaliana"
org.check <- clusters.go.annotated[, all(species == tolower(gsub(
  "^([[:alpha:]]).*[[:blank:]]+([[:alpha:]]).*", "\\1\\2",
  organism.shortName
)))]
if(!org.check) {
  stop("Phytozome script returned non-matching organism short names")
  quit("no", 1)
}


# set up gene ontology enrichment tests
GenerateMessage("Running GO enrichment tests")
RunGoForCluster <- function(cluster.species, cluster.number,
                            cluster.data = clusters.go.annotated){
  GenerateMessage(paste0(cluster.species, " Cluster ", cluster.number))
  # get a list of genes with annotations for cluster.species
  measured.genes <- cluster.data[
    species == cluster.species, unique(gene)]
  # make a vector for cluster.number
  go.vector <- cluster.data[measured.genes,
                            as.numeric(any(cluster == cluster.number)),
                            by = gene][, V1]
  names(go.vector) <- measured.genes
  go.vector[is.na(go.vector)] <- 0
  # set up length vector
  length.data <- cluster.data[measured.genes, length[1], by = gene][,V1]
  # extract category mappings
  category.mapping <- cluster.data[species == cluster.species, data.frame(.(
    gene = gene, category = ontologyAnnotations.ontologyTerm.identifier
  ))]
  # run GO analysis
  pwf <- nullp(DEgenes = go.vector, bias.data = length.data, plot.fit = FALSE)
  data.table(goseq(pwf = pwf, gene2cat = category.mapping),
             species = cluster.species, cluster = cluster.number)
}

# run gene ontology
cluster.go.enrichment <- rbindlist(
  lapply(clusters.go.annotated[, unique(species)], function(x)
    rbindlist(
      lapply(clusters.go.annotated[
        species == x & !is.na(cluster), 1:length(unique(cluster))], function(y)
          RunGoForCluster(x, y)
      ))))

GenerateMessage("Done")

# make output folder
out.dir <- paste0("output/merged/go")
GenerateMessage(paste("Saving output to", out.dir))
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}

# save output
saveRDS(cluster.go.enrichment, paste0(out.dir, "/cluster_go_enrichment.Rds"))
saveRDS(clusters.go.annotated, paste0(out.dir, "/clusters_go_annotated.Rds"))
write.table(cluster.go.enrichment, paste0(out.dir, "/cluster_go_all.tsv"),
            quote = FALSE, sep = "\t", na = "", row.names = FALSE)
write.table(cluster.go.enrichment[over_represented_pvalue < 0.1],
            paste0(out.dir, "/cluster_go_overrepresented.tsv"),
            quote = FALSE, sep = "\t", na = "", row.names = FALSE)

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
writeLines(sInf, paste0(out.dir, "/SessionInfo.txt"))

# exit 0
quit(save = "no", status = 0)
