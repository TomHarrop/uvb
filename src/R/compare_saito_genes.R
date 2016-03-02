library(data.table)
library(ggplot2)

# messages
GenerateMessage <- function(message.text){
  message(paste0("[ ", date(), " ]: ", message.text))
}

GenerateMessage("Compare expression of homologs between species")

# load combined DESeq2 results
GenerateMessage("Loading DESeq2 results table")
deseq.results.file <- "output/merged/deseq2/deseq_results.Rds"
if (!file.exists(deseq.results.file)) {
  stop("Couldn't find deseq_results.Rds")
}
deseq.results<- readRDS(deseq.results.file)

# load Saito et al. data (doi:10.1016/j.plaphy.2013.02.001, 
# http://dx.doi.org/10.1016/j.plaphy.2013.02.001)
GenerateMessage("Parsing data from saitoT1.csv")
raw.data.file <- "data/saitoT1.csv"
if (!file.exists(raw.data.file)) {
  stop("Couldn't find saitoT1.csv")
}
raw.data <- data.table(read.csv(raw.data.file, skip = 2, fill = TRUE,
                                stringsAsFactors = FALSE))

# retrieve the wonky rows
shifted.rows <- raw.data[, grep("\\[[[:digit:]]+\\]", Locus)]

# tidy up stepwise
raw.data[shifted.rows, `Ref.` := Locus]
raw.data[shifted.rows, Locus := Transparent.testa.mutant]
raw.data[shifted.rows, Gene.designation := Protein.function]
raw.data[shifted.rows, c("Protein.function", "Transparent.testa.mutant") := NA]

# the merged rows should be the ones above the wonky ones. Fill Protein.function
# from above
j <- which(names(raw.data) == "Protein.function")
for (i in shifted.rows) {
  set(raw.data, i, j, raw.data[i - 1, Protein.function])
}
rm(j)

# set blank cells to NA
for (j in names(raw.data)) {
  set(raw.data, which(raw.data[[j]] == ""), j, NA)
}

# uppercase the Locus column for phytozome searches
raw.data[, Locus := toupper(Locus)]

saito.t1 <- copy(raw.data)
setkey(saito.t1, "Locus")

######################
# PHYTOMINE HOMOLOGS #
######################

# write At loci to tempfile for python
tmp <- tempfile(fileext = ".txt")
write.table(saito.t1[, Locus], tmp, quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")

# call python script to retrieve homologs
GenerateMessage("Calling Phytomine script to retrieve homologs")
phytomine.results.file <- system(paste("uvb/retrieve_homologs.py --input", tmp),
                                 intern = TRUE)

# read the phytomine results
GenerateMessage("Reading Phytomine results")
phytomine.results <- data.table(
  read.table(phytomine.results.file, sep = "\t", header = TRUE,
             stringsAsFactors = FALSE, na.strings = "None"),
  key = "primaryIdentifier")

# merge phytozome results to saito table
GenerateMessage("Generating data for plotting")
flavonol.homologs <- phytomine.results[saito.t1, .(
  At.id = primaryIdentifier,
  At.function = Protein.function,
  Saito.ref = `Ref.`,
  homolog.species = gsub("^([[:upper:]]).*[[:blank:]]([[:lower:]]).*", "\\1\\2",
                         homolog.gene2.organism.shortName),
  homolog.id = homolog.gene2.primaryIdentifier,
  homolog.relationship
)]
setkey(flavonol.homologs, "At.id")

# add rows for at (at homologs are the original gene)
flavonol.homologs <- rbind(flavonol.homologs,
                           unique(flavonol.homologs[, .(
                             At.id, At.function, Saito.ref,
                             homolog.species = "At",
                             homolog.id = At.id,
                             homolog.relationship = "identical"
                           )]))

# remove rows where no homologs were found
flavonol.homologs <- flavonol.homologs[!is.na(homolog.id)]

# join expression data
GenerateMessage("Joining LFC and p-values")
setkey(flavonol.homologs, "homolog.id")
flavonol.expression <- deseq.results[flavonol.homologs, .(
  At.id, At.function, Saito.ref, homolog.species, homolog.id,
  homolog.relationship, UV, log2FoldChange, padj
)]

# go wide to plot hs vs bw
plot.data
reshape2::dcast(flavonol.expression, homolog.id ~ UV, value.var = "log2FoldChange")
flavonol.expression[is.na(UV)]

# visualise
ggplot(flavonol.expression, aes(x = ))
