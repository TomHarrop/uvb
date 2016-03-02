#!/usr/bin/Rscript

library(data.table)
library(DESeq2)

# messages
GenerateMessageText <- function(message.text){
  paste0("[ ", date(), " ]: ", message.text)
}

message(GenerateMessageText("Combine lists of differentially expressed genes"))

# load results
message(GenerateMessageText("Loading DESeq2 results"))

# parse the directory names
SplitPath <- function(x) {
  if (dirname(x)==x) {
    x
  } else {
    c(basename(x), SplitPath(dirname(x)))
  }
}
# load and process files
LoadDeseqResults <- function(pattern) {
  results.files <- list.files(path = "output", pattern = pattern,
                              full.names = TRUE, recursive = TRUE)  
  names(results.files) <- sapply(results.files, function(x) SplitPath(x)[3])
  results.tables <- lapply(results.files, function(x)
    data.table(data.frame(readRDS(x)), keep.rownames = TRUE, key = "rn"))
  rbindlist(results.tables, idcol = "Species")
}
# combine
results.lists <- lapply(list(bs = "bs.Rds", hw = "hw.Rds"), LoadDeseqResults)
deseq.results <- rbindlist(results.lists, idcol = "UV")
setkey(deseq.results, "rn")

# extract "significant" genes
message(GenerateMessageText("Extracting significant genes"))
p.cutoff <- 0.1
sig.genes <- deseq.results[padj <= p.cutoff, unique(rn)]
sig.results <- deseq.results[sig.genes]

# annotate with phytozome
message(GenerateMessageText("Calling Phytomine annotations script"))
# save the sig.genes into a temp file
tmp <- tempfile(fileext = ".txt")
write.table(sig.genes, tmp, quote = FALSE, row.names = FALSE, col.names = FALSE,
            sep = "\t")

# call the python script on the temp file and catch the name of the results file
phytomine.results.file <- system(paste("uvb/phytomine.py --input", tmp),
                                 intern = TRUE)

# read the results in
message(GenerateMessageText("Reading annotation results"))
phytomine.results <- data.table(
  read.table(phytomine.results.file, sep = "\t", header = TRUE,
             stringsAsFactors = FALSE, na.strings = "None"), key = "name")

# join annotations to sig.results
merged.results <- phytomine.results[sig.results, .(
  UV, Species, rn, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
  briefDescription, organism.shortName
)]
setkey(merged.results, "Species", "rn", "UV")

# write merged.results to tsv
message(GenerateMessageText("Saving output"))
out.dir <- "output/merged/deseq2"
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}
out.file <- paste0(out.dir, "/de_genes.tsv")
write.table(merged.results, file = out.file, quote = FALSE, row.names = FALSE,
            sep = "\t", na = "")

# save deseq results RDS
saveRDS(deseq.results, paste0(out.dir, "/deseq_results.Rds"))

# save logs
sInf <- c(paste("git branch:",system("git rev-parse --abbrev-ref HEAD",
                                     intern = TRUE)),
          paste("git hash:", system("git rev-parse HEAD", intern = TRUE)),
          capture.output(sessionInfo()))
logLocation <- paste0(out.dir, "/SessionInfo.de_genes.txt")
writeLines(sInf, logLocation)

message(GenerateMessageText("Done"))
quit(save = "no", status = 0)
