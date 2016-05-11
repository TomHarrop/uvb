library(data.table)
library(ggplot2)
library(gtable)

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

# the merged rows should be the ones above the wonky ones. Fill
# Protein.function from above
j <- which(names(raw.data) == "Protein.function")
for (i in shifted.rows) {
  set(raw.data, i, j, raw.data[i - 1, Protein.function])
}
rm(j)

# deal with PAL1/PAL2 (these genes have the same homologs)
raw.data[Gene.designation %in% c("PAL1", "PAL2"),
         Gene.designation := "PAL1/PAL2"]

# set blank cells to NA
for (j in names(raw.data)) {
  set(raw.data, which(raw.data[[j]] == ""), j, NA)
}

# uppercase the Locus column for phytozome searches
raw.data[, Locus := toupper(Locus)]

# label synthesis and modification
raw.data[grep("^Flav.*transferase$", Protein.function),
         pathway := "Flavonol modification"]
raw.data[grep("^Anthocyanin.*transferase", Protein.function),
         pathway := "Anthocyanin modification"]
raw.data[Gene.designation %in% c("GSTF12", "TT12", "AHA10"),
         pathway := "Translocation"]
raw.data[is.na(pathway), pathway := "Synthesis"]

# set gene order
raw.data[, Gene.designation := factor(
  Gene.designation, levels = unique(Gene.designation))]

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
phytomine.results.file <- system(
  paste("uvb/retrieve_homologs.py --input", tmp), intern = TRUE)

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
  At.designation = Gene.designation,
  Pathway = pathway,
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
                             At.id, At.function, At.designation, Pathway,
                             Saito.ref,
                             homolog.species = "At",
                             homolog.id = At.id,
                             homolog.relationship = "identical"
                           )]))

# set order of species
flavonol.homologs[, homolog.species := factor(
  homolog.species, levels = c("At", "Sl" ,"Os", "Sp", "Sm", "Pp", "Cr"))]

# order genes by species
flavonol.homologs[, s.order := as.numeric(homolog.species)]
setkey(flavonol.homologs, "s.order", "homolog.id")
flavonol.homologs[, homolog.id := factor(
  homolog.id, levels = unique(homolog.id))]
flavonol.homologs[, s.order := NULL]

# remove rows where no homologs were found
flavonol.homologs <- flavonol.homologs[!is.na(homolog.id)]
setkey(flavonol.homologs, "homolog.id")

# join expression data
GenerateMessage("Joining LFC and p-values")
flavonol.expression <- deseq.results[flavonol.homologs, .(
  At.id, At.function, At.designation, Pathway, Saito.ref, homolog.species,
  homolog.id, homolog.relationship, UV, log2FoldChange, lfcSE
)]

# order genes by species
flavonol.expression[, s.order := as.numeric(homolog.species)]
setkey(flavonol.expression, "s.order", "homolog.id")
flavonol.expression[, homolog.id := factor(
  homolog.id, levels = rev(unique(homolog.id)))]
flavonol.expression[, s.order := NULL]

setkey(flavonol.expression, "homolog.id", "UV")

# plot with a categorical x-axis
pd <- flavonol.expression[!is.na(log2FoldChange) &
                            Pathway == "Synthesis"]

# split into two columns
keep <- c("PAL1/PAL2", "C4H", "4CL3", "ACC1", "CHS", "CHI")
pd[At.designation %in% keep, col := 1]
pd[is.na(col), col := 2]

# fix strip labels (convert to expression syntax)
QuoteString <- function(x) {paste0("italic('", x, "')")}
pd[, At.designation := plyr::mapvalues(
  At.designation, from = levels(At.designation),
  to = QuoteString(levels(At.designation)))]
pd[, At.designation := plyr::mapvalues(
  At.designation,
  from = "italic('F3′H (CYP75B1)')",
  to = "italic('F3'*minute*'H')")]
pd[, At.designation := plyr::mapvalues(
  At.designation,
  from = unique(At.designation),
  to =  gsub("/", "')*/*italic('", unique(At.designation)))]


# fix UV labels
pd[, UV := toupper(UV)]

# function for plotting each column
PlotColumn <- function(colnum) {
  ggplot(pd[col == colnum],
         aes(y = homolog.id, x = log2FoldChange, colour = homolog.species)) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text.y = element_text(angle = 0, face = "italic", hjust = 0),
      axis.text.y = element_text(face = "italic"),
      axis.ticks.length	= grid::unit(0, "mm"),
      axis.text.x = element_text(),
      axis.title.x = element_text(),
      legend.text = element_text(face = "italic")) +
    scale_colour_brewer(palette = "Set1", guide = guide_legend(title=NULL)) +
    facet_grid(At.designation ~ UV, scales = "free_y", space = "free_y",
               drop = TRUE, labeller = label_parsed) +
    ylab(NULL) + xlab(expression(Log[2]*"-"*fold~change)) +
    xlim(-2.2, 2.2) +
    geom_vline(xintercept = 0, size = 0.2) +
    geom_vline(xintercept = c(-log(1.5, 2), log(1.5, 2)),
               size = 0.2, linetype = 2) +
    geom_errorbarh(aes(xmax = log2FoldChange + lfcSE,
                       xmin = log2FoldChange - lfcSE),
                   height = 0.3, size = 0.1, colour = "black") +
    geom_point(size = 0.5)
}

# plot left column with no legend or axis title
c1 <- PlotColumn(1)
c1grob <- ggplotGrob(c1)
legend <- gtable_filter(c1grob, "guide-box")
legend.loc <- c1grob$layout[c1grob$layout$name == "guide-box",c("l","r")]

xlab <- gtable_filter(c1grob, "xlab")
xlab.loc <- c1grob$layout[c1grob$layout$name == "xlab",c("t","b")]

c1.nolegend <- c1grob[-c(min(xlab.loc), max(xlab.loc)),
                      -c(min(legend.loc), max(legend.loc))]

# plot right column and extract legend (ggtable)
c2 <- PlotColumn(2)
c2grob <- ggplotGrob(c2)
legend.loc <- c2grob$layout[c2grob$layout$name == "guide-box",c("l","r")]
xlab.loc <- c2grob$layout[c2grob$layout$name == "xlab",c("t","b")]

c2.nolegend <- c2grob[-c(min(xlab.loc), max(xlab.loc)),
                      -c(min(legend.loc), max(legend.loc))]

# grid back together (left, right, legend)
plot.with.legend <- gridExtra::arrangeGrob(
  grobs = list(c1.nolegend, c2.nolegend, legend),
  widths = list(0.45, 0.45, 0.1), ncol = 3)

# add x-axis label
# add two lines, one for whitespace and one for the images
combined.figure.grob <- gtable_add_rows(
  plot.with.legend,
  grid::unit.c( # bit of whitespace?
    unit(12, "pt"), # 12pt for axis label 
    unit(6, "pt"))) # bottom margin
fig.with.xaxis <- gtable_add_grob(
  combined.figure.grob, xlab, t=2, b=2,l=1,r=2)

# draw
pdf("/home/tom/Dropbox/temp/homologs.pdf", width = 3.937 * 2, height = 9.843)
grid::grid.newpage()
grid::grid.draw(fig.with.xaxis)
dev.off()










#########
# SKIP #
#########

# some homologs not found in deseq results???
bof <- flavonol.expression[is.na(log2FoldChange), unique(homolog.id)]
bof %in% deseq.results[, rn]
deseq.results[rn %in% bof]

flavonol.expression[homolog.id == "Cre07.g322884"]
# -> no, they just weren't detected (0 reads)

flavonol.expression[homolog.species == "Os" & UV == "hw", ][23,]
flavonol.expression[homolog.id == "AT2G37040"]
flavonol.expression[homolog.id == "LOC_Os04g43800"]
flavonol.expression[padj < 0.1]

# go wide to plot hs vs bw
plot.data <- reshape2::dcast(
  unique(flavonol.expression), homolog.id + homolog.species + At.id +
    At.function + homolog.relationship ~ UV, value.var = "log2FoldChange")

# visualise
ggplot(plot.data, aes(x = bs, y = hw, label = homolog.id)) +
  facet_wrap(~homolog.species) +
  geom_point() +
  geom_text()


ggplot(plot.data, aes(x = bs, y = hw, label = homolog.id, colour = homolog.species)) +
  facet_wrap(~At.id) +
  scale_color_brewer(palette = "Set1") +
  geom_point()



pd2 <- plot.data[plot.data$homolog.relationship %in% c("one-to-one", "identical"), ]
pd2[is.na(pd2)] <- 0
ggplot(pd2, aes(x = bs, y = hw, label = homolog.id, colour = homolog.species)) +
  facet_wrap(~At.id) +
  scale_color_brewer(palette = "Set1") +
  geom_point()


# error bars

# pick two genes
bof.wide <- deseq.results[rn %in% c("LOC_Os04g43800", "AT3G51240")]
bof.long <- reshape2::melt(bof.wide, id.vars = c("rn", "UV"),
                           measure.vars = c("log2FoldChange", "lfcSE"))
bof.long[, variable := paste(UV, variable, sep = ".")]
bof <- reshape2::dcast(bof.long, rn ~ variable)

# plot
ggplot(bof, aes(x = bs.log2FoldChange, y =  hw.log2FoldChange)) +
  geom_point() +
  geom_errorbarh(aes(xmax = bs.log2FoldChange + bs.lfcSE,
                     xmin = bs.log2FoldChange - bs.lfcSE),
                 height = 0) +
  geom_errorbar(aes(ymax = hw.log2FoldChange + hw.lfcSE,
                    ymin = hw.log2FoldChange - hw.lfcSE),
                width = 0)

# add errorbars (DESeq2 lfcSE) to plot
# flavonol.expression <- deseq.results[flavonol.homologs, .(
#   At.id, At.function, Saito.ref, homolog.species, homolog.id,
#   homolog.relationship, UV, log2FoldChange, lfcSE
# )]
# setkey(flavonol.expression, "homolog.id", "UV")

# melt the l2fc and se columns
f.exp.long <- reshape2::melt(
  flavonol.expression,
  id.vars = c("At.id", "homolog.species", "homolog.id", "UV"),
  measure.vars = c("log2FoldChange", "lfcSE"))

# make new labels for variable
f.exp.long[, measure := paste(UV, variable, sep = ".")]

# cast by variable
f.exp.plot <- data.table(reshape2::dcast(
  f.exp.long,
  homolog.id + At.id + homolog.species ~ measure,
  value.var = "value"))

# plot as before
ggplot(f.exp.plot,
       aes(x = bs.log2FoldChange, y =  hw.log2FoldChange,
           colour = homolog.species)) +
  facet_wrap(~At.id, ncol = 7) +
  theme_grey(base_size = 8) +
  coord_fixed() +
  geom_errorbarh(aes(xmax = bs.log2FoldChange + bs.lfcSE,
                     xmin = bs.log2FoldChange - bs.lfcSE),
                 height = 0, size = 0.5, colour = "grey50") +
  geom_errorbar(aes(ymax = hw.log2FoldChange + hw.lfcSE,
                    ymin = hw.log2FoldChange - hw.lfcSE),
                width = 0, size = 0.5, colour = "grey50") +
  geom_point(alpha = 0.8, size = 1) +
  scale_colour_brewer(palette = "Set1")
ggsave("/home/tom/Dropbox/temp/homologs.pdf", width = 10, height = 7.5)


