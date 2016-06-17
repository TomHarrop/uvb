#!/usr/bin/Rscript

library(ggplot2)
library(data.table)
source("src/R/facet_adjust.R")

star.logs <- readRDS("output/mapping_stats/starLogs.Rds")

# tidy
star.logs[, Species := substr(Library, 1, 2)]
star.logs[, Sample := substr(Library, 4, 8)]
star.logs[, Sample := gsub("^uv", "", Sample)]
star.logs[!is.na(as.numeric(Sample)), Sample := paste(Sample, "nm")]
star.logs[, Species := factor(
  Species, levels = c("at", "sl", "os", "sp", "sm", "pp", "cr"))]
star.logs[, Species := plyr::mapvalues(
  Species,
  from = c("at", "sl" ,"os", "sp", "sm", "pp", "cr"),
  to = c("A. thaliana", "S. lycopersicum", "O. sativa",
         "S. polyrhiza", "S. moellendorffii",
         "P. patens", "C. reinhardtii"))]

# set colours
wespal <- wesanderson::wes_palette("Rushmore", 5)[-2]

panel.a <- ggplot(star.logs, aes(x = `Number of input reads` / 1000000,
                      y = `Uniquely mapped reads %`,
                      colour = Sample, group = Species)) +
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size = rel(1), face = "bold", hjust = 0),
        legend.position = c(7/8,1/4),
        legend.justification = c(0.5,0.5),
        axis.ticks.length = unit(0, "mm"),
        strip.text = element_text(face = "italic")) +
  xlab("Input reads (M)") + ylab("Uniquely mapped reads (%)") +
  facet_wrap(~Species, ncol = 4) +
  geom_smooth(method = "lm", se = FALSE, size = 0.25, alpha = 0.5,
              show.legend = FALSE, colour = "black") +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = wespal,
                     guide = guide_legend(title = "Treatment")) +
  ggtitle("A")

panel.b <- ggplot(star.logs, aes(x = `Number of input reads` / 1000000,
                             y = `Number of reads in genes` / 1000000,
                             colour = Sample, group = Species)) +
  xlab("Input reads (M)") + ylab("Reads mapped to genes (M)") +
  theme_minimal(base_size = 10) + 
  theme(plot.title = element_text(size = rel(1), face = "bold", hjust = 0),
        legend.position = c(7/8,1/4),
        legend.justification = c(0.5,0.5),
        axis.ticks.length = unit(0, "mm"),
        strip.text = element_text(face = "italic")) +
  facet_wrap(~Species, ncol = 4) +
  geom_hline(yintercept = 20, linetype = 5, size = 0.25, alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 0.25, alpha = 0.5,
              show.legend = FALSE, colour = "black") +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = wespal,
                     guide = guide_legend(title = "Treatment")) +
  ggtitle("B")

pdf("output/mapping_stats/Figure\ S1.pdf", width = 7.874, height = 9.843,
    pointsize = 10)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
facet_adjust(panel.a, "up", vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
facet_adjust(panel.b, "up", vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()
