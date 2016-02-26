#!/usr/bin/Rscript

library(ggplot2)
library(data.table)

star.logs <- readRDS("output/mapping_stats/starLogs.Rds")

star.logs[, Species := substr(Library, 1, 2)]
star.logs[, Sample := substr(Library, 4, 8)]

ggplot(star.logs, aes(x = `Number of input reads` / 1000000,
                      y = `Number of reads in genes` / 1000000,
                      colour = Sample, group = Species)) +
  coord_fixed() +
  xlab("Millions of input reads") + ylab("Millions of reads mapped to genes") +
  theme(legend.position = c(7/8,1/4),
        legend.justification = c(0.5,0.5)) +
  facet_wrap(~Species, ncol = 4) +
  geom_hline(yintercept = 20, linetype = 5, size = 0.25, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, size = 0.25, alpha = 0.5,
              show.legend = FALSE, colour = "black") +
  geom_point(alpha = 0.75) +
  scale_color_brewer(palette = "Set1", guide = guide_legend(title = NULL))

ggsave("output/mapping_stats/readsInGenes.pdf", width = 10, height = 7.5)

ggplot(star.logs, aes(x = `Number of input reads` / 1000000,
                      y = `Uniquely mapped reads %`,
                      colour = Sample, group = Species)) +
  theme(legend.position = c(7/8,1/4),
        legend.justification = c(0.5,0.5)) +
  xlab("Millions of input reads") + ylab("Uniquely mapped reads (%)") +
  facet_wrap(~Species, ncol = 4) +
  geom_smooth(method = "lm", se = TRUE, size = 0.25, alpha = 0.5,
              show.legend = FALSE, colour = "black") +
  geom_point(alpha = 0.75) +
  scale_color_brewer(palette = "Set1")

ggsave("output/mapping_stats/mappingPercentage.pdf", width = 10, height = 7.5)



