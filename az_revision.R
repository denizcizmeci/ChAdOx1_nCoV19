setwd("/Users/dz480/Documents/projects/covid/az/")

library(tidyverse)
library(stringr)
library("ggplot2")
library("ggpubr")
library("ggsci")
library(patchwork)
library(scales)
library(rstatix)
library(RColorBrewer)
library("viridis") 
library(ggprism)
library(systemsseRology)
library(DMwR)
library(pheatmap)


load("./data/sero.RData")


col.g1 <- "#a6cee3"
col.g2 <- "#b2df8a"
col.g3 <- "#beaed4"
col.g4 <- "#fdc086"


g.arms <- levels(sero$group)
col.arms <- c(col.g1, col.g2, col.g3, col.g4)


sero.s <- subset(sero, sero$group %in% c("Vaccinated\nCOVID-19-",
                                       "Vaccinated\nCOVID-19+"))

g1.f2a <- ggplot(sero.s, aes(log10(IgG1_WT.S), 
                   log10(FcR2A_WT.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g1.f2a.beta <- ggplot(sero.s, aes(log10(IgG1_B.1.351.S), 
                          log10(FcR2A_B.1.351.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g1.f3a <- ggplot(sero.s, aes(log10(IgG1_WT.S), 
                          log10(FcR3A_WT.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g1.f3a.beta <- ggplot(sero.s, aes(log10(IgG1_B.1.351.S), 
                               log10(FcR3A_B.1.351.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g1.f2b <- ggplot(sero.s, aes(log10(IgG1_WT.S), 
                             log10(FcR2B_WT.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g1.f2b.beta <- ggplot(sero.s, aes(log10(IgG1_B.1.351.S), 
                                  log10(FcR2B_B.1.351.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g1.f3b <- ggplot(sero.s, aes(log10(IgG1_WT.S), 
                             log10(FcR3B_WT.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g1.f3b.beta <- ggplot(sero.s, aes(log10(IgG1_B.1.351.S), 
                                  log10(FcR3B_B.1.351.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))


pdf("IgG1.pdf", height = 5, width = 10)
(g1.f2a | g1.f2b | g1.f3a | g1.f3b) /
  (g1.f2a.beta | g1.f2b.beta | g1.f3a.beta | g1.f3b.beta)
dev.off()




# igg3 --------------------------------------------------------------------


g3.f2a <- ggplot(sero.s, aes(log10(IgG3_WT.S), 
                             log10(FcR2A_WT.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g3.f2a.beta <- ggplot(sero.s, aes(log10(IgG3_B.1.351.S), 
                                  log10(FcR2A_B.1.351.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g3.f3a <- ggplot(sero.s, aes(log10(IgG3_WT.S), 
                             log10(FcR3A_WT.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g3.f3a.beta <- ggplot(sero.s, aes(log10(IgG3_B.1.351.S), 
                                  log10(FcR3A_B.1.351.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g3.f2b <- ggplot(sero.s, aes(log10(IgG3_WT.S), 
                             log10(FcR2B_WT.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g3.f2b.beta <- ggplot(sero.s, aes(log10(IgG3_B.1.351.S), 
                                  log10(FcR2B_B.1.351.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g3.f3b <- ggplot(sero.s, aes(log10(IgG3_WT.S), 
                             log10(FcR3B_WT.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))

g3.f3b.beta <- ggplot(sero.s, aes(log10(IgG3_B.1.351.S), 
                                  log10(FcR3B_B.1.351.S))) + 
  geom_point(size = 2, shape = 21, aes(fill = group))+ 
  stat_cor(method = "spearman", size = 3, label.x.npc = 0.1, aes(color = group))+
  scale_fill_manual(breaks = g.arms, 
                    values = col.arms) +
  scale_color_manual(breaks = g.arms, 
                     values = col.arms) +
  xlim(3,7)+ ylim(3,7)+
  theme_bw()  + 
  theme(aspect.ratio = 1, legend.position = "none", axis.title = element_text(size = 10))


pdf("IgG3.pdf", height = 5, width = 10)
(g3.f2a | g3.f2b | g3.f3a | g3.f3b) /
  (g3.f2a.beta | g3.f2b.beta | g3.f3a.beta | g3.f3b.beta)
dev.off()
