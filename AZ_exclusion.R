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

dem <- read.csv("./data/dem.csv")

load("/Users/dz480/Documents/projects/covid/az copy 2/data/sero.RData")
sero_old <- sero

load("/Users/dz480/Documents/projects/covid/az copy 2/excluded/sero.RData")
vaccinated <- sero

pbo <- subset(sero_old, sero_old$group %in% c("Placebo_control",    "Placebo_case"))

placebo <- subset(pbo, (pbo$infection.status == "case" & pbo$strain ==  "Beta") | 
                    (pbo$infection.status == "control" & is.na(pbo$cov_res) ) )
                    

pbo.dem <- merge(dem, placebo, by.x = "pid", by.y = "id.y")
pbo.dem$group <- as.vector(pbo.dem$group.y)
pbo.dem$group[which(pbo.dem$group.y == "Placebo_control")] <- "Placebo\nCOVID-19-"
pbo.dem$group[which(pbo.dem$group.y == "Placebo_case")] <- "Placebo\nCOVID-19+"
pbo.dem$group[which(pbo.dem$group.y == "Vaccine_control")] <- "Vaccinated\nCOVID-19-"
pbo.dem$group[which(pbo.dem$group.y == "Vaccine_case")] <- "Vaccinated\nCOVID-19+"
pbo.dem$group <- factor(pbo.dem$group, levels = c("Placebo\nCOVID-19-", "Placebo\nCOVID-19+",
                                                  "Vaccinated\nCOVID-19-","Vaccinated\nCOVID-19+" ))

pbo.dem <- pbo.dem[,c(1,233, 2:232)]


sero_corrected <- rbind(pbo.dem, vaccinated)

sero <- sero_corrected
sero <- sero[order(sero$group),]
save(sero, file = "./data/sero.RData")
write.csv(sero, file = "./data/AZ_systemsserology_091421.csv", row.names = FALSE)



col.g1 <- "#a6cee3"
col.g2 <- "#b2df8a"
col.g3 <- "#beaed4"
col.g4 <- "#fdc086"




# pbo and vac ------------------------------------------------------------------

g.arms <- levels(sero$group)
col.arms <- c(col.g1, col.g2, col.g3, col.g4)


stat.wil <- long %>%
  group_by(feature, assay, antigen, type) %>%
  wilcox_test(value ~ group)
library(tidyverse)
long <- gather(sero, "feature", "value", -c(1:57))

library(stringr)
df.feature <- data.frame(feature = colnames(sero)[58:ncol(sero)], 
                         assay = str_split_fixed(colnames(sero)[58:ncol(sero)], "_",2)[,1],
                         antigen = str_split_fixed(colnames(sero)[58:ncol(sero)], "_",2)[,2])

df.feature$antigen <- factor(df.feature$antigen, levels = c("WT.S", "B.1.1.7.S" ,"B.1.351.S", "P1.S" , "B.617.1.S" , "B.617.2.S",
                                                            "WT.RBD", "B.1.1.7.RBD" ,"B.1.351.RBD", "P1.RBD" , "B.617.1.RBD" , "B.1.617.2.RBD",
                                                            "WT.S1"     ,    "WT.S2",     "WT.NTD","WT.N"   , 
                                                            "HKU1.Spike", "OC43.Spike", "RSV", "CMV", "Ebola", "Influenza"))


df.feature$assay <- factor(df.feature$assay, levels = c("IgG1", "IgG3","IgM","IgA", "FcR2A", "FcR2B", "FcR3A", "FcR3B"))
df.feature$type <- "Titer"
df.feature$type[which(df.feature$assay %in% c("FcR2A", "FcR2B", "FcR3A", "FcR3B"))]  <- "FcR"
df.feature$type <- factor(df.feature$type, levels = c("Titer", "FcR"))

long <- merge(long, df.feature, by = "feature")


save(long, file = "long.RData") 


stat.wil$p.star <- symnum(stat.wil$p, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))
stat.wil$q <- p.adjust(stat.wil$p, method = "BH")
stat.wil$q.star <- symnum(stat.wil$q, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))


temp.s <- subset(long, long$antigen == "WT.S")

# df_text <- subset(stat.wil, stat.wil$antigen == "WT.S" & stat.wil$q < 0.05)[3,]
# df_text$y.position <- 6.1
# 
# df_text2 <- df_text[, c(1:7, 16, 17)]
# colnames(df_text2)[8] <- "p.signif"


pdf("TitersFcRS_WT.S_pboANDvaccine.pdf", width = 14, height = 7)
ggplot(temp.s, aes(group, log10(value)))+
  geom_violin(scale = "width", aes(fill = group),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          axis.text.y = element_text( color = "black",size = 14),
                          axis.title.y = element_text( color = "black",size = 16),
                          legend.position =  "none", aspect.ratio = 1) + 
  #add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
  #          colour = "black")+
  # add_pvalue(df_text2, label.size = 5)+
  # scale_x_discrete(labels = c("Vaccinated\nCOVID-19-_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19-_Male"=  "Vaccinated\ninfected",
  #                             "Vaccinated\nCOVID-19+_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19+_Male" = "Vaccinated\nuninfected"))+
  #stat_compare_means(label = "p.signif", size = 4,  label.x = 1.5)+
  ylim(3.3,6.3)+
  xlab(" ") +ylab("SARS-CoV-2 WT\nlog10 MFI")+ facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()


#  vac wt  ------------------------------------------------------------------

long.vac <- droplevels(subset(long, long$group %in% levels(sero$group)[3:4]))

g.arms <- levels(sero$group)[3:4]
col.arms <- c(col.g3, col.g4)


stat.wil <- long.vac %>%
  group_by(feature, assay, antigen, type) %>%
  wilcox_test(value ~ group)

stat.wil$p.star <- symnum(stat.wil$p, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))
stat.wil$q <- p.adjust(stat.wil$p, method = "BH")
stat.wil$q.star <- symnum(stat.wil$q, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))


temp.s <- subset(long.vac, long.vac$antigen == "WT.S")

df_text <- subset(stat.wil, stat.wil$antigen == "WT.S" & stat.wil$q < 0.05)
df_text$y.position <- 6.1
df_text2 <- df_text[, c("feature","assay","antigen","type",".y.","group1","group2",
                        "q.star","y.position")]
colnames(df_text2)[which(colnames(df_text2) == "q.star")] <- "p.signif"


pdf("Figure2B_TitersFcRS_WT.S.pdf", width = 10, height = 7)
ggplot(temp.s, aes(group, log10(value)))+
  geom_violin(scale = "width", aes(fill = group),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          axis.text.y = element_text( color = "black",size = 14),
                          axis.title.y = element_text( color = "black",size = 16),
                          legend.position =  "none", aspect.ratio = 1) + 
  add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,label.size = 4,
           colour = "black")+
  ylim(3.3,6.3)+
  xlab(" ") +ylab("SARS-CoV-2 WT\nlog10 MFI")+ facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()


temp.s <- subset(long.vac, long.vac$antigen == "B.1.351.S")

df_text <- subset(stat.wil, stat.wil$antigen == "B.1.351.S" & stat.wil$q < 0.05)
df_text$y.position <- 6.1
df_text2 <- df_text[, c("feature","assay","antigen","type",".y.","group1","group2",
                        "q.star","y.position")]
colnames(df_text2)[which(colnames(df_text2) == "q.star")] <- "p.signif"


pdf("Figure2C_TitersFcRS_Beta.pdf", width = 10, height = 7)
ggplot(temp.s, aes(group, log10(value)))+
  geom_violin(scale = "width", aes(fill = group),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          axis.text.y = element_text( color = "black",size = 14),
                          axis.title.y = element_text( color = "black",size = 16),
                          legend.position =  "none", aspect.ratio = 1) + 
  add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,label.size = 4,
             colour = "black")+
  ylim(3.3,6.3)+
  xlab(" ") +ylab("Beta B.1.351 S\nlog10 MFI")+ facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()

####RBD
long.vac <- droplevels(subset(long, long$group %in% levels(sero$group)[3:4]))

g.arms <- levels(sero$group)[3:4]
col.arms <- c(col.g3, col.g4)


stat.wil <- long.vac %>%
  group_by(feature, assay, antigen, type) %>%
  wilcox_test(value ~ group)

stat.wil$p.star <- symnum(stat.wil$p, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))
stat.wil$q <- p.adjust(stat.wil$p, method = "BH")
stat.wil$q.star <- symnum(stat.wil$q, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))


temp.s <- subset(long.vac, long.vac$antigen == "WT.RBD")

df_text <- subset(stat.wil, stat.wil$antigen == "WT.RBD" & stat.wil$q < 0.05)
df_text$y.position <- 6.1
df_text2 <- df_text[, c("feature","assay","antigen","type",".y.","group1","group2",
                        "q.star","y.position")]
colnames(df_text2)[which(colnames(df_text2) == "q.star")] <- "p.signif"


pdf("Figuresupp_2B_TitersFcRS_WT.RBD.pdf", width = 10, height = 7)
ggplot(temp.s, aes(group, log10(value)))+
  geom_violin(scale = "width", aes(fill = group),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          axis.text.y = element_text( color = "black",size = 14),
                          axis.title.y = element_text( color = "black",size = 16),
                          legend.position =  "none", aspect.ratio = 1) + 
  add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,label.size = 4,
             colour = "black")+
  ylim(3.3,6.3)+
  xlab(" ") +ylab("WT RBD\nlog10 MFI")+ facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()


temp.s <- subset(long.vac, long.vac$antigen == "B.1.351.RBD")

df_text <- subset(stat.wil, stat.wil$antigen == "B.1.351.RBD" & stat.wil$q < 0.05)
df_text$y.position <- 6.1
df_text2 <- df_text[, c("feature","assay","antigen","type",".y.","group1","group2",
                        "q.star","y.position")]
colnames(df_text2)[which(colnames(df_text2) == "q.star")] <- "p.signif"


pdf("Figuresupp_2C_TitersFcRRBD_Beta.pdf", width = 10, height = 7)
ggplot(temp.s, aes(group, log10(value)))+
  geom_violin(scale = "width", aes(fill = group),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          axis.text.y = element_text( color = "black",size = 14),
                          axis.title.y = element_text( color = "black",size = 16),
                          legend.position =  "none", aspect.ratio = 1) + 
  add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,label.size = 4,
             colour = "black")+
  ylim(3.3,6.3)+
  xlab(" ") +ylab("Beta B.1.351 RBD\nlog10 MFI")+ facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()




###






temp.s <- droplevels(subset(long.vac, long.vac$antigen %in% c("WT.S", "B.1.351.S") &
                   long.vac$assay == "IgG1"))

pdf("Figure1_IgG1.pdf", width = 6, height = 4)
ggplot(temp.s, aes(group, log10(value)))+
  geom_violin(scale = "width", aes(fill = group),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          axis.text.y = element_text( color = "black",size = 14),
                          axis.title.y = element_text( color = "black",size = 16),
                          legend.position =  "none", aspect.ratio = 1) + 
  ylim(3.3,6.3)+
  xlab(" ") +ylab("IgG1 log10 MFI")+ facet_wrap(.~antigen, ncol = 2, scales = "free_x",
                                                labeller=as_labeller(c("WT.S" = "SARS-CoV-2 WT",
                                                              "B.1.351.S" = "Beta B.1.351 S")))
dev.off()



# heatmap -----------------------------------------------------------------
sero.vac <- droplevels(subset(sero, sero$group %in% levels(sero$group)[3:4]))

my.df <- sero.vac[order(sero.vac$group),]

row.names(my.df) <- my.df$ids

X <- as.matrix(my.df[,c(58:ncol(my.df))])
log.ind <- which(grepl("ADCD|IgG|IgA|IgM|FcR", colnames(X)))
X[,log.ind] <- log10(X[,log.ind])
# log.ind <- which(grepl("IgG|IgA|IgM|FcR", colnames(X)))
# X[,log.ind] <- log10(X[,log.ind])
X <- X[, which(!grepl("RSV|CMV|Ebola|Influenza", colnames(X)))]
X <- X[, which(grepl("WT.S|B.1.351.S|WT.RBD|B.1.351.RBD", colnames(X)))]

X <- scale(X, center = TRUE, scale = TRUE)
#X <- knnImputation(X)


antigens <- c("WT.S", "B.1.351.S", "WT.RBD", "B.1.351.RBD", "WT.S1", "WT.S2")
as <-  c("IgG1", "IgG3", "IgM", "IgA", "FcR2A", "FcR3A", "FcR2B", "FcR3B")
col.labels <- paste0(rep(as, each = length(antigens)), "_", antigens)
#col.labels <- paste0(as, "_", rep(antigens, each = length(as)))


heat.m <- X[, col.labels]

ann_colors <- list(group = c("Vaccinated\nCOVID-19-" = col.g3,
                             "Vaccinated\nCOVID-19+" = col.g4))

row.names(heat.m) <- row.names(my.df)
myBreaks <- seq(min(heat.m, na.rm = TRUE), max(heat.m, na.rm = TRUE), length.out = 100)
#myBreaks <- seq(-4, 5, length.out = 100)
myColor <- inferno(99)
pdf("Figure2A_heatmap_allAZ_SandSAvarRBDS1S2_byassay.pdf", height = 20, width = 35)
print(pheatmap(heat.m,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               #show_colnames = FALSE,
               show_rownames = FALSE,
               #labels_row = c(),
               annotation_row = my.df[,c("group"), drop = FALSE],
               #annotation_names_row = FALSE,
               #annotation_legend = FALSE,
               annotation_colors = ann_colors,
               color=myColor, breaks=myBreaks,
               fontsize_row = 12,
               fontsize_col = 8,
               #angle_col = 45,
               gaps_col = seq(6, ncol(heat.m), 6),
               gaps_row = c(140),
               #gaps_row = c(2,4,6,8,10,12,14),
               cellheight = 6,
               cellwidth = 10
               #gaps_col = c(10,seq(20,ncol(heat),20)),
))
dev.off()



# variants ----------------------------------------------------------------

median.long.vac <- long.vac %>%
  group_by(feature, assay, antigen, type, group) %>%
  summarise(value = median(value, na.rm = TRUE))

median.long.vac <- as.data.frame(median.long.vac)

stat.group <- long.vac %>%
  group_by(antigen, assay) %>%
  wilcox_test(value ~ group)

stat.group$p.star <- symnum(stat.group$p, corr = FALSE, na = FALSE, 
                            cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                            symbols = c("***", "**", "*", " "))
stat.group$q <- p.adjust(stat.group$p, method = "BH")
stat.group$q.star <- symnum(stat.group$q, corr = FALSE, na = FALSE, 
                            cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                            symbols = c("***", "**", "*", " "))

long.vac <- long.vac[order(long.vac$pid),]
stat.antigen <- long.vac %>%
  group_by(assay, type, group) %>%
  wilcox_test(value ~ antigen, paired = TRUE)

stat.antigen$p.star <- symnum(stat.antigen$p, corr = FALSE, na = FALSE, 
                              cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                              symbols = c("***", "**", "*", " "))
stat.antigen$q <- p.adjust(stat.antigen$p, method = "BH")
stat.antigen$q.star <- symnum(stat.antigen$q, corr = FALSE, na = FALSE, 
                              cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                              symbols = c("***", "**", "*", " "))



antigens <- c("WT.S", "B.1.351.S", "WT.RBD", "B.1.351.RBD") 
  
  
temp <- droplevels(subset(long.vac, long.vac$antigen %in% antigens))
#temp.s$value <- log10(temp.s$value)

temp.median <- droplevels(subset(median.long.vac, median.long.vac$antigen %in% antigens ))

temp.median1 <- droplevels(subset(median.long.vac, median.long.vac$antigen %in% antigens[1:2] ))
temp.median2 <- droplevels(subset(median.long.vac, median.long.vac$antigen %in% antigens[3:4] ))

group_text <- subset(stat.group, stat.group$q < 0.05 &  stat.group$antigen %in% antigens)


antigen_text <- subset(stat.antigen, stat.antigen$q < 0.05 & 
                       stat.antigen$group1 %in% antigens & stat.antigen$group2 %in% antigens)

antigen_text <- subset(antigen_text, (antigen_text$group1 == "WT.S" & antigen_text$group2 == "B.1.351.S") |
                         (antigen_text$group1 == "B.1.351.S" & antigen_text$group2 == "WT.S") |
                         (antigen_text$group1 == "WT.RBD" & antigen_text$group2 == "B.1.351.RBD") |
                         (antigen_text$group1 == "B.1.351.RBD" & antigen_text$group2 == "WT.RBD") )

antigen_text$y.position <- 6.3

antigen_text <- antigen_text[,c("group","assay", "type",".y.", "group1", "group2",
                                "q.star", "y.position")]
colnames(antigen_text)[which(colnames(antigen_text) == "q.star") ] <- "p.signif"


pdf("Figure4A_univariates_S_RBD_wtbeta.pdf", width = 10, height = 7)
ggplot(temp, aes(antigen, log10(value))) + 
  geom_point( shape = 21, size = 2, aes(fill = group),position = position_dodge(width = 0.8))+
  geom_line(data = temp.median1, aes( group = group, color =  group))+
  geom_line(data = temp.median2, aes( group = group, color =  group))+
  geom_text(data = group_text, mapping= aes(antigen, y = 6, label = q.star),
                       hjust   = 0.5,
                       vjust   = 0)+
  add_pvalue(antigen_text, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
             colour = "group", label.size = 3.6)+
  scale_fill_manual(breaks = g.arms, values = col.arms)+
  scale_colour_manual(breaks = g.arms, values = col.arms)+
  #ylim(3,6.7)+  
  ylab("log10 MFI")+ xlab("")+ #coord_flip()+
  theme_bw() + theme( plot.title = element_text(hjust = 0.5, size = 8),
                      legend.position = "top",
                      #aspect.ratio = 1,
                      axis.text.x = element_text(angle = 45, hjust=1),
                      strip.background=element_rect(fill="white"))+
  facet_wrap( ~ assay, ncol = 4, scales = "free_x")
dev.off()


### Figure S3


antigens <- levels(long.vac$antigen)[1:6]


temp <- droplevels(subset(long.vac, long.vac$antigen %in% antigens))
#temp.s$value <- log10(temp.s$value)

temp.median <- droplevels(subset(median.long.vac, median.long.vac$antigen %in% antigens ))


group_text <- subset(stat.group, stat.group$q < 0.05 &  stat.group$antigen %in% antigens)


antigen_text <- subset(stat.antigen, stat.antigen$q < 0.05 & 
                         stat.antigen$group1 %in% antigens & stat.antigen$group2 %in% antigens)

antigen_text <- subset(antigen_text, (antigen_text$group1 == "WT.S") |
                         (antigen_text$group2 == "WT.S") )

antigen_text$y.position <- 6.3

antigen_text <- antigen_text[,c("group","assay", "type",".y.", "group1", "group2",
                                "q.star", "y.position")]
colnames(antigen_text)[which(colnames(antigen_text) == "q.star") ] <- "p.signif"


pdf("FigureS3A_univariates_S_variants.pdf", width = 10, height = 7)
ggplot(temp, aes(antigen, log10(value))) + 
  geom_point( shape = 21, size = 2, aes(fill = group),position = position_dodge(width = 0.8))+
  geom_line(data = temp.median, aes( group = group, color =  group))+
  geom_text(data = group_text, mapping= aes(antigen, y = 6, label = q.star),
            hjust   = 0.5,
            vjust   = 0)+
  add_pvalue(antigen_text, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
             colour = "group", label.size = 3.6)+
  scale_fill_manual(breaks = g.arms, values = col.arms)+
  scale_colour_manual(breaks = g.arms, values = col.arms)+
  #ylim(3,6.7)+  
  ylab("log10 MFI")+ xlab("")+ #coord_flip()+
  theme_bw() + theme( plot.title = element_text(hjust = 0.5, size = 8),
                      legend.position = "top",
                      #aspect.ratio = 1,
                      axis.text.x = element_text(angle = 45, hjust=1),
                      strip.background=element_rect(fill="white"))+
  facet_wrap( ~ assay, ncol = 4, scales = "free_x")
dev.off()




antigens <- levels(long.vac$antigen)[7:12]


temp <- droplevels(subset(long.vac, long.vac$antigen %in% antigens))
#temp.s$value <- log10(temp.s$value)

temp.median <- droplevels(subset(median.long.vac, median.long.vac$antigen %in% antigens ))


group_text <- subset(stat.group, stat.group$q < 0.05 &  stat.group$antigen %in% antigens)


antigen_text <- subset(stat.antigen, stat.antigen$q < 0.05 & 
                         stat.antigen$group1 %in% antigens & stat.antigen$group2 %in% antigens)

antigen_text <- subset(antigen_text, (antigen_text$group1 == "WT.RBD") |
                         (antigen_text$group2 == "WT.RBD") )

antigen_text$y.position <- 6.3

antigen_text <- antigen_text[,c("group","assay", "type",".y.", "group1", "group2",
                                "q.star", "y.position")]
colnames(antigen_text)[which(colnames(antigen_text) == "q.star") ] <- "p.signif"


pdf("FigureS3B_univariates_RBD_variants.pdf", width = 10, height = 7)
ggplot(temp, aes(antigen, log10(value))) + 
  geom_point( shape = 21, size = 2, aes(fill = group),position = position_dodge(width = 0.8))+
  geom_line(data = temp.median, aes( group = group, color =  group))+
  geom_text(data = group_text, mapping= aes(antigen, y = 6, label = q.star),
            hjust   = 0.5,
            vjust   = 0)+
  add_pvalue(antigen_text, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
             colour = "group", label.size = 3.6)+
  scale_fill_manual(breaks = g.arms, values = col.arms)+
  scale_colour_manual(breaks = g.arms, values = col.arms)+
  #ylim(3,6.7)+  
  ylab("log10 MFI")+ xlab("")+ #coord_flip()+
  theme_bw() + theme( plot.title = element_text(hjust = 0.5, size = 8),
                      legend.position = "top",
                      #aspect.ratio = 1,
                      axis.text.x = element_text(angle = 45, hjust=1),
                      strip.background=element_rect(fill="white"))+
  facet_wrap( ~ assay, ncol = 4, scales = "free_x")
dev.off()




# correlations ------------------------------------------------------------

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)

sero.vac.control <- subset(sero.vac, sero.vac$group == "Vaccinated\nCOVID-19-")
sero.vac.case <- subset(sero.vac, sero.vac$group == "Vaccinated\nCOVID-19+")

X.control <- sero.vac.control[,58:ncol(sero.vac.control)]
X.case <- sero.vac.case[,58:ncol(sero.vac.case)]

# X <- X[, which(!grepl("Ebola|CMV|RSV|Influenza", colnames(X)))]
# X <- X[,-c(6,24,42,60,78,96,114,132)]

X.control <- X.control[, which(!grepl("RSV|CMV|Ebola|Influenza|WT.N|S1|S2", colnames(X.control)))]
X.control <- X.control[, which(grepl("WT.S|B.1.351.S|WT.RBD|B.1.351.RBD", colnames(X.control)))]

X.case <- X.case[, which(!grepl("RSV|CMV|Ebola|Influenza|WT.N|S1|S2", colnames(X.case)))]
X.case <- X.case[, which(grepl("WT.S|B.1.351.S|WT.RBD|B.1.351.RBD", colnames(X.case)))]


library(corrplot)
library(RColorBrewer)
library(psych)


# X <- as.matrix(my.df[,c(41:ncol(my.df))])
# # X <- X[, which(!grepl("RSV|CMV|Ebola|Influenza|S1|S2", colnames(X)))]
# # X <- X[, which(grepl("WT.S|B.1.351.S|WT.RBD|B.1.351.RBD", colnames(X)))]
antigens <- c("WT.S", "B.1.351.S", "WT.RBD", "B.1.351.RBD")
as <-  c("IgG1", "IgG3", "IgM", "IgA", "FcR2A", "FcR3A", "FcR2B", "FcR3B")
col.labels <- paste0(rep(as, each = length(antigens)), "_", antigens)
X.control <- X.control[,col.labels]

a <- corr.test(X.control, y = NULL, use = "pairwise",method="spearman",adjust="BH", 
               alpha=.05,ci=TRUE,minlength=5)
pdf("Figure4B_Vaccine_control.pdf", width = 8, height = 8)
print(corrplot(a$r, p.mat = a$p, method = "color",  
               #type = "upper",
               sig.level = c(.001, .01, .05), 
               pch.cex = .8,
               diag = FALSE,
               insig = "label_sig", pch.col = "white", order = "original", 
               col=rev(brewer.pal(n=8, name="PuOr")),
               tl.cex = .5, tl.col="black", tl.srt=45))
dev.off()

X.case <- X.case[,col.labels]
a <- corr.test(X.case, y = NULL, use = "pairwise",method="spearman",adjust="BH", 
               alpha=.05,ci=TRUE,minlength=5)
pdf("Figure4C_Vaccine_case.pdf", width = 8, height = 8)
print(corrplot(a$r, p.mat = a$p, method = "color",  
               #type = "upper",
               sig.level = c(.001, .01, .05), 
               pch.cex = .8,
               diag = FALSE,
               insig = "label_sig", pch.col = "white", order = "original", 
               col=rev(brewer.pal(n=8, name="PuOr")),
               tl.cex = .5, tl.col="black", tl.srt=45))
dev.off()






# 
# res <- rcorr(as.matrix(X.case), type = "spearman")
# 
# resr <- as.data.frame(res$r)
# resr$selfeat <- row.names(resr)
# df.r <- resr
# 
# resp <- as.data.frame(res$P)
# resp$selfeat <- row.names(resp)
# df.p <- resp
# 
# 
# dfr <- df.r %>% gather(feature, rho, -c(selfeat))
# dfp <- df.p %>% gather(feature, p, -c(selfeat))
# dfp$q <- p.adjust(dfp$p, method = "BH")
# 
# 
# identical(dfr$selfeat, dfp$selfeat)
# identical(dfr$feature, dfp$feature)
# 
# df <- cbind(dfr, dfp[,c("p","q")])
# df$antigen <- as.vector(str_split_fixed(df$feature, "_", 2)[,2])
# 
# df.case <- df
# df.control <- df
# 
# 
# r.st.tp1 <- droplevels(subset(r.st, r.st$day == 1))
# 
# tp1.ptext <- spread(r.st.tp1[c("Antigen", "Assay", "group2","p.adjstar")], 
#                     Assay, p.adjstar)
# ## all non sig
# 
# tp1.fc <- spread(st.tp1[c("Antigen", "Assay", "log2FC.PCV13_PPSV23")], 
#                  Assay, log2FC.PCV13_PPSV23)
# row.names(tp2.fc) <- tp2.fc$Antigen
# tp2.fc <- tp2.fc[,-1]
# 
# pdf("./figures/tp2.fc.pdf", width = 7, height = 4)
# pheatmap(tp2.fc,
#          cluster_rows = FALSE, 
#          color=brewer.pal(n=8, name="PuOr"),
#          breaks = seq(-max(abs(tp2.fc), na.rm = TRUE), max(abs(tp2.fc), na.rm = TRUE),
#                       length.out = 9),
#          cluster_cols = FALSE, 
#          display_numbers = tp2.ptext, 
#          gaps_col = c(8, 12), 
#          angle_col = "45",
#          cellheight=18, cellwidth = 18,
#          number_color = "grey30", fontsize_number = 12,
#          legend_breaks = c(-4, -2, 0, 2, 4, max(abs(tp2.fc), na.rm = TRUE)), 
#          main = "29-43 Days \nPost-vaccination", legend_labels = c("-4", "-2", "0", "2", "4", "log2 (PCV13/PPSV23)"),
#          na_col = "grey")
# dev.off()
# 
# 

# gender ------------------------------------------------------------------

col.d1 <- "#efedf5"
col.d2 <- "#756bb1"
col.d3 <- "#fff7bc"
col.d4 <- "#d95f0e"

temp <- long

temp$comp <- paste0(temp$group, "\n", temp$sc_gender)
temp$comp <- factor(temp$comp)

g.arms <- levels(temp$comp)
col.arms <- c(col.d1, col.d2, col.d3, col.d4)


stat.wil <- temp %>%
  group_by(feature, assay, antigen, type) %>%
  wilcox_test(value ~ comp)

stat.wil$p.star <- symnum(stat.wil$p, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))
stat.wil$q <- p.adjust(stat.wil$p, method = "BH")
stat.wil$q.star <- symnum(stat.wil$q, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))


temp.s <- subset(temp, temp$antigen == "WT.S")

df_text <- subset(stat.wil, stat.wil$antigen == "WT.S" & stat.wil$q < 0.05)[3,]
df_text$y.position <- 6.1

df_text2 <- df_text[, c(1:7, 16, 17)]
colnames(df_text2)[8] <- "p.signif"


pdf("TitersFcRS_WT.S_gender.pdf", width = 14, height = 7)
ggplot(temp.s, aes(comp, log10(value)))+
  geom_violin(scale = "width", aes(fill = comp),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          legend.position =  "none", aspect.ratio = 1) + 
  #add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
  #          colour = "black")+
  # add_pvalue(df_text2, label.size = 5)+
  # scale_x_discrete(labels = c("Vaccinated\nCOVID-19-_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19-_Male"=  "Vaccinated\ninfected",
  #                             "Vaccinated\nCOVID-19+_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19+_Male" = "Vaccinated\nuninfected"))+
  #stat_compare_means(label = "p.signif", size = 4,  label.x = 1.5)+
  ylim(3.3,6.3)+
  xlab(" ") +ylab(" ")+ facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()


# bmi ------------------------------------------------------------------

col.d1 <- "#f2f0f7"
col.d2 <- "#cbc9e2"
col.d3 <- "#9e9ac8"
col.d4 <- "#6a51a3"
col.d5 <- "#ffffd4"
col.d6 <- "#fed98e"
col.d7 <- "#fe9929"
col.d8 <- "#cc4c02"

temp <- long

temp$comp <- paste0(temp$group, "\n", temp$bmi)
temp$comp <- factor(temp$comp, levels = c("Vaccinated\nCOVID-19-\nUnderweight",
                                          "Vaccinated\nCOVID-19-\nNormal",
                                          "Vaccinated\nCOVID-19-\nOverweight",
                                          "Vaccinated\nCOVID-19-\nObese",
                                          "Vaccinated\nCOVID-19+\nUnderweight",
                                          "Vaccinated\nCOVID-19+\nNormal",
                                          "Vaccinated\nCOVID-19+\nOverweight",
                                          "Vaccinated\nCOVID-19+\nObese"))

g.arms <- levels(temp$comp)
col.arms <- c(col.d1, col.d2, col.d3, col.d4, col.d5, col.d6, col.d7, col.d8)


stat.wil <- temp %>%
  group_by(feature, assay, antigen, type) %>%
  wilcox_test(value ~ comp)

stat.wil$p.star <- symnum(stat.wil$p, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))
stat.wil$q <- p.adjust(stat.wil$p, method = "BH")
stat.wil$q.star <- symnum(stat.wil$q, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))


temp.s <- subset(temp, temp$antigen == "WT.S")

df_text <- subset(stat.wil, stat.wil$antigen == "WT.S" & stat.wil$q < 0.05)
df_text$y.position <- 6.1

df_text2 <- df_text[, c(1:7, 16, 17)]
colnames(df_text2)[8] <- "p.signif"


pdf("TitersFcRS_WT.S_bmi.pdf", width = 20, height = 12)
ggplot(temp.s, aes(comp, log10(value)))+
  geom_violin(scale = "width", aes(fill = comp),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          legend.position =  "none", aspect.ratio = 1) + 
  #add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
  #          colour = "black")+
  # add_pvalue(df_text2, label.size = 5)+
  # scale_x_discrete(labels = c("Vaccinated\nCOVID-19-_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19-_Male"=  "Vaccinated\ninfected",
  #                             "Vaccinated\nCOVID-19+_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19+_Male" = "Vaccinated\nuninfected"))+
  #stat_compare_means(label = "p.signif", size = 4,  label.x = 1.5)+
  ylim(3.3,6.3)+
  xlab(" ") +ylab(" ")+ facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()



# pca ---------------------------------------------------------------------

data <- sero.vac

X <- data[,58:ncol(data)]

X <- X[, which(!grepl("Ebola|CMV|RSV|Influenza", colnames(X)))]
X <- X[,-c(6,24,42,60,78,96,114,132)]

#X <- X[, which(!grepl("RSV|CMV|Ebola|Influenza|WT.N|S1|S2", colnames(X)))]
#X <- X[, which(grepl("WT.S|B.1.351.S|WT.RBD|B.1.351.RBD", colnames(X)))]


log.ind <- which(grepl("ADCD|IgG|IgA|IgM|FcR", colnames(X)))
X[,log.ind] <- log10(X[,log.ind])

X <- scale(X, center = TRUE, scale = TRUE)
#X <- knnImputation(X)
y <- factor(data$group)



# PCA using "prcomp" functions
pca.temp <- prcomp(X , scale = FALSE, center = FALSE)
pcs <- as.data.frame(pca.temp$x)

# labels to use to colour and shape samples 
#pcs$protection <- y
pcs$sex <- data$sc_gender

col.d1 <- "#ffff99"
#col.d2 <- "#756bb1"
col.d3 <- "#386cb0"
#col.d4 <- "#d95f0e"

col.sex <- c(col.d1, col.d3)
g.sex <- c("Female","Male")


# Calculate percentage of variance
PoV <- pca.temp$sdev^2/sum(pca.temp$sdev^2)

# col.nonprotected <- "blue"
# col.protected <- "yellow"

scores <- ggplot(pcs, aes(x=PC1, y=PC2)) +
  geom_point(aes( fill = sex), shape = 21, size = 3)+
  #scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = col.sex, breaks = g.sex) +
  scale_color_manual(values = col.arms, breaks = g.arms) +
  geom_hline(yintercept = 0, linetype = 2, colour = "slategray")+
  geom_vline(xintercept = 0, linetype = 2, colour = "slategray")+
  theme_classic()+
  xlab(paste0("PC1 (",as.character(signif(PoV[1]*100,3)),"%)")) +
  ylab(paste0("PC2 (",as.character(signif(PoV[2]*100,3)),"%)")) + ggtitle(" ")+
  theme(aspect.ratio=1, legend.title=element_blank())
pdf("PCA_sex.pdf", width = 5, height = 5)
scores
dev.off()


pcs$age <- data$age_calc
pdf("PCA_age.pdf", width = 5, height = 5)
ggplot(pcs, aes(x=PC1, y=PC2)) +
  geom_point(aes( fill = age), shape = 21, size = 3)+
  #scale_shape_manual(values = c(21, 24))+
  #scale_fill_manual(values = col.sex, breaks = g.sex) +
  #scale_color_manual(values = col.arms, breaks = g.arms) +
  geom_hline(yintercept = 0, linetype = 2, colour = "slategray")+
  geom_vline(xintercept = 0, linetype = 2, colour = "slategray")+
  theme_classic()+
  xlab(paste0("PC1 (",as.character(signif(PoV[1]*100,3)),"%)")) +
  ylab(paste0("PC2 (",as.character(signif(PoV[2]*100,3)),"%)")) + ggtitle(" ")+
  theme(aspect.ratio=1)
dev.off()



pcs$bmi <- data$bmi
pdf("PCA_bmi.pdf", width = 5, height = 5)
ggplot(pcs, aes(x=PC1, y=PC2)) +
  geom_point(aes( fill = bmi), shape = 21, size = 3)+
  #scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
  #scale_color_manual(values = col.arms, breaks = g.arms) +
  geom_hline(yintercept = 0, linetype = 2, colour = "slategray")+
  geom_vline(xintercept = 0, linetype = 2, colour = "slategray")+
  theme_classic()+
  xlab(paste0("PC1 (",as.character(signif(PoV[1]*100,3)),"%)")) +
  ylab(paste0("PC2 (",as.character(signif(PoV[2]*100,3)),"%)")) + ggtitle(" ")+
  theme(aspect.ratio=1, legend.title=element_blank())
dev.off()


pcs$race <- data$sc_race
pdf("PCA_race.pdf", width = 5, height = 5)
ggplot(pcs, aes(x=PC1, y=PC2)) +
  geom_point(aes( fill = race), shape = 21, size = 3)+
  #scale_shape_manual(values = c(21, 24))+
  scale_fill_manual(values = c("#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4")) +
  #scale_color_manual(values = col.arms, breaks = g.arms) +
  geom_hline(yintercept = 0, linetype = 2, colour = "slategray")+
  geom_vline(xintercept = 0, linetype = 2, colour = "slategray")+
  theme_classic()+
  xlab(paste0("PC1 (",as.character(signif(PoV[1]*100,3)),"%)")) +
  ylab(paste0("PC2 (",as.character(signif(PoV[2]*100,3)),"%)")) + ggtitle(" ")+
  theme(aspect.ratio=1, legend.title=element_blank())
dev.off()



temp <- subset(long.vac, long.vac$feature == "FcR3B_WT.S")

col.d2 <- "#ffff99"
col.d1 <- "#e6f5c9"
col.d4 <- "#386cb0"
col.d3 <- "#b3cde3"

temp$comp <- paste0(temp$group, "\n", temp$sc_gender)
temp$comp <- factor(temp$comp)

g.arms <- levels(temp$comp)
col.arms <- c(col.d1, col.d2, col.d3, col.d4)



pdf("FcR3B_WT.S_gender.pdf", width = 5, height = 4)
ggplot(temp, aes(comp, log10(value)))+
  geom_violin(scale = "width", aes(fill = comp),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          legend.position =  "none", aspect.ratio = 1) + 
  #add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
  #          colour = "black")+
  # add_pvalue(df_text2, label.size = 5)+
  # scale_x_discrete(labels = c("Vaccinated\nCOVID-19-_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19-_Male"=  "Vaccinated\ninfected",
  #                             "Vaccinated\nCOVID-19+_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19+_Male" = "Vaccinated\nuninfected"))+
  #stat_compare_means(label = "p.signif", size = 4,  label.x = 1.5)+
  ylim(3.3,6.3)+
  xlab(" ")+
  ylab("FcR3B_WT.S log10 MFI")
# + 
#   facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()



pdf("FcR3B_WT.S_age.pdf", width = 5, height = 4)
ggplot(temp, aes(age_calc, log10(value)))+ geom_point()+
  #geom_violin(scale = "width", aes(fill = comp),alpha = 0.8, draw_quantiles = c(0.5))+ 
  #geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          legend.position =  "none", aspect.ratio = 1) + 
  #add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
  #          colour = "black")+
  # add_pvalue(df_text2, label.size = 5)+
  # scale_x_discrete(labels = c("Vaccinated\nCOVID-19-_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19-_Male"=  "Vaccinated\ninfected",
  #                             "Vaccinated\nCOVID-19+_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19+_Male" = "Vaccinated\nuninfected"))+
  #stat_compare_means(label = "p.signif", size = 4,  label.x = 1.5)+
  ylim(3.3,6.3)+
  xlab(" ")+
  ylab("FcR3B_WT.S log10 MFI")
# + 
#   facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()


pdf("FcR3B_WT.S_bmi.pdf", width = 5, height = 4)
ggplot(temp, aes(bmi, log10(value)))+
  geom_violin(scale = "width", aes(fill = bmi),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          legend.position =  "none", aspect.ratio = 1) + 
  #add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
  #          colour = "black")+
  # add_pvalue(df_text2, label.size = 5)+
  # scale_x_discrete(labels = c("Vaccinated\nCOVID-19-_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19-_Male"=  "Vaccinated\ninfected",
  #                             "Vaccinated\nCOVID-19+_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19+_Male" = "Vaccinated\nuninfected"))+
  #stat_compare_means(label = "p.signif", size = 4,  label.x = 1.5)+
  ylim(3.3,6.3)+
  xlab(" ")+
  ylab("FcR3B_WT.S log10 MFI")
# + 
#   facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()




pdf("FcR3B_WT.S_race.pdf", width = 5, height = 4)
ggplot(temp, aes(sc_race, log10(value)))+
  geom_violin(scale = "width", aes(fill = sc_race),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = c("#b3e2cd", "#fdcdac", "#cbd5e8", "#f4cae4")) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          legend.position =  "none", aspect.ratio = 1) + 
  #add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
  #          colour = "black")+
  # add_pvalue(df_text2, label.size = 5)+
  # scale_x_discrete(labels = c("Vaccinated\nCOVID-19-_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19-_Male"=  "Vaccinated\ninfected",
  #                             "Vaccinated\nCOVID-19+_Female" = "Vaccinated\nuninfected",
  #                             "Vaccinated\nCOVID-19+_Male" = "Vaccinated\nuninfected"))+
  #stat_compare_means(label = "p.signif", size = 4,  label.x = 1.5)+
  ylim(3.3,6.3)+
  xlab(" ")+
  ylab("FcR3B_WT.S log10 MFI")
# + 
#   facet_wrap(.~assay, ncol = 4, scales = "free_x")
dev.off()
