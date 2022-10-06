
library(systemsseRology)
library(ggplot2)
library(DMwR)
library(patchwork)

setwd("/Users/dz480/Documents/projects/covid/az")

load("./data/sero.RData")
sero.vac <- droplevels(subset(sero, sero$group %in% levels(sero$group)[3:4]))

long.vac <- droplevels(subset(long, long$group %in% levels(long$group)[3:4]))
median.df <- long.vac %>%
  group_by(feature, assay, antigen, type, group) %>%
  summarise(value = median(value, na.rm = TRUE))
median.df <- as.data.frame(median.df)

diff.df <- spread(median.df, "group", "value")
diff.df$med_diff <- diff.df$`Vaccinated
COVID-19+` - diff.df$`Vaccinated
COVID-19-`




col.g1 <- "#ffeda0"
col.g2 <- "#f03b20"
col.g3 <- "#beaed4"
col.g4 <- "#fdc086"

g.arms <- levels(sero$group)[3:4]
col.arms <- c(col.g3, col.g4)

df <- sero.vac

my_colors <- list(group = c("Vaccinated\nCOVID-19-" = col.g3,
                            "Vaccinated\nCOVID-19+" = col.g4))


data <- df

X <- data[,58:ncol(data)]

# X <- X[, which(!grepl("Ebola|CMV|RSV|Influenza", colnames(X)))]
# X <- X[,-c(6,24,42,60,78,96,114,132)]

X <- X[, which(!grepl("RSV|CMV|Ebola|Influenza|WT.N|S1|S2", colnames(X)))]
X <- X[, which(grepl("WT.S|B.1.351.S|WT.RBD|B.1.351.RBD", colnames(X)))]


log.ind <- which(grepl("ADCD|IgG|IgA|IgM|FcR", colnames(X)))
X[,log.ind] <- log10(X[,log.ind])

X <- scale(X, center = TRUE, scale = TRUE)
#X <- knnImputation(X)
y <- factor(data$group)

df_features <- data.frame(name = colnames(X),
                          label = colnames(X))
opts_plot <- list(LV_ind = c(1,2), # which LVs to plot
                  colors = my_colors,
                  y_name = "group",
                  level = 0.75) 

# Feature selection
opts_sel <- list(threshold = 0.9, n_trials = 100, return_count = TRUE)
out <- select_repeat(X, y, selector = select_lasso, options = opts_sel)

df_count <- data.frame(features = names(out$feature_count), 
                       name = names(out$feature_count), 
                       selected = out$feature_count*100/opts_sel$n_trials,
                       mark = NA)
df_count <- df_count[which(df_count$selected > 0),]
df_count <- df_count[order(-df_count$selected),]
df_count$features <- df_features$label[match(df_count$features, df_features$name)]
df_count$features <- factor(df_count$features, levels = df_count$features)
# annotation where feature is enriched
for (ind_feat in 1:nrow(df_count)) {
  tmp_mean <- rep(NA, length = nlevels(y))
  for (ind_class in 1:nlevels(y)) {
    tmp_mean[ind_class] <- mean(X[which(y == levels(y)[ind_class]),
                                  which(colnames(X) == df_count$name[ind_feat])])
  }
  df_count$mark[ind_feat] <- levels(y)[which.max(tmp_mean)]
}
df_count$mark  <- factor(df_count$mark, levels = levels(y))

plt_bar <- ggplot(data = df_count, aes(x = features, y = selected, fill = mark)) +
  scale_fill_manual(values = my_colors$group) +
  geom_bar(stat = "identity", color = "black") +
  xlab("") + geom_hline(yintercept = opts_sel$threshold*100) +
  ylab("selected (%)") +
  labs(fill = "enriched in") +
  theme_classic() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 10, angle = 80, hjust = 1, vjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
print(plt_bar)

pdf("feature_selection_repetition_figure3.pdf", width = 10, height = 4)
print(plt_bar)
dev.off()

set.seed(123)
opts_sel <- list(threshold = 0.90, n_trials = 100, return_count = FALSE)
sel_features <- select_repeat(X, y, selector = select_lasso, options = opts_sel)


select <- function(X, y) { return(select_repeat(X, y, selector = select_lasso, options = opts_sel)) }

method = list(select = select, 
              train = train_ropls,
              predict = predict_ropls,
              score = score_accuracy)


opts_val <- list(n_folds = 5, rf_trials = 1, pt_trials = 1) # repeat with higher numbers
vals_Figure3_v2 <- validate_repeat(X, y, method, opts_val, n_trials = 100)
save(vals_Figure3_v2, file = "vals_Figure3_v2.RData")

cv.acc <- 0
for (i in 1:100){
  cv.acc <- cv.acc + vals_Figure3_v2[[i]]$cv_score
}
cv.acc <- cv.acc/100  ## 95.6% cv accuracy
#

plt_val <- visualize_validate(vals_Figure3_v2, options = list(y_label = "accuracy"))


pdf("vals_Figure3_val.pdf", width = 3, height = 3)
plt_val
dev.off()



# # 
# plt_val <- visualize_validate(vals, options = list(y_label = "accuracy"))
# 
# pdf("val.pdf", width = 3, height = 3)
# plt_val
# dev.off()

X_sel <- X[, sel_features]
opts_model <- list(n_LV = 2)
model <- train_ropls(X_sel, y, options = opts_model)

p.scores <- visualize_ropls_scores(model, y, options = opts_plot) +
  xlab("Scores on LV1") + ylab("Scores on LV2")+
  theme(legend.title = element_blank(), legend.position = "top",
        legend.text = element_text( size = 8))#+guides(fill=guide_legend(nrow=2,byrow=TRUE))



loadsm <- as.data.frame(model@loadingMN)
vips <- as.data.frame(model@vipVn)

loads <- cbind(loadsm, vips)
loads <- loads[order(loads$`model@vipVn`),]
loads$features <- row.names(loads)
loads$features <- factor(loads$features , levels = loads$features)

loads$group <- levels(y)[2]
loads$group[which(loads$p1 < 0)] <- levels(y)[1]

p.loads <- ggplot(data=loads, aes(x=features, y=p1))+
  #geom_bar(stat="identity", colour = "black") +coord_flip()+
  geom_bar(stat="identity", aes(fill = group), colour = "black") +coord_flip()+
  scale_fill_manual(breaks = g.arms , values = col.arms)+
  theme_classic() + ylab("LV1 loadings") + xlab("") +theme(legend.position = "none")


pdf("Figure3AB_plsda_lasso90.pdf", width = 7, height = 5)
p.scores  + p.loads 
#+ plt_val
dev.off()


# p.loads2 <- ggplot(data=loads, aes(x=features, y=p2))+
#   #geom_bar(stat="identity", colour = "black") +coord_flip()+
#   geom_bar(stat="identity", 
#            #aes(fill = "slategray"), 
#            fill = "slategray",
#            colour = "black") +coord_flip()+
#   scale_fill_manual(breaks = g.arms , values = col.arms)+
#   theme_classic() + ylab("LV2 loadings") + xlab("") +theme(legend.position = "none")
# 
# 
# pdf("plsda_lasso90_wtrbd_2v2.pdf", width = 7, height = 5)
# p.scores  + p.loads2 
# #+ plt_val
# dev.off()


df.vip <- data.frame(features = names(model@vipVn), VIP = model@vipVn)

median.sel <- subset(diff.df, diff.df$feature %in% df.vip$features)
median.sel$enrichment <- "Vaccinated\nCOVID-19+"
median.sel$enrichment[which(median.sel$med_diff < 0)] <- "Vaccinated\nCOVID-19-"

df.vip <- merge(df.vip, median.sel, by.x = "features", by.y = "feature")
df.vip <- df.vip[order(df.vip$VIP),]
df.vip$features <- factor(df.vip$features , levels = df.vip$features)


p.vip <- ggplot(data=df.vip, aes(x=features, y=VIP))+
  #geom_bar(stat="identity", colour = "black") +coord_flip()+
  geom_bar(stat="identity", aes(fill = enrichment), colour = "black") +coord_flip()+
  scale_fill_manual(breaks = g.arms , values = col.arms)+
  theme_classic() + ylab("VIP") + xlab("") +theme(legend.position = "none")




pdf("Figure3AB_plsda_lasso90_wtrbd_vip.pdf", width = 7, height = 5)
p.scores  + p.vip
#+ plt_val
dev.off()




library(DMwR)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(patchwork)
library(paletteer)
library(circlize)
library(RColorBrewer)
library(Hmisc)
library(tidyverse)
library(ComplexHeatmap)
library(network)
library(stringr)



pred <- X

res <- rcorr(as.matrix(pred), type = "spearman")

resr <- as.data.frame(res$r)
resr$selfeat <- row.names(resr)
df.r <- subset(resr, resr$selfeat %in% sel_features)

resp <- as.data.frame(res$P)
resp$selfeat <- row.names(resp)
df.p <- subset(resp, resr$selfeat %in% sel_features)


dfr <- df.r %>% gather(feature, rho, -c(selfeat))
dfp <- df.p %>% gather(feature, p, -c(selfeat))
dfp$q <- p.adjust(dfp$p, method = "BH")


identical(dfr$selfeat, dfp$selfeat)
identical(dfr$feature, dfp$feature)

df <- cbind(dfr, dfp[,c("p","q")])
df$antigen <- as.vector(str_split_fixed(df$feature, "_", 2)[,2])

df.sub <- subset(df, df$q < 0.01 & abs(df$rho) > 0.8)
write.csv(df.sub, file = "df.sub.csv", row.names = FALSE)

link.widths <- abs(df.sub$rho)*5
#link.widths <- ifelse(abs(df.sub$rho) > 0.5, 5, 1)

link.clr <- ifelse(df.sub$rho > 0, "#e7d4e8", "#c7eae5")

#col_fun = colorRamp2(c(-1, -0.5, 0, 0.5, 1), c("#003c30","#c7eae5", "white", "#e7d4e8", "#40004b"))


net <- network(df.sub, matrix.type = "edgelist",multiple = TRUE, ignore.eval = FALSE, directed = FALSE)

node.cols <- network.vertex.names(net)
fc.ind <- which(grepl("FcR", node.cols))
node.cols[fc.ind] <- "#7fc97f"
titer.ind <- which(grepl("Ig", node.cols))
#node.cols[titer.ind] <- "#beaed4"
node.cols[titer.ind] <- "yellow"
func.ind <- which(grepl("ADCD|ADCP|ADNP|NK|MIP|IFN|CD107", node.cols))
node.cols[func.ind] <- "#fdc086"

pdf(paste0("Figure3D_cornet.pdf"), width = 14, height = 14)
print(plot.network(net, label = network.vertex.names(net), 
                   #vertex.col = "slategray",
                   vertex.col = as.color(node.cols),
                   #mode = "circle",
                   edge.col = link.clr,
                   edge.lwd = link.widths,
                   pad = 2,
                   label.cex = 0.8))
dev.off()



# lasso selected univariates ----------------------------------------------

temp <- subset(long, long$group %in% levels(long$group)[3:4])

stat.wil <- temp %>%
  group_by(feature, assay, antigen, type) %>%
  wilcox_test(value ~ group)

stat.wil$p.star <- symnum(stat.wil$p, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))
stat.wil$q <- p.adjust(stat.wil$p, method = "BH")
stat.wil$q.star <- symnum(stat.wil$q, corr = FALSE, na = FALSE, 
                          cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                          symbols = c("***", "**", "*", " "))


temp <- subset(long, long$feature %in% sel_features & long$group %in% levels(long$group)[3:4])

df_text <- subset(stat.wil, stat.wil$feature %in% sel_features & stat.wil$q < 0.05)
df_text$y.position <- 6.1

df_text2 <- df_text[, c(1:7, 14, 15)]
colnames(df_text2)[8] <- "p.signif"

levels(loads$features)

temp$feature <- factor(temp$feature, levels = rev(levels(loads$features)))
temp <- temp[order(temp$feature),]

df_text2$feature <- factor(df_text2$feature, levels = rev(levels(loads$features)))

pdf("Figure3C_lasso_uni_srbd.pdf", width = 10, height = 10)
ggplot(temp, aes(group, log10(value)))+
  geom_violin(scale = "width", aes(fill = group),alpha = 0.8, draw_quantiles = c(0.5))+ 
  geom_jitter(height = 0, width = 0.1, size = 0.5 , color = "gray32")+ 
  scale_fill_manual(values = col.arms, breaks = g.arms) +
  theme_classic() + theme(axis.text.x = element_text( color = "black",size = 8),
                          legend.position =  "bottom", aspect.ratio = 1) + 
  # add_pvalue(df_text2, step.increase = 0.08, step.group.by = "assay", tip.length = 0,
  #            colour = "vaccine")+
  add_pvalue(df_text2, label.size = 5)+
  scale_x_discrete(labels = c("Vaccine_control" = "Vaccinated\nuninfected",
                              "Vaccine_case"=  "Vaccinated\ninfected"))+
  #stat_compare_means(label = "p.signif", size = 4,  label.x = 1.5)+
  ylim(3.3,6.3)+
  xlab(" ") +ylab("log10 MFI")+ facet_wrap(.~feature, ncol = 5, scales = "free_x")
dev.off()



