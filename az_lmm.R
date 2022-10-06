setwd("/Users/dz480/Documents/projects/covid/az/")
load("sero.RData")

library(systemsseRology)
library(ggplot2)
library(DMwR)
library(patchwork)

data <- sero.vac

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




library(lme4)
library(effectsize)

data.mixed.linear.model <- cbind(data[,1:57], X)

#   
feature.name <- colnames(X)
pvalue <- c()
coefficient <- data.frame()
t.value <- c()
beta <- c()
effect.size <- data.frame()

#Per measurement
for(feat.index in 1:length(feature.name)){
  #prepare the dataset.
  response <- feature.name[feat.index]
  print(response)
  temp.data.mlm <- data.mixed.linear.model[, c(1:57, feat.index + 57)]
  colnames(temp.data.mlm)[58] <- 'response'
  
  #Basic model without considering Pat.clinical.label variable.
  #model.null <- lmer(response ~ 1   + sc_gender + age_calc+ bmi + (1|pid) , data = temp.data.mlm, REML = FALSE)
  model.null <- lm(response ~ 1   + sc_gender + age_calc+ bmi + social_smoking + social_alcohol + sc_healthcare+
                     sc_race, data = temp.data.mlm)
  #saveRDS(summary(model.null), file.path(saved.dic, paste(response, 'null.model.RDS')))
  #saveRDS(coef(model.null)$Pat.ID, file.path(saved.dic, paste(response, 'null.mode.coeff.RDS')))
  
  #Alternative model with Pat.clinical.label variable.
  #model.alternative <- lmer(response ~ 1 + sc_gender + age_calc + bmi + group + (1|pid) , data = temp.data.mlm, REML = FALSE)
  model.alternative <- lm(response ~ 1 + sc_gender + age_calc + bmi + social_smoking + social_alcohol + sc_healthcare+
                            sc_race+ group , data = temp.data.mlm)
  #saveRDS(summary(model.alternative), file.path(saved.dic, paste(response, 'alternative.model.RDS')))
  #saveRDS(coef(model.alternative)$Pat.ID, file.path(saved.dic, paste(response, 'alternative.mode.coeff.RDS')))
  
  #model.alternative Summary
  summary.model.alternative <- summary(model.alternative)
  coefficient <- rbind(coefficient, summary.model.alternative$coefficients[, 1])
  beta <- c(beta, summary.model.alternative$coefficients["groupVaccinated\nCOVID-19+", 1])
  #t.value <- rbind(t.value, summary.model.alternative$coefficients[, 3])
  t.value <- c(t.value, summary.model.alternative$coefficients["groupVaccinated\nCOVID-19+", 3])
  
  #Calculate the effect size on full model
  cal.effect.clinical.label <- effectsize::effectsize(model.alternative)
  #saveRDS(cal.effect.clinical.label, file.path(saved.dic, paste(response, 'alternative.model.effectsize.RDS')))
  effect.size <- rbind(effect.size, cal.effect.clinical.label$Std_Coefficient)
  
  #Likelihood ratio test
  LRT <- anova(model.alternative, model.null)
  #saveRDS(LRT, file.path(saved.dic, paste(response, 'LRT.RDS')))
  #pvalue <- c(pvalue, LRT$`Pr(>Chisq)`[2])
  pvalue <- c(pvalue, LRT$`Pr(>F)`[2])
}


# library(lmtest)
# LRT <- lrtest(model.alternative, model.null)

labels <- feature.name
#labels[7:48] <- paste0("Ig", labels[7:48])

covid <- rep("Vaccinated\nCOVID-19+", length(t.value))
covid[which(t.value < 0)] <- "Vaccinated\nCOVID-19-"

# type <- c(rep("Function",6), rep("IgG1", 7), rep("IgG2", 7), rep("IgG3", 7), rep("IgG4", 7),
#           rep("IgA", 7), rep("IgM", 7), rep("FcR", 28))
# 
# assays <- c("Function","IgG1","IgG2","IgG3","IgG4",
#             "IgA", "IgM","FcR")
# col.assays <- c("pink", rgb(68,147,80,max=255),
#                 rgb(47,104,61,max=255),
#                 rgb(90,179,120,max=255),
#                 rgb(78,165,156,max=255),
#                 rgb(233,153,64,max=255),
#                 rgb(223,100,60,max=255),
#                 "slategray")

library(ggrepel)
# lmm.df <- data.frame(feature = feature.name, pvalue = pvalue, tvalue = t.value, coefficient = beta,
#                      lab = labels, assay = type, group = tb)
lmm.df <- data.frame(feature = feature.name, pvalue = pvalue, tvalue = t.value, coefficient = beta,
                     lab = labels,  group = covid)
#lmm.df$lab[which(lmm.df$pvalue > 0.05)] <- ""

col.g3 <- "#beaed4"
col.g4 <- "#fdc086"

g.arms <- c("Vaccinated\nCOVID-19-", "Vaccinated\nCOVID-19+" )
col.arms <- c(col.g3, col.g4)

p <- ggplot(lmm.df,
            aes(coefficient, -log10(pvalue), label = lab)) +
  geom_hline(yintercept=-log10(0.05), colour = "slategray", linetype="dashed")+
  #geom_vline(xintercept =0, colour = "slategray", linetype="dashed")+
  scale_fill_manual(breaks = g.arms, values = col.arms)+
  geom_text_repel(size = 2,  show.legend = FALSE) + 
  #scale_color_manual(breaks = assays, values = col.assays)+
  geom_point(size = 1.5, aes(fill = covid), shape = 21) + theme_classic()+
  #xlab("Coefficient (No Recurrent TB vs Recurrent TB)")+
  theme(legend.title = element_blank(), legend.position = "top")
p

pdf("lm_az.pdf", height = 5, width = 7)
p
dev.off()


