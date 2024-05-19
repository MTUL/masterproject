setwd("/home/projects/cu_10184/data/projects/groenbaek/data/data_raw/EVI2_EPICarray/")
library(tidyverse)
library(xtable)
library(ggpubr)
library(lme4)
library(LabApplStat)
library(lattice)

theme_set(theme_bw())

# Functions to be used:
inverse_logit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

# Loading dataset
load("trial_data_720sites.Rdata")

#### Table 1: ####
table1data <- dlong %>% select(patient_id, treatment_arm, TET2_mut) %>% 
      distinct()

table1data %>% group_by(treatment_arm) %>% tally() %>% xtable()

xtable(addmargins(table(table1data$treatment_arm, table1data$TET2_mut)))

#### Descriptives ####

## Descriptive Figure 1:
# Panel A: beta value distributions split by mutation
beta_hist <- dlong %>% 
      filter(visit == 1) %>% 
      ggplot(aes(x = beta, fill = TET2_mut)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 100) +
      facet_wrap(~ TET2_mut) +
      theme(legend.position = "none")

# Panel B: Four representative CpG sites 
foursites <- c("cg12275410", "cg01498832", "cg13747899", "cg00782708")
foursites_hist <- dlong %>% 
      filter(probeID %in% foursites) %>% 
      ggplot(aes(x = beta, fill = TET2_mut, group = TET2_mut)) +
      geom_histogram(position = "dodge", bins = 20) +
      facet_wrap(~ probeID) +
      theme(legend.position = "none")

# Panel C: Distributions of pairwise Pearson correlation coefficients between baseline beta values in TET2 mut and WT subsets
beta_df.base.mut <- dlong %>% 
      filter(visit == 1 & TET2_mut == "mut") %>% 
      select(arrayID, probeID, beta) %>% 
      pivot_wider(values_from = beta, names_from = probeID) %>% 
      column_to_rownames("arrayID")

beta_df.base.wt <- dlong %>% 
      filter(visit == 1 & TET2_mut == "WT") %>% 
      select(arrayID, probeID, beta) %>% 
      pivot_wider(values_from = beta, names_from = probeID) %>% 
      column_to_rownames("arrayID")

cor.base.mut <- cor(beta_df.base.mut)
cor.base.mut[!upper.tri(cor.base.mut)] <- NA
cor.base.mut.df <- data.frame(cor = as.vector(cor.base.mut),
                              TET2_mut = "mut") %>% 
      filter(!is.na(cor))

cor.base.wt <- cor(beta_df.base.wt)
cor.base.wt[!upper.tri(cor.base.wt)] <- NA
cor.base.wt.df <- data.frame(cor = as.vector(cor.base.wt),
                             TET2_mut = "WT") %>% 
      filter(!is.na(cor))

cor_hist <- bind_rows(cor.base.mut.df, cor.base.wt.df) %>% 
      mutate(TET2_mut = factor(TET2_mut, levels = c("WT", "mut"))) %>% 
      ggplot(aes(x = cor, fill = TET2_mut)) +
      geom_histogram(bins = 100, aes(y = after_stat(density))) +
      facet_wrap( ~ TET2_mut) +
      theme(legend.position = "none") +
      xlim(c(-1,1))

# Panel D: Distributions of pairwise Pearson correlation coefficients between changes in beta values from baseline to end of study
#  in TET2 mut and WT subsets
beta_df.end.mut <- dlong %>% 
      filter(visit == 5 & TET2_mut == "mut") %>% 
      select(arrayID, probeID, beta) %>% 
      pivot_wider(values_from = beta, names_from = probeID) %>% 
      column_to_rownames("arrayID")

beta_df.end.wt <- dlong %>% 
      filter(visit == 5 & TET2_mut == "WT") %>% 
      select(arrayID, probeID, beta) %>% 
      pivot_wider(values_from = beta, names_from = probeID) %>% 
      column_to_rownames("arrayID")

z_df.mut <- beta_df.end.mut - beta_df.base.mut
z_df.wt <- beta_df.end.wt - beta_df.base.wt

cor.z.mut <- cor(z_df.mut)
cor.z.mut[!upper.tri(cor.z.mut)] <- NA
cor.z.mut.df <- data.frame(cor = as.vector(cor.z.mut),
                              TET2_mut = "mut") %>% 
      filter(!is.na(cor))

cor.z.wt <- cor(z_df.wt)
cor.z.wt[!upper.tri(cor.z.wt)] <- NA
cor.z.wt.df <- data.frame(cor = as.vector(cor.z.wt),
                             TET2_mut = "WT") %>% 
      filter(!is.na(cor))

cor.z_hist <- bind_rows(cor.z.mut.df, cor.z.wt.df) %>% 
      mutate(TET2_mut = factor(TET2_mut, levels = c("WT", "mut"))) %>% 
      ggplot(aes(x = cor, fill = TET2_mut)) +
      geom_histogram(bins = 100, aes(y = after_stat(density))) +
      facet_wrap( ~ TET2_mut) +
      theme(legend.position = "none") +
      xlim(c(-1,1))

pdf("figures/Descriptives1.pdf")
ggarrange(beta_hist, foursites_hist, cor_hist, cor.z_hist, 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()

# Extra: Correlations between individuals instead of CpG sites:
t.beta_df.base.mut <- t(beta_df.base.mut)
t.beta_df.base.wt <- t(beta_df.base.wt)
cor.base.mut <- cor(t.beta_df.base.mut)
cor.base.mut[!upper.tri(cor.base.mut)] <- NA
cor.base.mut.df <- data.frame(cor = as.vector(cor.base.mut),
                              TET2_mut = "mut") %>% 
      filter(!is.na(cor))

cor.base.wt <- cor(t.beta_df.base.wt)
cor.base.wt[!upper.tri(cor.base.wt)] <- NA
cor.base.wt.df <- data.frame(cor = as.vector(cor.base.wt),
                             TET2_mut = "WT") %>% 
      filter(!is.na(cor))

cor_hist.t <- bind_rows(cor.base.mut.df, cor.base.wt.df) %>% 
      mutate(TET2_mut = factor(TET2_mut, levels = c("WT", "mut"))) %>% 
      ggplot(aes(x = cor, fill = TET2_mut)) +
      geom_histogram(bins = 100, aes(y = after_stat(density))) +
      facet_wrap( ~ TET2_mut) +
      theme(legend.position = "none") +
      xlim(c(-1,1))
# Maybe show this figure as well.

#### Factor diagram ####

## For some reason the DD function fails when given the probe x mut product factor
# Diagram 1: No product factors
t0 <- Sys.time()
factordiagram1 <- plot(LabApplStat::DD(Mval ~ TET2_mut + visit + treatment_constrained, 
                                      random =  ~ probeID + patient_id, 
                                      data = dlong))
t1 <- Sys.time()
(time_spent1 <- difftime(t1, t0, units = "mins"))


# Diagram 2: Including the TET2 x visit product factor
t0 <- Sys.time()
factordiagram2 <- plot(LabApplStat::DD(Mval ~ TET2_mut * visit + treatment_constrained, 
                                       random =  ~ probeID + patient_id, 
                                       data = dlong))
t1 <- Sys.time()
(time_spent2 <- difftime(t1, t0, units = "mins"))
factordiagram2

svg("figures/factordiagram1.svg", height = 6, width = 8)
factordiagram1
dev.off()

pdf("figures/factordiagram1.pdf", height = 6, width = 8)
factordiagram1
dev.off()

svg("figures/factordiagram2.svg", height = 6, width = 8)
factordiagram2
dev.off()

pdf("figures/factordiagram2.pdf", height = 6, width = 8)
factordiagram2
dev.off()

#### Simple fixed-effects models for one CpG site at a time #### 
dlong1 <- dlong %>% filter(probeID == foursites[1])
fixmod1 <- lm(Mval ~ TET2_mut * visit + treatment_constrained, data = dlong1)
dlong2 <- dlong %>% filter(probeID == foursites[2])
fixmod2 <- lm(Mval ~ TET2_mut * visit + treatment_constrained, data = dlong2)
dlong3 <- dlong %>% filter(probeID == foursites[3])
fixmod3 <- lm(Mval ~ TET2_mut * visit + treatment_constrained, data = dlong3)
dlong4 <- dlong %>% filter(probeID == foursites[4])
fixmod4 <- lm(Mval ~ TET2_mut * visit + treatment_constrained, data = dlong4)

newdata.fixed <- dlong1 %>% select(TET2_mut, visit, treatment_arm, treatment_constrained) %>% distinct()
newdata.fixed$cg12275410 <- inverse_logit(predict(fixmod1, newdata = newdata.fixed))
newdata.fixed$cg01498832 <- inverse_logit(predict(fixmod2, newdata = newdata.fixed))
newdata.fixed$cg13747899 <- inverse_logit(predict(fixmod3, newdata = newdata.fixed))
newdata.fixed$cg00782708 <- inverse_logit(predict(fixmod4, newdata = newdata.fixed))

pdf("figures/fixed_effects.onesite.pdf", width = 6, height = 4) 
newdata.fixed %>% 
      select(-treatment_constrained) %>%
      mutate(visit = if_else(visit == "1", "baseline", "end")) %>% 
      ggplot(aes(x = visit, y = cg13747899, 
                 group = interaction(TET2_mut, treatment_arm),
                 color = TET2_mut,
                 linetype = treatment_arm)) +
      geom_line() +
      ylim(c(0,0.7))
dev.off()


pdf("figures/fixed_effects.foursites.pdf")
newdata.fixed %>% 
      select(-treatment_constrained) %>% 
      pivot_longer(cols = c(-TET2_mut, -visit, -treatment_arm)) %>% 
      ggplot(aes(x = visit, y = value, 
                 group = interaction(TET2_mut, treatment_arm),
                 color = TET2_mut,
                 linetype = treatment_arm)) +
      geom_line() +
      facet_wrap(~name) +
      ylim(c(0,1))
dev.off()

#### Attempts at simulation of full dataset: Selection of random effects ####


## Models under the null hypothesis of no treatment effect
# mnull1 <- lmer(Mval ~ factor(TET2_mut) * visit + (1|patient_id) + (1|probeID), data = dlong)
# mnull2 <- lmer(Mval ~ factor(TET2_mut) * visit + (1|patient_visit) + (1|probeID), data = dlong)
# mnull3 <- lmer(Mval ~ factor(TET2_mut) * visit + (1|patient_id) + (1|probe_mut), data = dlong)
# mnull4 <- lmer(Mval ~ factor(TET2_mut) * visit + (1|patient_visit) + (1|probe_mut), data = dlong)
# mnull5 <- lmer(Mval ~ factor(TET2_mut) * visit + (1|patient_id) + (1|patient_visit) + (1|probeID) + (1|probe_mut), data = dlong)
# save(mnull1, mnull2, mnull3, mnull4, mnull5, file = "null_models.Rdata")

load("null_models.Rdata")

mtreat1 <- lmer(Mval ~ treatment_constrained + factor(TET2_mut) * visit + (1|patient_id) + (1|probeID), data = dlong)
summary(mtreat1)
confint(mtreat1, parm = "treatment_constrainedvitaminC")

mtreat4 <- lmer(Mval ~ treatment_constrained + factor(TET2_mut) * visit + (1|patient_visit) + (1|probe_mut), data = dlong)
summary(mtreat4)
confint(mtreat4, parm = "treatment_constrainedvitaminC")
# Interesting difference between mtreat1 and mtreat4!!

# Function to create similar plots to the descriptives but with simulated datasets
modelplot <- function(modelobject){
      # modelobject = mnull1
      dsim <- dlong %>% 
            mutate(Mval = simulate(modelobject)$sim_1,
                   beta = inverse_logit(Mval))
      beta_hist <- dsim %>% 
            filter(visit == 1) %>% 
            ggplot(aes(x = beta, fill = TET2_mut)) +
            geom_histogram(aes(y = after_stat(density)),
                           bins = 100) +
            facet_wrap(~ TET2_mut) +
            theme(legend.position = "none")
      
      # Panel B: Four representative CpG sites 
      foursites <- c("cg12275410", "cg01498832", "cg13747899", "cg00782708")
      foursites_hist <- dsim %>% 
            filter(probeID %in% foursites) %>% 
            ggplot(aes(x = beta, fill = TET2_mut, group = TET2_mut)) +
            geom_histogram(position = "dodge", bins = 20) +
            facet_wrap(~ probeID) +
            theme(legend.position = "none")
      
      # Panel C: Distributions of pairwise Pearson correlation coefficients between baseline beta values in TET2 mut and WT subsets
      beta_df.base.mut <- dsim %>% 
            filter(visit == 1 & TET2_mut == "mut") %>% 
            select(arrayID, probeID, beta) %>% 
            pivot_wider(values_from = beta, names_from = probeID) %>% 
            column_to_rownames("arrayID")
      
      beta_df.base.wt <- dsim %>% 
            filter(visit == 1 & TET2_mut == "WT") %>% 
            select(arrayID, probeID, beta) %>% 
            pivot_wider(values_from = beta, names_from = probeID) %>% 
            column_to_rownames("arrayID")
      
      cor.base.mut <- cor(beta_df.base.mut)
      cor.base.mut[!upper.tri(cor.base.mut)] <- NA
      cor.base.mut.df <- data.frame(cor = as.vector(cor.base.mut),
                                    TET2_mut = "mut") %>% 
            filter(!is.na(cor))
      
      cor.base.wt <- cor(beta_df.base.wt)
      cor.base.wt[!upper.tri(cor.base.wt)] <- NA
      cor.base.wt.df <- data.frame(cor = as.vector(cor.base.wt),
                                   TET2_mut = "WT") %>% 
            filter(!is.na(cor))
      
      cor_hist <- bind_rows(cor.base.mut.df, cor.base.wt.df) %>% 
            mutate(TET2_mut = factor(TET2_mut, levels = c("WT", "mut"))) %>% 
            ggplot(aes(x = cor, fill = TET2_mut)) +
            geom_histogram(bins = 100, aes(y = after_stat(density))) +
            facet_wrap( ~ TET2_mut) +
            theme(legend.position = "none") +
            xlim(c(-1,1))
      
      # Panel D: Distributions of pairwise Pearson correlation coefficients between changes in beta values from baseline to end of study
      #  in TET2 mut and WT subsets
      beta_df.end.mut <- dsim %>% 
            filter(visit == 5 & TET2_mut == "mut") %>% 
            select(arrayID, probeID, beta) %>% 
            pivot_wider(values_from = beta, names_from = probeID) %>% 
            column_to_rownames("arrayID")
      
      beta_df.end.wt <- dsim %>% 
            filter(visit == 5 & TET2_mut == "WT") %>% 
            select(arrayID, probeID, beta) %>% 
            pivot_wider(values_from = beta, names_from = probeID) %>% 
            column_to_rownames("arrayID")
      
      z_df.mut <- beta_df.end.mut - beta_df.base.mut
      z_df.wt <- beta_df.end.wt - beta_df.base.wt
      
      cor.z.mut <- cor(z_df.mut)
      cor.z.mut[!upper.tri(cor.z.mut)] <- NA
      cor.z.mut.df <- data.frame(cor = as.vector(cor.z.mut),
                                 TET2_mut = "mut") %>% 
            filter(!is.na(cor))
      
      cor.z.wt <- cor(z_df.wt)
      cor.z.wt[!upper.tri(cor.z.wt)] <- NA
      cor.z.wt.df <- data.frame(cor = as.vector(cor.z.wt),
                                TET2_mut = "WT") %>% 
            filter(!is.na(cor))
      
      cor.z_hist <- bind_rows(cor.z.mut.df, cor.z.wt.df) %>% 
            mutate(TET2_mut = factor(TET2_mut, levels = c("WT", "mut"))) %>% 
            ggplot(aes(x = cor, fill = TET2_mut)) +
            geom_histogram(bins = 100, aes(y = after_stat(density))) +
            facet_wrap( ~ TET2_mut) +
            theme(legend.position = "none") +
            xlim(c(-1,1))
      t.beta_df.base.mut <- t(beta_df.base.mut)
      t.beta_df.base.wt <- t(beta_df.base.wt)
      cor.base.mut <- cor(t.beta_df.base.mut)
      cor.base.mut[!upper.tri(cor.base.mut)] <- NA
      cor.base.mut.df <- data.frame(cor = as.vector(cor.base.mut),
                                    TET2_mut = "mut") %>% 
            filter(!is.na(cor))
      
      cor.base.wt <- cor(t.beta_df.base.wt)
      cor.base.wt[!upper.tri(cor.base.wt)] <- NA
      cor.base.wt.df <- data.frame(cor = as.vector(cor.base.wt),
                                   TET2_mut = "WT") %>% 
            filter(!is.na(cor))
      
      cor_hist.t <- bind_rows(cor.base.mut.df, cor.base.wt.df) %>% 
            mutate(TET2_mut = factor(TET2_mut, levels = c("WT", "mut"))) %>% 
            ggplot(aes(x = cor, fill = TET2_mut)) +
            geom_histogram(bins = 100, aes(y = after_stat(density))) +
            facet_wrap( ~ TET2_mut) +
            theme(legend.position = "none") +
            xlim(c(-1,1))
      
      return(list(beta_hist, foursites_hist, cor_hist, cor.z_hist, cor_hist.t))
}

plot <- modelplot(mnull1)
pdf(paste0("figures/", "nullmodel1.simulationplot", ".pdf"))
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()
plot[[5]]

plot <- modelplot(mtreat1)
pdf(paste0("figures/", "treatmodel1.simulationplot", ".pdf"))
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()
plot[[5]]

plot <- modelplot(mnull2)
pdf(paste0("figures/", "nullmodel2.simulationplot", ".pdf"))
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()

plot <- modelplot(mnull3)
pdf(paste0("figures/", "nullmodel3.simulationplot", ".pdf"))
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()

plot <- modelplot(mnull4)
pdf(paste0("figures/", "nullmodel4.simulationplot", ".pdf"))
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()

set.seed(9)
plot <- modelplot(mtreat4)
pdf(paste0("figures/", "treatmodel4.simulationplot", ".pdf"))
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()


mtreat4 <- lmer(Mval ~ treatment_constrained + factor(TET2_mut) * visit + (1|patient_visit) + (1|probe_mut), data = dlong)
plot <- modelplot(mtreat5)
pdf(paste0("figures/", "treatmodel5.simulationplot", ".pdf"))
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()


## Also, the point about mean differences between mutation groups in each CpG site:
meandiffs1 <- dlong %>% group_by(TET2_mut, probeID) %>% 
      summarise(mean = mean(beta)) %>% 
      ungroup() %>% group_by(probeID) %>% 
      summarise(meandiff = mean[TET2_mut == "mut"] - mean[TET2_mut == "WT"])
h1 <- meandiffs1 %>%  ggplot(aes(x = meandiff)) +
      geom_histogram()

meandiffs2 <- dlong %>% 
      mutate(Mval = simulate(mtreat1)$sim_1,
             beta = inverse_logit(Mval)) %>% 
      group_by(TET2_mut, probeID) %>% 
      summarise(mean = mean(beta)) %>% 
      ungroup() %>% group_by(probeID) %>% 
      summarise(meandiff = mean[TET2_mut == "mut"] - mean[TET2_mut == "WT"])
h2 <- meandiffs2 %>% ggplot(aes(x = meandiff)) +
      geom_histogram()
meandiffs3 <- dlong %>% 
      mutate(Mval = simulate(mtreat4)$sim_1,
             beta = inverse_logit(Mval)) %>% 
      group_by(TET2_mut, probeID) %>% 
      summarise(mean = mean(beta)) %>% 
      ungroup() %>% group_by(probeID) %>% 
      summarise(meandiff = mean[TET2_mut == "mut"] - mean[TET2_mut == "WT"])
h3 <- meandiffs3 %>% ggplot(aes(x = meandiff)) +
      geom_histogram()

ggarrange(h1, h2, h3)

#### Table: Coefficients and 95% CIs in three models ####
st1 <- summary(mtreat1)
st4 <- summary(mtreat4)
sn4 <- summary(mnull4)

# cit1 <- confint(mtreat1)
# cit4 <- confint(mtreat4)
# cin4 <- confint(mnull4)
save(cit1, cit4, cin4, file = "confidence_intervals.Rdata")
load("confidence_intervals.Rdata")

st1$sigma
print(VarCorr(mtreat4), comp = "Variance")
as.data.frame(VarCorr(mtreat1))[,5]
as.data.frame(VarCorr(mtreat1))[,1]

tabt1 <- data.frame(est_mtreat1 = st1$coefficients[,1]) %>% 
      rownames_to_column() %>% 
      bind_rows(data.frame(rowname = as.data.frame(VarCorr(mtreat1))[,1],
                           est_mtreat1 = as.data.frame(VarCorr(mtreat1))[,4]))
tabt1$rowname[tabt1$rowname == "probeID"] <- "tau2.sq"
tabt1$rowname[tabt1$rowname == "patient_id"] <- "tau1.sq"
tabt1$rowname[tabt1$rowname == "Residual"] <- "sigma.sq"
citabt1 <- cit1 %>% data.frame() %>% 
      rownames_to_column()
citabt1$X2.5..[citabt1$rowname == ".sig01"] <- citabt1$X2.5..[citabt1$rowname == ".sig01"]^2
citabt1$X97.5..[citabt1$rowname == ".sig01"] <- citabt1$X97.5..[citabt1$rowname == ".sig01"]^2
citabt1$rowname[citabt1$rowname == ".sig01"] <- "tau2.sq"
citabt1$X2.5..[citabt1$rowname == ".sig02"] <- citabt1$X2.5..[citabt1$rowname == ".sig02"]^2
citabt1$X97.5..[citabt1$rowname == ".sig02"] <- citabt1$X97.5..[citabt1$rowname == ".sig02"]^2
citabt1$rowname[citabt1$rowname == ".sig02"] <- "tau1.sq"
citabt1$X2.5..[citabt1$rowname == ".sigma"] <- citabt1$X2.5..[citabt1$rowname == ".sigma"]^2
citabt1$X97.5..[citabt1$rowname == ".sigma"] <- citabt1$X97.5..[citabt1$rowname == ".sigma"]^2
citabt1$rowname[citabt1$rowname == ".sigma"] <- "sigma.sq"
fullt1 <- tabt1 %>% left_join(citabt1)

tabt4 <- data.frame(est_mtreat4 = st4$coefficients[,1]) %>% 
      rownames_to_column() %>% 
      bind_rows(data.frame(rowname = as.data.frame(VarCorr(mtreat4))[,1],
                           est_mtreat4 = as.data.frame(VarCorr(mtreat4))[,4]))
tabt4$rowname[tabt4$rowname == "probe_mut"] <- "tau2.sq"
tabt4$rowname[tabt4$rowname == "patient_visit"] <- "tau1.sq"
tabt4$rowname[tabt4$rowname == "Residual"] <- "sigma.sq"
citabt4 <- cit4 %>% data.frame() %>% 
      rownames_to_column()
citabt4$X2.5..[citabt4$rowname == ".sig01"] <- citabt4$X2.5..[citabt4$rowname == ".sig01"]^2
citabt4$X97.5..[citabt4$rowname == ".sig01"] <- citabt4$X97.5..[citabt4$rowname == ".sig01"]^2
citabt4$rowname[citabt4$rowname == ".sig01"] <- "tau2.sq"
citabt4$X2.5..[citabt4$rowname == ".sig02"] <- citabt4$X2.5..[citabt4$rowname == ".sig02"]^2
citabt4$X97.5..[citabt4$rowname == ".sig02"] <- citabt4$X97.5..[citabt4$rowname == ".sig02"]^2
citabt4$rowname[citabt4$rowname == ".sig02"] <- "tau1.sq"
citabt4$X2.5..[citabt4$rowname == ".sigma"] <- citabt4$X2.5..[citabt4$rowname == ".sigma"]^2
citabt4$X97.5..[citabt4$rowname == ".sigma"] <- citabt4$X97.5..[citabt4$rowname == ".sigma"]^2
citabt4$rowname[citabt4$rowname == ".sigma"] <- "sigma.sq"
fullt4 <- tabt4 %>% left_join(citabt4)

tabn4 <- data.frame(est_mnull4 = sn4$coefficients[,1]) %>% 
      rownames_to_column() %>% 
      bind_rows(data.frame(rowname = as.data.frame(VarCorr(mnull4))[,1],
                           est_mnull4 = as.data.frame(VarCorr(mnull4))[,4]))
tabn4$rowname[tabn4$rowname == "probe_mut"] <- "tau2.sq"
tabn4$rowname[tabn4$rowname == "patient_visit"] <- "tau1.sq"
tabn4$rowname[tabn4$rowname == "Residual"] <- "sigma.sq"
citabn4 <- cin4 %>% data.frame() %>% 
      rownames_to_column()
citabn4$X2.5..[citabn4$rowname == ".sig01"] <- citabn4$X2.5..[citabn4$rowname == ".sig01"]^2
citabn4$X97.5..[citabn4$rowname == ".sig01"] <- citabn4$X97.5..[citabn4$rowname == ".sig01"]^2
citabn4$rowname[citabn4$rowname == ".sig01"] <- "tau2.sq"
citabn4$X2.5..[citabn4$rowname == ".sig02"] <- citabn4$X2.5..[citabn4$rowname == ".sig02"]^2
citabn4$X97.5..[citabn4$rowname == ".sig02"] <- citabn4$X97.5..[citabn4$rowname == ".sig02"]^2
citabn4$rowname[citabn4$rowname == ".sig02"] <- "tau1.sq"
citabn4$X2.5..[citabn4$rowname == ".sigma"] <- citabn4$X2.5..[citabn4$rowname == ".sigma"]^2
citabn4$X97.5..[citabn4$rowname == ".sigma"] <- citabn4$X97.5..[citabn4$rowname == ".sigma"]^2
citabn4$rowname[citabn4$rowname == ".sigma"] <- "sigma.sq"
fulln4 <- tabn4 %>% left_join(citabn4)

joint1 <- fullt1 %>% 
      transmute(rowname = rowname,
                est.ci.mtreat1 = paste0(round(est_mtreat1,3), " ",
                               "[", round(X2.5.., 3), ",",
                               round(X97.5.., 3), "]"))
joint4 <- fullt4 %>% 
      transmute(rowname = rowname,
                est.ci.mtreat4 = paste0(round(est_mtreat4,3), " ",
                                "[", round(X2.5.., 3), ",",
                                round(X97.5.., 3), "]"))
joinn4 <- fulln4 %>% 
      transmute(rowname = rowname,
                est.ci.null4 = paste0(round(est_mnull4,3), " ",
                                "[", round(X2.5.., 3), ",",
                                round(X97.5.., 3), "]"))
fulltable <- joint1 %>% left_join(joint4) %>% left_join(joinn4)
xtable(fulltable)

#### Model diagnostics ####
dlong$mnull4.resid <- summary(mnull4)$residuals
dlong$mnull4.fitted <- fitted(mnull4)

resid.plot <- dlong %>% 
      ggplot(aes(x = mnull4.fitted, y = mnull4.resid)) +
            geom_point(alpha = 0.1) +
            geom_smooth() +
      xlab("Fitted values") +
      ylab("Residuals")

qq <- dlong %>% ggplot(aes(sample = mnull4.resid)) +
      stat_qq() + 
      stat_qq_line() +
      xlab("Theoretical quantiles") +
      ylab("Sample quantiles")

ranef1 <- qqmath(ranef(mnull4)$probe_mut[,1], main="Random intercept by probe x mut")
ranef2 <- qqmath(ranef(mnull4)$patient_visit[,1], main="Random intercept by patient x visit")

pdf("figures/Diagnostics.pdf", width = 6, height = 3)
ggarrange(resid.plot, qq, 
          nrow = 1, ncol = 2, labels = c("A", "B"))
dev.off()

#### Visualization of simulated treatment effects ####
dlong$fitted_mnull4 <- getME(mnull4, name = "X") %*% fixef(mnull4) 

beta <- fixef(mnull4) 

small <- dlong %>% select(visit, fitted_mnull4, treatment_arm, treatment_constrained, TET2_mut) %>% 
      data.frame() %>% 
      distinct()
small.X <- model.matrix( ~ treatment_constrained + TET2_mut * visit, data = small)

possible_effects <- 1-exp(seq(from = 0, to = 1, length.out = 10))

possible_effects <- possible_effects[1:8]

small_predictions <- matrix(nrow = nrow(small), ncol = length(possible_effects))
for(i in 1:length(possible_effects)){
      beta_i <- c(beta[1], 0, beta[2:4])
      beta_i[2] <- possible_effects[i]
      pred_i <- small.X %*% beta_i
      small_predictions[, i] <- pred_i
}
colnames(small_predictions) <- paste("effect", 0:7)

pdf("figures/simulated_treatment_effects.pdf", width = 6, height = 3)
small_predictions %>% data.frame() %>% 
      bind_cols(small) %>% 
      pivot_longer(cols = starts_with("effect.")) %>% 
      mutate(fitted_beta = inverse_logit(value),
             effect = gsub("effect.", "Effect ", name)) %>% 
      filter(treatment_arm == "vitaminC" | name == "pred_0") %>%
      ggplot(aes(x = visit, y = fitted_beta, group = effect, color = effect)) +
      geom_line() +
      ylim(c(0,1)) +
      facet_wrap(~ TET2_mut)
dev.off()

# Numbers:
preddat <- small_predictions %>% data.frame() %>% 
      bind_cols(small) %>% 
      pivot_longer(cols = starts_with("effect.")) %>% 
      mutate(fitted_beta = inverse_logit(value),
             effect = gsub("effect.", "Effect ", name))
# Baseline means in mutation groups:
preddat %>% filter(visit  == 1) %>% distinct(TET2_mut, fitted_beta)
# Absolute difference: 0.329

# Reduction in end-of-study-methylation:
preddat %>% filter(visit == 5 & effect == "Effect 7") %>% 
      distinct(TET2_mut, treatment_arm, fitted_beta)
# In TET2 mutated arm, effect 8 reduces mean at end from 0.787 to 0.532 -> -0.255
## -> 68% of the entire mut vs WT difference
# In WT arm, effect 8 reduces mean at end from 0.443 to 0.197 -> -0.246
## -> 75%

preddat %>% filter(visit == 5 & effect == "Effect 1") %>% 
      distinct(TET2_mut, treatment_arm, fitted_beta)

#
#### Power as a function of different simulated effect sizes ####
remove_non_numeric <- function(x) {
      str_extract(x, "[0-9.]+") %>%
            as.numeric()
}

filenames <- 
      data.frame(filename = list.files(path = "simulations", 
                                       pattern = "perm.results")) %>% 
      filter(!grepl("strat", filename)) %>% 
      mutate(runname = gsub("perm.results", "", filename),
             runname = gsub(".csv", "", runname),
             setup = gsub("eff[0-9]_", "", runname)) %>% 
      separate(runname, 
               sep = "_", 
               into = c("samplesize", "cpgnumber", "effect", "nsim", "nperm", "resdist"),
               remove = F) %>% 
      mutate_at(vars(-filename, -runname, -setup), ~remove_non_numeric(.)) %>% 
      filter(samplesize != 100)


# Proportion of simulations that yielded a significant result:
filenames$g1_prop <- NA
filenames$g1_meanaatehat <- NA
filenames$g2_prop <- NA
filenames$g2_meanaatehat <- NA
filenames$g3_prop <- NA
filenames$g3_meanaatehat <- NA

# Reading data and calculating proportions:
p_cutoff <- 0.05
for(i in 1:nrow(filenames)){
      di <- read.csv(paste0("simulations/", filenames$filename[i]))
      props.i <- apply(di[,c(1,3,5)], 2, function(x) sum(x < p_cutoff) / length(x))
      filenames$g1_prop[i] <- props.i[1]
      filenames$g2_prop[i] <- props.i[2]
      filenames$g3_prop[i] <- props.i[3]
      
      filenames$g1_meanaatehat[i] <- mean(di$g1_AATEhat)
      filenames$g2_meanaatehat[i] <- mean(di$g2_AATEhat)
      filenames$g3_meanaatehat[i] <- mean(di$g3_AATEhat)
}

# Plotting: 
plotdata.1 <- filenames %>% 
      filter(samplesize == 0 & cpgnumber == 0) %>% 
      select(effect, g1_prop, g2_prop, g3_prop, resdist) %>% 
      pivot_longer(-c(effect, resdist)) %>% 
      mutate(`Residual distribution` = factor(resdist, 
                                              levels = c(0,1), 
                                              labels = c("Normal", "Empirical")),
             `g function` = gsub("_prop", "", name),
             effect.num = 1-exp((effect-1)/9))

# 80% power to detect: (using g1 and empirical residuals)
plotdata.1 %>% filter(effect %in% c(6,7) & resdist == 1 & `g function` == "g1")
# At effect -0.743, the power is 0.769, and at effect -0.948 the power is 0.92. 
x1 = -0.743
y1 = 0.769
x2 = -0.948
y2 = 0.92

y = 0.8
x = ((y - y1)*(x2-x1)/(y2-y1)) + x1

pdf("figures/Simulation_results_1.pdf", height = 4, width = 6)
plotdata.1 %>% ggplot(aes(x = effect.num, 
                          y = value, 
                          group = interaction(`g function`, `Residual distribution`), 
                          color = `g function`, 
                          linetype = `Residual distribution`)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      ylab("Proportion of rejected null hypotheses") +
      xlab("Simulated treatment effect") + 
      geom_segment(aes(x = -1.40, y = 0.8, xend = x, yend = 0.8), color = "gray", linetype = "dashed") +
      geom_segment(aes(x = x, y = -1, xend = x, yend = 0.8), color = "gray", linetype = "dashed") +
      coord_cartesian(xlim = c(-1.2, 0),
                      ylim = c(0,1))
dev.off()

# What would the effect size x = -0.7850861 look like in reality?

possible_effects <- c(0, -0.7850861)
small_predictions <- matrix(nrow = nrow(small), ncol = length(possible_effects))
for(i in 1:length(possible_effects)){
      beta_i <- c(beta[1], 0, beta[2:4])
      beta_i[2] <- possible_effects[i]
      pred_i <- small.X %*% beta_i
      small_predictions[, i] <- pred_i
}
colnames(small_predictions) <- paste("effect", 0:1)

pdf("figures/simulated_treatment_effects_at_80_power.pdf", width = 6, height = 3)
small_predictions %>% data.frame() %>% 
      bind_cols(small) %>% 
      pivot_longer(cols = starts_with("effect.")) %>% 
      filter(treatment_arm == "vitaminC" | name == "pred_0") %>%
      mutate(fitted_beta = inverse_logit(value),
             effect = if_else(name == "effect.0", "Effect 0", "Effect x")) %>% 
      ggplot(aes(x = visit, y = fitted_beta, group = effect, color = TET2_mut, linetype = effect)) +
      geom_line() +
      ylim(c(0,1)) +
      facet_wrap(~ TET2_mut)
dev.off()

# Numbers: 
small_predictions %>% data.frame() %>% 
      bind_cols(small) %>% 
      pivot_longer(cols = starts_with("effect.")) %>% 
      filter(treatment_arm == "vitaminC" | name == "pred_0") %>%
      mutate(fitted_beta = inverse_logit(value),
             effect = if_else(name == "effect.0", "Effect 0", "Effect x")) %>% 
      filter(visit == 5)
# WT: 0.443 - 0.267 = 0.176
# mut: 0.787 - 0.628 = 0.159

plotdata.1 %>% ggplot(aes(x = effect, 
                          y = value, 
                          group = interaction(`g function`, `residual distribution`), 
                          color = `g function`, 
                          linetype = `residual distribution`)) +
      geom_line() +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      facet_wrap(~name)

# Comparing different CpG numbers, g1 only:
filenames %>% 
      filter(cpgnumber != 0 & resdist == 1 & nsim == 1000) %>% 
      distinct(setup)
plotdata.2 <- filenames %>% 
      filter(cpgnumber != 0 & resdist == 1 & nsim == 1000) %>% 
      select(effect, g1_prop, cpgnumber) %>% 
      pivot_longer(-c(effect, cpgnumber)) %>% 
      mutate(cpgnumber = factor(cpgnumber),
             effect.num = 1-exp((effect-1)/9))

pdf("figures/simulations_CpGnumbers.pdf", width = 6, height = 5)
plotdata.2 %>% 
      mutate(`Number of CpG sites` = cpgnumber) %>% 
      ggplot(aes(x = effect.num, 
                 y = value, 
                 group = `Number of CpG sites`, 
                 color = `Number of CpG sites`)) +
      geom_line() +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      ylab("Proportion of rejected null hypotheses") +
      xlab("Simulated treatment effect")
dev.off()

# Comparing old vs new definition of CpG-specific effect
filenames %>% 
      filter((cpgnumber == 0 | cpgnumber == 720) & resdist == 1 & nsim == 1000) %>% 
      distinct(setup)
plotdata.3 <- filenames %>% 
      filter((cpgnumber == 0 | cpgnumber == 720) & resdist == 1 & nsim == 1000) %>% 
      select(effect, g1_prop, g2_prop, g3_prop, cpgnumber) %>% 
      pivot_longer(-c(effect, cpgnumber)) %>% 
      mutate(cpgnumber = factor(cpgnumber))

plotdata.3 %>% ggplot(aes(x = effect, 
                          y = value, 
                          group = interaction(name, cpgnumber), 
                          color = cpgnumber, 
                          linetype = name)) +
      geom_line() +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      facet_wrap(~ name)

#### Computation times ####

cpt.filenames <- 
      data.frame(filename = list.files(path = "simulations", 
                                       pattern = "perm.comptime")) %>% 
      mutate(runname = gsub("perm.comptime_", "", filename),
             runname = gsub(".csv", "", runname),
             setup = gsub("eff[0-9]_", "", runname)) %>% 
      separate(runname, 
               sep = "_", 
               into = c("samplesize", "cpgnumber", "effect", "nsim", "nperm", "resdist"),
               remove = F) %>% 
      mutate_at(vars(-filename, -runname, -setup), ~remove_non_numeric(.)) %>% 
      filter(samplesize != 100)

cpt.filenames$comptime <- NA

# Reading data and calculating proportions:
for(i in 1:nrow(cpt.filenames)){
      di <- read.csv(paste0("simulations/", cpt.filenames$filename[i]))
      cpt.filenames$comptime[i] <- di$x[1]
}

# First experiments:
cpt.filenames %>% 
      filter(samplesize == 0 & cpgnumber == 0) %>% 
      summarise(meantime = mean(comptime))
# 3545 minutes = 59 hours

# CpG number experiments:
cpt.filenames %>% 
      filter(cpgnumber != 0 & resdist == 1 & nsim == 1000) %>% 
      group_by(cpgnumber) %>% 
      summarise(meantime = mean(comptime)) %>% 
      mutate(meantime_hours = meantime / 60)


#### Modeling CpG effect ####

cpg_numbers <- c(10, 50, 100, 300, 720, 1500, 3000, 5000, 10000)

mnc <- vroom::vroom("mnc_results_commonprobes.csv")
gran <- vroom::vroom("granulocyte_results_commonprobes.csv")

real.estimates <- {}
pos.estimates <- {}
for(i in 1:length(cpg_numbers)){
      K = cpg_numbers[i]
      topKmnc <- mnc %>% filter(estimate > 0) %>% 
            slice_tail(n = K) %>% 
            pull(probe)
      topKgran <- gran %>% 
            filter(probe %in% topKmnc) %>% 
            mutate(estimate0 = if_else(estimate < 0, 0, estimate)) 
      topKgran.real <- topKgran %>% pull(estimate)
      topKgran.pos <- topKgran %>% pull(estimate0)
      topKgran.real.padded <- c(topKgran.real, rep(NA, 10000 - K))
      topKgran.pos.padded <- c(topKgran.pos, rep(NA, 10000 - K))
      real.estimates <- bind_cols(real.estimates, topKgran.real.padded)
      pos.estimates <- bind_cols(scaled.estimates, topKgran.pos.padded)
      names(real.estimates)[i] <- paste0("CpG.", K)
      names(pos.estimates)[i] <- paste0("CpG.", K)
}

pdf("figures/CpG_gran_effects.pdf", height = 4, width = 6)
real.estimates %>% 
      select(c(1, 2,4,5,6,9)) %>% 
      pivot_longer(everything()) %>% 
      mutate(`Number of CpG sites` = as.numeric(gsub("CpG.", "", name)) %>% factor()) %>% 
      ggplot(aes(x = value, group = `Number of CpG sites`, color = `Number of CpG sites`)) +
      geom_density() +
      xlab("TET2 mutation effect in granulocyte dataset")
dev.off()

pos.estimates %>% 
      pivot_longer(everything()) %>% 
      mutate(n_cpgs = as.numeric(gsub("CpG.", "", name)) %>% factor()) %>% 
      ggplot(aes(x = value, group = n_cpgs, color = n_cpgs)) +
      geom_density()

# Selecting positively associated probes with lowest P values in MNC dataset
topmnc <- vroom::vroom("mnc_results_commonprobes.csv") %>% 
      arrange(desc(pval)) %>% 
      filter(estimate > 0) %>% 
      slice_tail(n = n_cpg) %>% 
      pull(probe)
# Finding these probes' coefficients from EWAS with TET2 mutation status in the granulocyte dataset
# These coefficients are then scaled with a factor 1/1.82, since 1.79 is the mean of the granulocyte coefficients
# for the 500 top probes from the mnc dataset.
# Estimate < 0 are set to 0
site_effect <- vroom::vroom("granulocyte_results_commonprobes.csv") %>% 
      filter(probe %in% topmnc) %>% 
      transmute(coef_scaled = (estimate > 0) * estimate / 1.82,
                probeID = factor(1:n_cpg))

# Altering the treatment_constrained variable s.t. the design matrix will instead contain the scaled coefficients as treatment effects
dsim <- dsim %>% 
      left_join(site_effect) %>% 
      mutate(treatment_constrained = if_else(treatment_constrained == "baseline", 0, coef_scaled))


#### Stratified models for simulations ####

# Models: mutated subset
dlong.mut <- dlong %>% filter(TET2_mut == "mut")

mnull4.mut <- lmer(Mval ~ visit + (1|patient_visit) + (1|probeID), data = dlong.mut)
mtreat4.mut <- lmer(Mval ~ treatment_constrained + visit + (1|patient_visit) + (1|probeID), data = dlong.mut)

summary(mnull4.mut)
summary(mtreat4.mut)

dlong.mut$mnull4.mut.resid <- summary(mnull4.mut)$residuals
dlong.mut$mnull4.mut.fitted <- fitted(mnull4.mut)

resid.plot <- dlong.mut %>% 
      ggplot(aes(x = mnull4.mut.fitted, y = mnull4.mut.resid)) +
      geom_point(alpha = 0.1) +
      geom_smooth() +
      xlab("Fitted values") +
      ylab("Residuals")

qq <- dlong.mut %>% ggplot(aes(sample = mnull4.mut.resid)) +
      stat_qq() + 
      stat_qq_line() +
      xlab("Theoretical quantiles") +
      ylab("Sample quantiles")

pdf("figures/Diagnostics_mut.stratified.pdf", width = 6, height = 3)
ggarrange(resid.plot, qq, 
          nrow = 1, ncol = 2, labels = c("A", "B"))
dev.off()

# Models: WT subset
dlong.wt <- dlong %>% filter(TET2_mut == "WT")
mnull4.wt <- lmer(Mval ~ visit + (1|patient_visit) + (1|probeID), data = dlong.wt)
mtreat4.wt <- lmer(Mval ~ treatment_constrained + visit + (1|patient_visit) + (1|probeID), data = dlong.wt)

summary(mnull4.wt)
summary(mtreat4.wt)

dlong.wt$mnull4.wt.resid <- summary(mnull4.wt)$residuals
dlong.wt$mnull4.wt.fitted <- fitted(mnull4.wt)

resid.plot <- dlong.wt %>% 
      ggplot(aes(x = mnull4.wt.fitted, y = mnull4.wt.resid)) +
      geom_point(alpha = 0.1) +
      geom_smooth() +
      xlab("Fitted values") +
      ylab("Residuals")

qq <- dlong.wt %>% ggplot(aes(sample = mnull4.resid)) +
      stat_qq() + 
      stat_qq_line() +
      xlab("Theoretical quantiles") +
      ylab("Sample quantiles")

pdf("figures/Diagnostics_WT.stratified.pdf", width = 6, height = 3)
ggarrange(resid.plot, qq, 
          nrow = 1, ncol = 2, labels = c("A", "B"))
dev.off()

# Making the same modelplot based on simulations from each of the two models:

# modelobject = mnull1
dsim.mut <- dlong.mut %>% 
      mutate(Mval = simulate(mnull4.mut)$sim_1,
             beta = inverse_logit(Mval))
dsim.wt <- dlong.wt %>% 
      mutate(Mval = simulate(mnull4.wt)$sim_1,
             beta = inverse_logit(Mval))

dsim <- bind_rows(dsim.mut, dsim.wt)

beta_hist <- dsim %>%
      filter(visit == 1) %>%
      ggplot(aes(x = beta, fill = TET2_mut)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 100) +
      facet_wrap(~ TET2_mut) +
      theme(legend.position = "none")

# Panel B: Four representative CpG sites
foursites <- c("cg12275410", "cg01498832", "cg13747899", "cg00782708")
foursites_hist <- dsim %>%
      filter(probeID %in% foursites) %>%
      ggplot(aes(x = beta, fill = TET2_mut, group = TET2_mut)) +
      geom_histogram(position = "dodge", bins = 20) +
      facet_wrap(~ probeID) +
      theme(legend.position = "none")

# Panel C: Distributions of pairwise Pearson correlation coefficients between baseline beta values in TET2 mut and WT subsets
beta_df.base.mut <- dsim %>%
      filter(visit == 1 & TET2_mut == "mut") %>%
      select(arrayID, probeID, beta) %>%
      pivot_wider(values_from = beta, names_from = probeID) %>%
      column_to_rownames("arrayID")

beta_df.base.wt <- dsim %>%
      filter(visit == 1 & TET2_mut == "WT") %>%
      select(arrayID, probeID, beta) %>%
      pivot_wider(values_from = beta, names_from = probeID) %>%
      column_to_rownames("arrayID")

cor.base.mut <- cor(beta_df.base.mut)
cor.base.mut[!upper.tri(cor.base.mut)] <- NA
cor.base.mut.df <- data.frame(cor = as.vector(cor.base.mut),
                              TET2_mut = "mut") %>%
      filter(!is.na(cor))

cor.base.wt <- cor(beta_df.base.wt)
cor.base.wt[!upper.tri(cor.base.wt)] <- NA
cor.base.wt.df <- data.frame(cor = as.vector(cor.base.wt),
                             TET2_mut = "WT") %>%
      filter(!is.na(cor))

cor_hist <- bind_rows(cor.base.mut.df, cor.base.wt.df) %>%
      mutate(TET2_mut = factor(TET2_mut, levels = c("WT", "mut"))) %>%
      ggplot(aes(x = cor, fill = TET2_mut)) +
      geom_histogram(bins = 100, aes(y = after_stat(density))) +
      facet_wrap( ~ TET2_mut) +
      theme(legend.position = "none") +
      xlim(c(-1,1))

# Panel D: Distributions of pairwise Pearson correlation coefficients between changes in beta values from baseline to end of study
#  in TET2 mut and WT subsets
beta_df.end.mut <- dsim %>%
      filter(visit == 5 & TET2_mut == "mut") %>%
      select(arrayID, probeID, beta) %>%
      pivot_wider(values_from = beta, names_from = probeID) %>%
      column_to_rownames("arrayID")

beta_df.end.wt <- dsim %>%
      filter(visit == 5 & TET2_mut == "WT") %>%
      select(arrayID, probeID, beta) %>%
      pivot_wider(values_from = beta, names_from = probeID) %>%
      column_to_rownames("arrayID")

z_df.mut <- beta_df.end.mut - beta_df.base.mut
z_df.wt <- beta_df.end.wt - beta_df.base.wt

cor.z.mut <- cor(z_df.mut)
cor.z.mut[!upper.tri(cor.z.mut)] <- NA
cor.z.mut.df <- data.frame(cor = as.vector(cor.z.mut),
                           TET2_mut = "mut") %>%
      filter(!is.na(cor))

cor.z.wt <- cor(z_df.wt)
cor.z.wt[!upper.tri(cor.z.wt)] <- NA
cor.z.wt.df <- data.frame(cor = as.vector(cor.z.wt),
                          TET2_mut = "WT") %>%
      filter(!is.na(cor))

cor.z_hist <- bind_rows(cor.z.mut.df, cor.z.wt.df) %>%
      mutate(TET2_mut = factor(TET2_mut, levels = c("WT", "mut"))) %>%
      ggplot(aes(x = cor, fill = TET2_mut)) +
      geom_histogram(bins = 100, aes(y = after_stat(density))) +
      facet_wrap( ~ TET2_mut) +
      theme(legend.position = "none") +
      xlim(c(-1,1))
t.beta_df.base.mut <- t(beta_df.base.mut)
t.beta_df.base.wt <- t(beta_df.base.wt)
cor.base.mut <- cor(t.beta_df.base.mut)
cor.base.mut[!upper.tri(cor.base.mut)] <- NA
cor.base.mut.df <- data.frame(cor = as.vector(cor.base.mut),
                              TET2_mut = "mut") %>%
      filter(!is.na(cor))

cor.base.wt <- cor(t.beta_df.base.wt)
cor.base.wt[!upper.tri(cor.base.wt)] <- NA
cor.base.wt.df <- data.frame(cor = as.vector(cor.base.wt),
                             TET2_mut = "WT") %>%
      filter(!is.na(cor))

cor_hist.t <- bind_rows(cor.base.mut.df, cor.base.wt.df) %>%
      mutate(TET2_mut = factor(TET2_mut, levels = c("WT", "mut"))) %>%
      ggplot(aes(x = cor, fill = TET2_mut)) +
      geom_histogram(bins = 100, aes(y = after_stat(density))) +
      facet_wrap( ~ TET2_mut) +
      theme(legend.position = "none") +
      xlim(c(-1,1))

plot <- list(beta_hist, foursites_hist, cor_hist, cor.z_hist, cor_hist.t)

pdf(paste0("figures/", "stratified_simulationplot_mnull4.strat", ".pdf"))
ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], 
          nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
dev.off()





#### Stratified simulations ####
filenames.strat <- 
      data.frame(filename = list.files(path = "simulations", 
                                       pattern = "perm.results")) %>% 
      filter(grepl("strat", filename)) %>% 
      mutate(runname = gsub("perm.results", "", filename),
             runname = gsub(".csv", "", runname),
             setup = gsub("eff[0-9]_", "", runname)) %>% 
      separate(runname, 
               sep = "_", 
               into = c("stratum", "cpgnumber", "effect", "nsim", "nperm"),
               remove = F) %>% 
      mutate_at(vars(-filename, -runname, -setup, -stratum), ~remove_non_numeric(.)) %>% 
      mutate(stratum = gsub("strat", "", stratum))


# Proportion of simulations that yielded a significant result:
filenames.strat$g1_prop <- NA
filenames.strat$g1_meanaatehat <- NA
filenames.strat$g2_prop <- NA
filenames.strat$g2_meanaatehat <- NA
filenames.strat$g3_prop <- NA
filenames.strat$g3_meanaatehat <- NA

# Reading data and calculating proportions:
p_cutoff <- 0.05
for(i in 1:nrow(filenames.strat)){
      di <- read.csv(paste0("simulations/", filenames.strat$filename[i]))
      props.i <- apply(di[,c(1,3,5)], 2, function(x) sum(x < p_cutoff) / length(x))
      filenames.strat$g1_prop[i] <- props.i[1]
      filenames.strat$g2_prop[i] <- props.i[2]
      filenames.strat$g3_prop[i] <- props.i[3]
      
      filenames.strat$g1_meanaatehat[i] <- mean(di$g1_AATEhat)
      filenames.strat$g2_meanaatehat[i] <- mean(di$g2_AATEhat)
      filenames.strat$g3_meanaatehat[i] <- mean(di$g3_AATEhat)
}

# Plotting: 
plotdata.1 <- filenames.strat %>% 
      filter(cpgnumber == 50 & nsim == 1000) %>% 
      select(effect, g1_prop, stratum) %>% 
      pivot_longer(-c(effect, stratum)) %>% 
      mutate(`g function` = gsub("_prop", "", name),
             effect.num = 1-exp((effect-1)/9))

# 80% power to detect:
plotdata.1 %>% filter(effect %in% c(7,8) & `g function` == "g1")
# TET2 mut: 
x1 = plotdata.1 %>% filter(effect == 7 & `g function` == "g1" & stratum == "mut") %>% pull(effect.num)
y1 = plotdata.1 %>% filter(effect == 7 & `g function` == "g1" & stratum == "mut") %>% pull(value)
x2 = plotdata.1 %>% filter(effect == 8 & `g function` == "g1" & stratum == "mut") %>% pull(effect.num)
y2 = plotdata.1 %>% filter(effect == 8 & `g function` == "g1" & stratum == "mut") %>% pull(value)

y.mut = 0.8
x.mut = ((y.mut - y1)*(x2-x1)/(y2-y1)) + x1
# 80% power at x = -1.14

powerplot.mut <- plotdata.1 %>% 
      filter(stratum == "mut") %>% 
      ggplot(aes(x = effect.num, 
                 y = value )) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      ylab("Proportion of rejected null hypotheses") +
      xlab("Simulated treatment effect") + 
      geom_segment(aes(x = -1.40, y = 0.8, xend = x.mut, yend = 0.8), color = "gray", linetype = "dashed") +
      geom_segment(aes(x = x.mut, y = -1, xend = x.mut, yend = 0.8), color = "gray", linetype = "dashed") +
      coord_cartesian(xlim = c(-1.2, 0),
                      ylim = c(0,1))

# TET2 WT: 
# For now, the power is simply 0.8 directly at effect = -0.948
x1 = plotdata.1 %>% filter(effect == 7 & `g function` == "g1" & stratum == "WT") %>% pull(effect.num)
y1 = plotdata.1 %>% filter(effect == 7 & `g function` == "g1" & stratum == "WT") %>% pull(value)
x2 = plotdata.1 %>% filter(effect == 8 & `g function` == "g1" & stratum == "WT") %>% pull(effect.num)
y2 = plotdata.1 %>% filter(effect == 8 & `g function` == "g1" & stratum == "WT") %>% pull(value)
y.wt = 0.8
x.wt = ((y.wt - y1)*(x2-x1)/(y2-y1)) + x1
# 80% power at -1.05

powerplot.wt <- plotdata.1 %>% 
      filter(stratum == "WT") %>% 
      ggplot(aes(x = effect.num, 
                 y = value)) +
      geom_line() +
      geom_point() +
      geom_hline(yintercept = 0.05, linetype = "dashed") +
      ylab("Proportion of rejected null hypotheses") +
      xlab("Simulated treatment effect") + 
      geom_segment(aes(x = -1.40, y = 0.8, xend = x.wt, yend = 0.8), color = "gray", linetype = "dashed") +
      geom_segment(aes(x = x.wt, y = -1, xend = x.wt, yend = 0.8), color = "gray", linetype = "dashed") +
      coord_cartesian(xlim = c(-1.2, 0),
                      ylim = c(0,1))

pdf("figures/Simulation_results_stratified.pdf", height = 4, width = 6)
ggarrange(powerplot.wt, powerplot.mut, 
          nrow = 1, ncol = 2, labels = c("A", "B"))
dev.off()

#### Effects at beta scale at 80% power in stratified analyses ####

# WT patients:
possible_effects <- c(0, x.wt, -0.7850861)
small.wt <- dlong %>% 
      filter(TET2_mut == "WT") %>% 
      select(visit, treatment_arm, treatment_constrained) %>% 
      data.frame() %>% 
      distinct()
small.wt.X <- model.matrix( ~ treatment_constrained + visit, data = small.wt)
beta.wt <- summary(mnull4.wt)$coefficients[,1]

small_predictions.wt <- matrix(nrow = nrow(small.wt), ncol = length(possible_effects))
for(i in 1:length(possible_effects)){
      beta_i <- c(beta.wt[1], 0, beta.wt[2])
      beta_i[2] <- possible_effects[i]
      pred_i <- small.wt.X %*% beta_i
      small_predictions.wt[, i] <- pred_i
}
colnames(small_predictions.wt) <- paste("effect", 0:2)

effectplot.wt <- small_predictions.wt %>% 
      data.frame() %>% 
      bind_cols(small.wt) %>% 
      pivot_longer(cols = starts_with("effect.")) %>% 
      filter(treatment_arm == "vitaminC") %>%
      mutate(fitted_beta = inverse_logit(value),
             effect = case_when(name == "effect.0" ~ "Placebo",
                                name == "effect.1" ~ "-1.05",
                                name == "effect.2" ~ "-0.79")) %>% 
      ggplot(aes(x = visit, y = fitted_beta, group = effect, color = effect)) +
      geom_line() +
      ylim(c(0,1)) +
      ylab(NULL)

# mutated patients:
# WT patients:
possible_effects <- c(0, x.mut, -0.7850861)
small.mut <- dlong %>% 
      filter(TET2_mut == "mut") %>% 
      select(visit, treatment_arm, treatment_constrained) %>% 
      data.frame() %>% 
      distinct()
small.mut.X <- model.matrix( ~ treatment_constrained + visit, data = small.mut)
beta.mut <- summary(mnull4.mut)$coefficients[,1]

small_predictions.mut <- matrix(nrow = nrow(small.mut), ncol = length(possible_effects))
for(i in 1:length(possible_effects)){
      beta_i <- c(beta.mut[1], 0, beta.mut[2])
      beta_i[2] <- possible_effects[i]
      pred_i <- small.mut.X %*% beta_i
      small_predictions.mut[, i] <- pred_i
}
colnames(small_predictions.mut) <- paste("effect", 0:2)

effectplot.mut <- small_predictions.mut %>% 
      data.frame() %>% 
      bind_cols(small.mut) %>% 
      pivot_longer(cols = starts_with("effect.")) %>% 
      filter(treatment_arm == "vitaminC") %>%
      mutate(fitted_beta = inverse_logit(value),
             effect = case_when(name == "effect.0" ~ "Placebo",
                                name == "effect.1" ~ "-1.14",
                                name == "effect.2" ~ "-0.79")) %>% 
      ggplot(aes(x = visit, y = fitted_beta, group = effect, color = effect)) +
      geom_line() +
      ylim(c(0,1))  +
      ylab(NULL)


pdf("figures/simulated_treatment_effects_at_80_power_stratified.pdf", width = 6, height = 4)
ggarrange(effectplot.wt, effectplot.mut, 
          nrow = 1, ncol = 2, labels = c("A", "B"))
dev.off()

# Numbers: 
small_predictions %>% data.frame() %>% 
      bind_cols(small) %>% 
      pivot_longer(cols = starts_with("effect.")) %>% 
      filter(treatment_arm == "vitaminC" | name == "pred_0") %>%
      mutate(fitted_beta = inverse_logit(value),
             effect = if_else(name == "effect.0", "Effect 0", "Effect x")) %>% 
      filter(visit == 5)
# WT: 0.443 - 0.267 = 0.176
# mut: 0.787 - 0.628 = 0.159


###### Conclusions: AATE and permutation sampling ####
# Now restricting to top 50 CpG sites:
load(file = "trial_data_top50mncsites.Rdata")

#### Main analyses (entire dataset, top50 mnc probes) ####
g <- function(x) x
calculate_AATE <- function(data, varname, g_function = g){
      visit1_vitamin_aids <- data %>% filter(visit == 1 & treatment_arm == "vitaminC") %>% select(arrayID) %>% distinct() %>% pull(arrayID)
      visit5_vitamin_aids <- data %>% filter(visit == 5 & treatment_arm == "vitaminC") %>% select(arrayID) %>% distinct() %>% pull(arrayID)
      visit1_placebo_aids <- data %>% filter(visit == 1 & treatment_arm == "placebo") %>% select(arrayID) %>% distinct() %>% pull(arrayID)
      visit5_placebo_aids <- data %>% filter(visit == 5 & treatment_arm == "placebo") %>% select(arrayID) %>% distinct() %>% pull(arrayID)
      
      dwide <- data %>% select(all_of(c("probeID", "arrayID", varname))) %>% 
            pivot_wider(values_from = all_of(varname), names_from = probeID) %>% 
            column_to_rownames(var = "arrayID") %>% 
            as.matrix()
      
      Y_1_vitamin <- dwide[visit1_vitamin_aids, ]
      Y_1_placebo <- dwide[visit1_placebo_aids, ]
      Y_2_vitamin <- dwide[visit5_vitamin_aids, ]
      Y_2_placebo <- dwide[visit5_placebo_aids, ]
      
      Z_vitamin <- Y_2_vitamin - Y_1_vitamin
      Z_placebo <- Y_2_placebo - Y_1_placebo
      
      Z_vitamin_means <- apply(Z_vitamin, 2, mean)
      Z_placebo_means <- apply(Z_placebo, 2, mean)
      
      ATE_k <- Z_vitamin_means - Z_placebo_means
      ATE_k <- g_function(ATE_k)
      return(mean(ATE_k))
}
AATE_hat <- calculate_AATE(data = dlong.top50, "beta", g_function = g)

nperm <- 10000
df.AATE_perm <- data.frame(AATE_perm = rep(NA, nperm))

dlong.top50.perm <- dlong.top50
N <- length(unique(dlong.top50$patient_id))
N_vitc <- length(unique(dlong.top50$patient_id[dlong.top50$treatment_arm == "vitaminC"]))
pids <- unique(dlong.top50$patient_id)
for(b in 1:nperm){
      sample_b <- sample(x = 1:N, size = N_vitc, replace = F)
      dlong.top50.perm$treatment_arm <- if_else(dlong.top50.perm$patient_id %in% pids[sample_b], "vitaminC", "placebo")
      df.AATE_perm$AATE_perm[b] <- calculate_AATE(data = dlong.top50.perm, "beta", g_function = g)
}
save(AATE_hat, df.AATE_perm, file = "main_AATE_permutation.Rdata")
load(file = "main_AATE_permutation.Rdata")

pdf("figures/AATE_permutation.pdf", height = 4, width = 6)
df.AATE_perm %>% ggplot(aes(x = AATE_perm)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 30) +
      geom_vline(xintercept = AATE_hat, color = "red") +
      xlab("AATE")
dev.off()
# One-sided P: 0.58
mean(df.AATE_perm$AATE_perm < AATE_hat)
# Two-sided P: 0.82
mean(abs(df.AATE_perm$AATE_perm) > AATE_hat)


#### Stratified: TET2 mut ####
dlong.top50.mut <- dlong.top50 %>% filter(TET2_mut == "mut")
dlong.top50.mut.perm <- dlong.top50.mut
N.mut <- length(unique(dlong.top50.mut$patient_id))
N.mut_vitc <- length(unique(dlong.top50.mut$patient_id[dlong.top50.mut$treatment_arm == "vitaminC"]))
pids.mut <- unique(dlong.top50.mut$patient_id)

AATE_hat.mut <- calculate_AATE(data = dlong.top50.mut, "beta", g_function = g)
df.AATE_perm.mut <- data.frame(AATE_perm = rep(NA, nperm))
for(b in 1:nperm){
      sample_b <- sample(x = 1:N.mut, size = N.mut_vitc, replace = F)
      dlong.top50.mut.perm$treatment_arm <- if_else(dlong.top50.mut.perm$patient_id %in% pids[sample_b], "vitaminC", "placebo")
      df.AATE_perm.mut$AATE_perm[b] <- calculate_AATE(data = dlong.top50.mut.perm, "beta", g_function = g)
}
save(AATE_hat.mut, df.AATE_perm.mut, file = "mut_AATE_permutation.Rdata")

load(file = "mut_AATE_permutation.Rdata")
permplot.mut <- df.AATE_perm.mut %>% ggplot(aes(x = AATE_perm)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 30) +
      geom_vline(xintercept = AATE_hat.mut, color = "red") +
      xlab("AATE")
# One-sided P: 0.47
mean(df.AATE_perm.mut$AATE_perm < AATE_hat.mut)
# Two-sided P: 0.74
mean(abs(df.AATE_perm.mut$AATE_perm) > AATE_hat.mut)


#### Stratified: TET2 WT ####
dlong.top50.wt <- dlong.top50 %>% filter(TET2_mut == "WT")
dlong.top50.wt.perm <- dlong.top50.wt
N.wt <- length(unique(dlong.top50.wt$patient_id))
N.wt_vitc <- length(unique(dlong.top50.wt$patient_id[dlong.top50.wt$treatment_arm == "vitaminC"]))
pids.wt <- unique(dlong.top50.wt$patient_id)

AATE_hat.wt <- calculate_AATE(data = dlong.top50.wt, "beta", g_function = g)
df.AATE_perm.wt <- data.frame(AATE_perm = rep(NA, nperm))
for(b in 1:nperm){
      sample_b <- sample(x = 1:N.wt, size = N.wt_vitc, replace = F)
      dlong.top50.wt.perm$treatment_arm <- if_else(dlong.top50.wt.perm$patient_id %in% pids[sample_b], "vitaminC", "placebo")
      df.AATE_perm.wt$AATE_perm[b] <- calculate_AATE(data = dlong.top50.wt.perm, "beta", g_function = g)
}
save(AATE_hat.wt, df.AATE_perm.wt, file = "WT_AATE_permutation.Rdata")
load(file = "WT_AATE_permutation.Rdata")
permplot.WT <- df.AATE_perm.wt %>% ggplot(aes(x = AATE_perm)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 30) +
      geom_vline(xintercept = AATE_hat.wt, color = "red") +
      xlab("AATE")
# One-sided P: 0.13
mean(df.AATE_perm.wt$AATE_perm < AATE_hat.wt)
# Two-sided P: 0.97
mean(abs(df.AATE_perm.wt$AATE_perm) > AATE_hat.wt)

pdf("figures/AATE_permutation_stratified.pdf", height = 4, width = 6)
ggarrange(permplot.WT, permplot.mut, 
          nrow = 1, ncol = 2, labels = c("A", "B"))

dev.off()

# 
#### For discussion: Correlations between CpG sites in different models ####

summary(mnull1)

summary(mnull4)
summary(mnull5)

# Random effects:
# Groups     Name        Variance Std.Dev.
# probe_mut  (Intercept) 0.15887  0.3986  
# probeID    (Intercept) 1.09427  1.0461  
# arrayID    (Intercept) 0.09337  0.3056  
# patient_id (Intercept) 0.85649  0.9255  
# Residual               0.38490  0.6204

# Variances:
v.probe_mut  <- 0.3986
v.probeID    <- 1.0461
v.arrayID    <- 0.3056
v.patient_id <- 0.9255
v.residual   <- 0.6204

# Correlation between two M values at the same CpG site in different times in the same patient:
(v.probe_mut + v.probeID + v.patient_id) / (v.probe_mut + v.probeID + v.arrayID + v.patient_id + v.residual)
# = 0.72

# Correlation between two M values at DIFFERENT CpG sites in the same time in the same patient:
(v.patient_id + v.arrayID) / (v.probe_mut + v.probeID + v.arrayID + v.patient_id + v.residual)
# = 0.37

# Correlation between two M values at DIFFERENT CpG sites in different times in the same patient:
(v.patient_id) / (v.probe_mut + v.probeID + v.arrayID + v.patient_id + v.residual)
# = 0.28

# Correlation between two M values at the same site in different patients with different mutation status
(v.probeID) / (v.probe_mut + v.probeID + v.arrayID + v.patient_id + v.residual)
# = 0.32

# Correlation between two M values at the same site in different patients with same mutation status
(v.probeID + v.probe_mut) / (v.probe_mut + v.probeID + v.arrayID + v.patient_id + v.residual)
# = 0.44

#### ####





































