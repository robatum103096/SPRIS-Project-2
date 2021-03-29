rm(list=ls())

library(tidyverse)
library(dplyr)
library(nlme)
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(qwraps2)
dpi <- 800

secdata <- read.csv("secdata.csv")

carryover <- c()
for (i in 1:nrow(secdata)) {
  if (secdata$sequence[i] == 1) {
    if (secdata$period[i] == "Period 1") carryover <- c(carryover, "None")
    if (secdata$period[i] == "Period 2") carryover <- c(carryover, "Pill A")
    if (secdata$period[i] == "Period 3") carryover <- c(carryover, "Gel B")
  }
  if (secdata$sequence[i] == 2) {
    if (secdata$period[i] == "Period 1") carryover <- c(carryover, "None")
    if (secdata$period[i] == "Period 2") carryover <- c(carryover, "Gel C")
    if (secdata$period[i] == "Period 3") carryover <- c(carryover, "Pill A")
  }
  if (secdata$sequence[i] == 3) {
    if (secdata$period[i] == "Period 1") carryover <- c(carryover, "None")
    if (secdata$period[i] == "Period 2") carryover <- c(carryover, "Gel B")
    if (secdata$period[i] == "Period 3") carryover <- c(carryover, "Gel C")
  }
  if (secdata$sequence[i] == 4) {
    if (secdata$period[i] == "Period 1") carryover <- c(carryover, "None")
    if (secdata$period[i] == "Period 2") carryover <- c(carryover, "Gel B")
    if (secdata$period[i] == "Period 3") carryover <- c(carryover, "Pill A")
  }
  if (secdata$sequence[i] == 5) {
    if (secdata$period[i] == "Period 1") carryover <- c(carryover, "None")
    if (secdata$period[i] == "Period 2") carryover <- c(carryover, "Pill A")
    if (secdata$period[i] == "Period 3") carryover <- c(carryover, "Gel C")
  }
  if (secdata$sequence[i] == 6) {
    if (secdata$period[i] == "Period 1") carryover <- c(carryover, "None")
    if (secdata$period[i] == "Period 2") carryover <- c(carryover, "Gel C")
    if (secdata$period[i] == "Period 3") carryover <- c(carryover, "Gel B")
  }
}
secdata$carryover <- factor(carryover, levels = c("None", "Pill A", "Gel B", "Gel C"))

secdata$AEsum <- as.numeric(secdata$AEsum > 0)
secdata$treatment <- factor(secdata$treatment, levels = c("Pill A", "Gel B", "Gel C"))
secdata$regimen <- ifelse(secdata$ treatment == "Gel B", "Three Times", "Once")
secdata$pill <- ifelse(secdata$treatment == "Pill A", "Yes", "No")

secdata$gender <- as.factor(secdata$gender)
levels(secdata$gender) <- c("Male", "Female")

df_mat <- model.matrix(~., data = secdata[,5:13])
df_cor <- cor(df_mat[,-1])
png("./corplot.png", width=8*dpi, height=6*dpi, res=dpi)
pheatmap(df_cor, display_numbers = T, angle_col = "45")
dev.off()

secdata_long <- rbind(data.frame(select(secdata, -viral_skin, -viral_blood),
                                 viral = secdata$viral_skin,
                                 position = rep("Skin", length(secdata$viral_skin))),
                      data.frame(select(secdata, -viral_skin, -viral_blood),
                                 viral = secdata$viral_blood,
                                 position = rep("Blood", length(secdata$viral_blood))))

# PK distributions by treatment and period
png("pk_period.png", width = 7*dpi, height = 4*dpi, res=dpi)
ggplot(data = secdata_long, aes(x = period, y = viral)) +
  geom_violin(aes(fill = treatment), position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(group=interaction(treatment,period)),
               width=0.1, fill="white", position = position_dodge(width = 0.8),
               outlier.shape=NA, outlier.color=NA) +
  facet_grid(~position) +
  ylab("Viral Load Change") +
  xlab("Period") +
  guides(fill=guide_legend(title="Treatment"))
dev.off()

# PK distributions by treatment, adherence
png("pk_adherence.png", width = 7*dpi, height = 4*dpi, res=dpi)
ggplot(data = secdata_long, aes(x = Adheresum, y = viral)) +
  geom_point(aes(color = treatment), alpha = 0.5) +
  geom_smooth(aes(color = treatment), method = lm, se = F) +
  facet_grid(~position) +
  ylab("Viral Load Change") +
  xlab("Adherence") +
  guides(color=guide_legend(title="Treatment"))
dev.off()

# PK distributions by treatment, ae
png("pk_ae.png", width = 7*dpi, height = 4*dpi, res=dpi)
ggplot(data = secdata_long) +
  geom_violin(aes(x = factor(AEsum), y = viral, fill = treatment), 
              position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = factor(AEsum), y = viral,
                   group=interaction(factor(AEsum),treatment)),
               width=0.1, fill="white", position = position_dodge(width = 0.8),
               outlier.shape=NA, outlier.color=NA) +
  geom_smooth(aes(x = AEsum+1, y = viral, color = treatment), method = lm, se = F) +
  facet_grid(~position) +
  ylab("Viral Load Change") +
  xlab("Adverse Event") +
  guides(color=guide_legend(title="Treatment"),
         fill=guide_legend(title="Treatment"))
dev.off()

# PK distributions by treatment, age
png("pk_age.png", width = 7*dpi, height = 4*dpi, res=dpi)
ggplot(data = secdata_long, aes(x = age, y = viral)) +
  geom_point(aes(color = treatment), alpha = 0.5) +
  geom_smooth(aes(color = treatment), method = lm, se = F) +
  facet_grid(~position) +
  ylab("Viral Load Change") +
  xlab("Age") +
  guides(color=guide_legend(title="Treatment"))
dev.off()

# adherence by treatment, age
png("adherence_age.png", width = 7*dpi, height = 4*dpi, res=dpi)
ggplot(data = secdata_long) +
  geom_point(aes(x = age, y = Adheresum, color = treatment), alpha = 0.5) +
  geom_smooth(aes(x = age, y = Adheresum, color = treatment), alpha = 0.5) +
  ylab("Adherence") +
  xlab("Age") +
  guides(color=guide_legend(title="Treatment"))
dev.off()

# adherence by treatment, gender
png("adherence_gender.png", width = 7*dpi, height = 4*dpi, res=dpi)
ggplot(data = secdata_long) +
  geom_violin(aes(x = gender, y = viral, fill = treatment), 
              position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = gender, y = viral,
                   group=interaction(gender,treatment)),
               width=0.1, fill="white", position = position_dodge(width = 0.8),
               outlier.shape=NA, outlier.color=NA) +
  facet_grid(~position) +
  ylab("Adherence") +
  xlab("Gender") +
  guides(color=guide_legend(title="Treatment"),
         fill=guide_legend(title="Treatment"))
dev.off()

# adherence by treatment, race
png("adherence_race.png", width = 7*dpi, height = 4*dpi, res=dpi)
ggplot(data = secdata_long) +
  geom_violin(aes(x = race, y = viral, fill = treatment), 
              position = position_dodge(width = 0.8)) +
  geom_boxplot(aes(x = race, y = viral,
                   group=interaction(race,treatment)),
               width=0.1, fill="white", position = position_dodge(width = 0.8),
               outlier.shape=NA, outlier.color=NA) +
  facet_grid(~position) +
  ylab("Adherence") +
  xlab("race") +
  guides(color=guide_legend(title="Treatment"),
         fill=guide_legend(title="Treatment"))
dev.off()

 
LMM1 <- lme(viral_skin ~ treatment + AEsum + Adheresum + age + race + gender + 
              + treatment*AEsum + treatment*age, 
            random = ~1| ptid, data = secdata) 
# summary(LMM1)

LMM2 <- lme(viral_blood ~ treatment +
              AEsum + Adheresum + age + race + gender + 
              treatment*Adheresum  + treatment*AEsum, 
            random = ~1| ptid, data = secdata) 
# summary(LMM2)

LMM3 <- lme(Adheresum ~ treatment + age + treatment*age, 
            random = ~1| ptid, data = secdata, method='REML') 
# summary(LMM3)

LMM4 <- lme(Adheresum ~ pill + regimen + age + age*pill + age*regimen, 
            random = ~1| ptid, data = secdata, method='REML') 
# summary(LMM4)

