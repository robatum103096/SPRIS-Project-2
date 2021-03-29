rm(list=ls())

library(tidyverse)
library(dplyr)

# baseline data
bl.data <- read.csv("baseline.csv")[,-1] %>%
  select(ptid, period1, period2, period3, 
         age, race, gender) %>%
  mutate(sequence = ifelse(period1 == "Pill A" & period2 == "Gel B" & period3 == "Gel C", 1, 
                           ifelse(period1 == "Gel C" & period2 == "Pill A" & period3 == "Gel B", 2,
                                  ifelse(period1 == "Gel B" & period2 == "Gel C" & period3 == "Pill A", 3,
                                         ifelse(period1 == "Gel B" & period2 == "Pill A" & period3 == "Gel C", 4,
                                                ifelse(period1 == "Pill A" & period2 == "Gel C" & period3 == "Gel B", 5, 6)))))) %>%
  reshape(direction='long', varying = c("period1", "period2", "period3"), 
          idvar = "ptid",
          timevar = c("period"), 
          times=c("Period 1", "Period 2", "Period 3"), 
          v.names = "treatment")
rownames(bl.data) <- NULL
# pk data
pk.data <- read.csv("baseline.csv")[,-1] %>%
  select(-period1, -period2, -period3, -age, -race, -gender)
# endpoint data
endpt.data <- read.csv("endpoints.csv")[,-(2:4)]


# adverse event pillA
endpt.data.ae.pillA <- endpt.data[,c(1,2:5)] %>%
  reshape(direction='long', 
          varying = c("AE_pillA_week1", "AE_pillA_week2", "AE_pillA_week3", 
                      "AE_pillA_week4"), 
          idvar = "ptid",
          timevar = c("week"), 
          times=c("week1", "week2", "week3", "week4"), 
          v.names = c("AE"), new.row.names = NULL) 
rownames(endpt.data.ae.pillA) <- NULL
endpt.data.ae.pillA$treatment <- rep("Pill A", nrow(endpt.data.ae.pillA))


endpt.data.ae.gelB <- endpt.data[,c(1,6:9)] %>%
  reshape(direction='long', 
          varying = c("AE_gelB_week1", "AE_gelB_week2", "AE_gelB_week3", 
                      "AE_gelB_week4"), 
          idvar = "ptid",
          timevar = c("week"), 
          times=c("week1", "week2", "week3", "week4"), 
          v.names = c("AE"), new.row.names = NULL) 
rownames(endpt.data.ae.gelB) <- NULL
endpt.data.ae.gelB$treatment <- rep("Gel B", nrow(endpt.data.ae.gelB))


endpt.data.ae.gelC <- endpt.data[,c(1,10:13)] %>%
  reshape(direction='long', 
          varying = c("AE_gelC_week1", "AE_gelC_week2", "AE_gelC_week3", 
                      "AE_gelC_week4"), 
          idvar = "ptid",
          timevar = c("week"), 
          times=c("week1", "week2", "week3", "week4"), 
          v.names = c("AE"), new.row.names = NULL) 
rownames(endpt.data.ae.gelC) <- NULL
endpt.data.ae.gelC$treatment <- rep("Gel C", nrow(endpt.data.ae.gelC))

endpt.data.ae.long <- rbind(endpt.data.ae.pillA, endpt.data.ae.gelB, endpt.data.ae.gelC)


endpt.data.ad.pillA <- endpt.data[,c(1,14:17)] %>%
  reshape(direction='long', 
          varying = c("Adhere_pillA_week1", "Adhere_pillA_week2", "Adhere_pillA_week3", 
                      "Adhere_pillA_week4"), 
          idvar = "ptid",
          timevar = c("week"), 
          times=c("week1", "week2", "week3", "week4"), 
          v.names = c("Adhere"), new.row.names = NULL) 
rownames(endpt.data.ad.pillA) <- NULL
endpt.data.ad.pillA$treatment <- rep("Pill A", nrow(endpt.data.ad.pillA))


endpt.data.ad.gelB <- endpt.data[,c(1,18:21)] %>%
  reshape(direction='long', 
          varying = c("Adhere_gelB_week1", "Adhere_gelB_week2", "Adhere_gelB_week3", 
                      "Adhere_gelB_week4"), 
          idvar = "ptid",
          timevar = c("week"), 
          times=c("week1", "week2", "week3", "week4"), 
          v.names = c("Adhere"), new.row.names = NULL) 
rownames(endpt.data.ad.gelB) <- NULL
endpt.data.ad.gelB$treatment <- rep("Gel B", nrow(endpt.data.ad.gelB))


endpt.data.ad.gelC <- endpt.data[,c(1,22:25)] %>%
  reshape(direction='long', 
          varying = c("Adhere_gelC_week1", "Adhere_gelC_week2", "Adhere_gelC_week3", 
                      "Adhere_gelC_week4"), 
          idvar = "ptid",
          timevar = c("week"), 
          times=c("week1", "week2", "week3", "week4"), 
          v.names = c("Adhere"), new.row.names = NULL) 
rownames(endpt.data.ad.gelC) <- NULL
endpt.data.ad.gelC$treatment <- rep("Gel C", nrow(endpt.data.ad.gelC))

endpt.data.ad.long <- rbind(endpt.data.ad.pillA, endpt.data.ad.gelB, endpt.data.ad.gelC)
endpt.data.long <- left_join(endpt.data.ae.long, endpt.data.ad.long, by = c("ptid", "week", "treatment"))

prim.data <- left_join(bl.data[, c("ptid", "sequence", "period", "treatment", "age", "race", "gender")], 
                      endpt.data.long, by = c("ptid", "treatment")) %>%
  select(ptid, sequence, period, treatment, week, AE, Adhere, everything())

write.csv(prim.data, file = "prim_data.csv")


prim.cumweek.data <- prim.data %>%
  group_by(ptid, sequence, period, treatment, age, race, gender) %>%
  summarise(AEsum = sum(AE), Adheresum = sum(Adhere)) %>%
  ungroup()

write.csv(prim.cumweek.data, file = "prim_cumweek_data.csv")


pk.data.blood.long <- rbind(data.frame(ptid = pk.data$ptid, period = rep("Period 1", length(pk.data$ptid)),
           viral_blood = pk.data$bviral1 - pk.data$bviral0),
      data.frame(ptid = pk.data$ptid, period = rep("Period 2", length(pk.data$ptid)),
                 viral_blood = pk.data$bviral3 - pk.data$bviral2),
      data.frame(ptid = pk.data$ptid, period = rep("Period 3", length(pk.data$ptid)),
                 viral_blood = pk.data$bviral5 - pk.data$bviral4))
pk.data.skin.long <- rbind(data.frame(ptid = pk.data$ptid, period = rep("Period 1", length(pk.data$ptid)),
                                       viral_skin = pk.data$sviral1 - pk.data$sviral0),
                            data.frame(ptid = pk.data$ptid, period = rep("Period 2", length(pk.data$ptid)),
                                       viral_skin = pk.data$sviral3 - pk.data$sviral2),
                            data.frame(ptid = pk.data$ptid, period = rep("Period 3", length(pk.data$ptid)),
                                       viral_skin = pk.data$sviral5 - pk.data$sviral4))


sec.data <- left_join(prim.cumweek.data, pk.data.skin.long, by = c("ptid", "period")) %>%
  left_join(pk.data.blood.long, by = c("ptid", "period"))

write.csv(sec.data, file = "secdata.csv")