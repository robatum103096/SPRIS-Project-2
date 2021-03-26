rm(list=ls())

library(tidyverse)
library(dplyr)

# baseline data
bl.data <- read.csv("baseline.csv")[,-1] %>%
  reshape(direction='long', varying = c("period1", "period2", "period3"), 
          idvar = "ptid",
          timevar = c("period"), 
          times=c("period1", "period2", "period3"), 
          v.names = "treatment")
rownames(bl.data) <- NULL
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

all.data <- left_join(bl.data, endpt.data.long, by = c("ptid", "treatment")) %>%
  select(ptid, period, treatment, week, AE, Adhere, everything())

write.csv(all.data, file = "alldata.csv")

pk.data1 <- all.data %>%
  gather(key = viral_time, value = viral_blood, bviral0, bviral1, bviral2,
         bviral3, bviral4, bviral5, bviral6) %>%
  select(-sviral0, -sviral1, -sviral2,
         -sviral3, -sviral4, -sviral5, -sviral6)
pk.data1$viral_time <- factor(pk.data1$viral_time, 
                              levels = c("bviral0", "bviral1", "bviral2", "bviral3", 
                                         "bviral4", "bviral5", "bviral6"))
levels(pk.data1$viral_time) <- 0:6

pk.data2 <- all.data %>%
  gather(key = viral_time, value = viral_skin, sviral0, sviral1, sviral2,
         sviral3, sviral4, sviral5, sviral6) %>%
  select(-bviral0, -bviral1, -bviral2,
         -bviral3, -bviral4, -bviral5, -bviral6)
pk.data2$viral_time <- factor(pk.data2$viral_time, 
                              levels = c("sviral0", "sviral1", "sviral2", "sviral3", 
                                         "sviral4", "sviral5", "sviral6"))
levels(pk.data2$viral_time) <- 0:6

pk.data <- full_join(pk.data1, pk.data2, 
                     by = c("ptid", "period", "treatment", "week", "viral_time",
                            "AE", "Adhere", "age", "race", "gender")) %>%
  gather(key = position, value = viral, viral_skin, viral_blood)
pk.data$position <- factor(pk.data$position, levels = c("viral_blood", "viral_skin"))
levels(pk.data$position) <- c("blood", "skin")

write.csv(pk.data, file = "pkdata.csv")