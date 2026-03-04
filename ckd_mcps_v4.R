library(tidyverse);library(survival);library(nephro);library(flextable);library(gtsummary); library(officer)
library(Epi);library(lubridate);library(ggpubr);library(rms);library(dcurves); library(gt); library(adjustedCurves)
library(survminer); library(cmprsk)

####Data set management####
setwd("~/Mi unidad (obello@facmed.unam.mx)/Mexico City Cohort Study/")
#setwd("~/Downloads/MCPS_2023")

#Load dataset
mcps<-read.csv("Data/2021-004 MCPS BASELINE.csv") %>% left_join(read.csv("Data/2022-012 MCPS MORTALITY.csv"), by="PATID") %>%
  left_join(read.csv("Data/2022-012 MCPS RESURVEY.csv"), by="PATID") %>% 
  left_join(read.csv("Data/2022-012_OmarBelloChavolla_NMR_Second_Release-Corrected/2022-012 MCPS BASE_NMR_DATA_RECALIB-corrected.csv"), by="PATID") %>%
  left_join(read_csv("Data/2022-012_OmarBelloChavolla_NMR_Second_Release-Corrected/2022-012 MCPS RESURVEY_NMR_DATA_RECALIB-corrected.csv"), by="PATID")

#Function to estimate mortality rate
get.death.rate <- function(x){
  a <- x %>% summarise(cases=n(), time=sum(PERSON_YEARS,na.rm = T))
  time <- sum(a$time); cases <- a$cases[2]
  lambda<-(cases/time)*1000; lambda.se<-(sqrt(cases/time^2))*1000
  LI<-lambda-(qnorm(.975)*lambda.se); LS<-lambda+(qnorm(.975)*lambda.se)
  cases.L<-LI*time/1000; cases.U<-LS*time/1000
  data.frame(lambda, cases, LI, LS, time)} #Crude death rates

un.ag.st <- function(x){
  a <- data.frame(); for(i in 1:length(table(x$age_cat2))){
    b <- x %>% dplyr::filter(age_cat2==i) %>% get.death.rate()
    a <- rbind(a, b)}; c <- a %>% cbind(
      "age"=x$age_cat2 %>% levels); c} #Uniformly age and sex standardized

un.ag.st.CI <- function(x){
  x %>% na.omit %>% mutate("w"=1/length(age)) %>%
    mutate("var"=(((w)^2*lambda)/time)) %>% summarise(
      lambda=mean(lambda), var=sum(var), lambda.se=sqrt(sum(var)),
      cases=sum(cases), cases.L=sum(LI), cases.U=sum(LS)) %>% 
    mutate("Lower"=lambda+(cases.L-cases)*sqrt(var/cases),
           "Upper"=lambda+(cases.U-cases)*sqrt(var/cases))} #Standardized confidence intervals

setwd("~/Mi unidad (obello@facmed.unam.mx)/Mexico City Cohort Study/Proyectos/2_Proyectos ENSANUT/CKD_EPI MX")

####Variables####
#Creatinine
mcps$creat<-(mcps$Creatinine.x) / 88.42

#Creatinine resurvey
mcps$creat_resurvey<-(mcps$Creatinine.y) / 88.42

#CKD-EPI creatinine-based 2009 equation
mcps$EFGR.09 <- mcps %>% with(CKDEpi.creat(
  creatinine=creat, sex=MALE, age=AGE, ethnicity = rep(0,nrow(mcps))))

mcps$EGFR_cat.09<-cut(mcps$EFGR.09, breaks=c(-Inf,15,30,45,60,90,Inf),
                      right=F, labels=c("G5","G4","G3b","G3a","G2","G1"))

mcps$CKD.09<-ifelse(mcps$EFGR.09 < 60, 1, 0)
mcps$CKD.09 <- factor(mcps$CKD.09, labels = c("No CKD", "CKD"))
mcps <- mcps %>% 
  mutate(CKD.09.plus = case_when(is.na(mcps$EFGR.09) ~ NA,
                                 EFGR.09 < 60 | mcps$BASE_CKD == 1 ~ "CKD",
                                 EFGR.09 >= 60 & mcps$BASE_CKD == 0 ~ "No CKD"))

mcps$R_Age <- with(mcps, (
  YEAR_RESURVEY + MONTH_RESURVEY/12) - (YEAR_DOB + MONTH_DOB/12)) %>% floor

#CKD-EPI creatinine-based 2009 equation resurvey
mcps$EFGR.09_R <- mcps %>% with(CKDEpi.creat(
  creatinine=creat_resurvey, sex=MALE, age=R_Age, ethnicity = rep(0,nrow(mcps))))

mcps$EGFR_cat.09_R<-cut(mcps$EFGR.09_R, breaks=c(-Inf,15,30,45,60,90,Inf),
                      right=F, labels=c("G5","G4","G3b","G3a","G2","G1"))

mcps$CKD.09_R<-ifelse(mcps$EFGR.09_R < 60, 1, 0)
mcps$CKD.09_R <- factor(mcps$CKD.09_R, labels = c("No CKD", "CKD"))
mcps <- mcps %>% 
  mutate(CKD.09.plus_R = case_when(is.na(mcps$EFGR.09_R) ~ NA,
                                 EFGR.09_R < 60 | mcps$R_HXCKD == 1 ~ "CKD",
                                 EFGR.09_R >= 60 & mcps$R_HXCKD == 0 ~ "No CKD"))

#CKD-EPI 2021 creatinine-based 2021 equation
mcps$EFGR.21 <- mcps %>% with(CKDEpi2021.creat(
  creatinine = creat, sex = MALE, age = AGE ))

mcps$EGFR_cat.21<-cut(mcps$EFGR.21, breaks = c(-Inf,15,30,45,60,90,Inf),
                      right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))

mcps$CKD.21<-ifelse(mcps$EFGR.21 < 60, 1, 0)
mcps$CKD.21 <- factor(mcps$CKD.21, labels = c("No CKD", "CKD"))
mcps <- mcps %>% 
  mutate(CKD.21.plus = case_when(is.na(mcps$EFGR.21) ~ NA,
                                 EFGR.21 < 60 | mcps$BASE_CKD == 1 ~ "CKD",
                                 EFGR.21 >= 60 & mcps$BASE_CKD == 0 ~ "No CKD"))

#CKD-EPI 2021 creatinine-based 2021 equation
mcps$EFGR.21_R <- mcps %>% with(CKDEpi2021.creat(
  creatinine = creat_resurvey, sex = MALE, age = R_Age ))

mcps$EGFR_cat.21_R<-cut(mcps$EFGR.21_R, breaks = c(-Inf,15,30,45,60,90,Inf),
                      right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))

mcps$CKD.21_R<-ifelse(mcps$EFGR.21_R < 60, 1, 0)
mcps$CKD.21_R <- factor(mcps$CKD.21_R, labels = c("No CKD", "CKD"))
mcps <- mcps %>% 
  mutate(CKD.21.plus_R = case_when(is.na(mcps$EFGR.21_R) ~ NA,
                                 EFGR.21_R < 60 | mcps$R_HXCKD == 1 ~ "CKD",
                                 EFGR.21_R >= 60 & mcps$R_HXCKD == 0 ~ "No CKD"))

## Reclassified

mcps$reclassified<-ifelse(((mcps$CKD.09 == "CKD" & mcps$CKD.21 == "No CKD")|(mcps$CKD.09 == "No CKD" & mcps$CKD.21 == "CKD"))==T,"Reclassified","Not reclassified")

#Age categories
mcps$age_cat <- cut(mcps$AGE, breaks = c(1, 35, 45, 55, 65, 75, Inf), right = FALSE)
mcps$age_cat2<-factor(mcps$age_cat, labels = c(1:5))

#Sex
mcps$MALE <- factor(mcps$MALE, labels = c("Female", "Male"))

#Place of residence
mcps$COYOACAN <- factor(mcps$COYOACAN, labels = c("Iztapalapa", "Coyoacan"))

#Educational level
mcps$EDU_LEVEL <- factor(mcps$EDU_LEVEL, labels = c("University", "High school", "Elementary", "Other"))

#Diabetes
mcps$diab <- 0
mcps$diab[mcps$BASE_HBA1C >= 6.5 | mcps$BASE_DIABETES == 1 | 
          mcps$DRUG_D1 == 1 | mcps$DRUG_D2 == 1 | mcps$DRUG_D3 == 1 | mcps$DRUG_D1 == 1] <- 1
mcps$diab <- factor(mcps$diab, labels = c("No diabetes", "Diabetes"))

#Smoking
mcps$smoke[mcps$SMOKEGP == 1] <-  1
mcps$smoke[mcps$SMOKEGP == 2] <-  2
mcps$smoke[mcps$SMOKEGP == 3 | mcps$SMOKEGP == 4 | mcps$SMOKEGP == 5] <-  3
mcps$smoke<-factor(mcps$smoke,labels = c("Never", "Former", "Current"))

#History of cardiovascular disease
mcps$cvd <- 0
mcps$cvd[mcps$BASE_ANGINA == 1 | mcps$BASE_HEARTATTACK == 1 | mcps$BASE_STROKE == 1] <- 1
mcps$cvd <- factor(mcps$cvd, labels = c("No CVD", "CVD"))

#Cardiovascular mortality
mcps$cv_mort <- if_else(mcps$D003 == 1 | mcps$D008 == 1 | mcps$D015 == 1, 1, 0)
mcps$cv_mort2<-mcps$cv_mort
mcps$cv_mort2[mcps$D000==1 & mcps$cv_mort==0]<-2
mcps$cv_mort2 <- factor(mcps$cv_mort2, 0:2, labels=c("censor", "cvd", "death"))

#Kidney mortality
mcps$kd_mort <- if_else(mcps$D016 == 1 | mcps$D017 == 1, 1, 0)
mcps$kd_mort2<-mcps$kd_mort
mcps$kd_mort2[mcps$D000==1 & mcps$kd_mort==0]<-2
mcps$kd_mort2 <- factor(mcps$kd_mort2, 0:2, labels=c("censor", "kdd", "death"))

## Reclassification
mcps <- mcps %>% 
  mutate(reclass_egfr = if_else(EGFR_cat.09 != EGFR_cat.21, "Upward reclassification", "No reclassification"),
         reclass_ckd = if_else(CKD.09 == CKD.21,"No reclassification","CKD reclassification")) %>%
  mutate(reclass_ckd=relevel(factor(reclass_ckd, ordered = F), ref="No reclassification"),
         reclass_egfr_ckd = case_when(reclass_egfr=="Upward reclassification" & reclass_ckd=="No reclassification"~"eGFR reclassification",
                                      reclass_egfr=="Upward reclassification" & reclass_ckd=="CKD reclassification"~"CKD reclassification",
                                      reclass_egfr=="No reclassification" & reclass_ckd=="No reclassification"~"No reclassification",
                                      reclass_egfr=="No reclassification" & reclass_ckd=="CKD reclassification"~"CKD reclassification")) %>%
  mutate(reclass_ckd=relevel(factor(reclass_ckd, ordered = F), ref="No reclassification"), reclass_egfr_ckd=relevel(factor(reclass_egfr_ckd, ordered = F), ref="No reclassification"))
levels(mcps$reclass_egfr_ckd)<-c(levels(mcps$reclass_egfr_ckd), "CKD (No reclassification)")
mcps$reclass_egfr_ckd[mcps$CKD.09 == "CKD" & mcps$CKD.21 =="CKD"]<-"CKD (No reclassification)"

####Filter dataset####
#Filter 1: missing data on mortality or creatinine
mcps1 <- mcps %>% 
  filter(creat > 0) %>% 
  select(PATID,PERSON_YEARS, D000, D003, D008, D015, cv_mort,cv_mort2,kd_mort,kd_mort2, EFGR.09, EFGR.21, CKD.09, CKD.21, CKD.09.plus, CKD.21.plus, EGFR_cat.09, EGFR_cat.21, 
         AGE, age_cat,age_cat2, MALE, COYOACAN, EDU_LEVEL, diab, smoke, cvd, SBP1, BASE_DIABETES, BASE_HYPERTENSION, BASE_CKD, reclass_egfr,reclass_ckd,reclass_egfr_ckd,
         BASE_HBA1C, Total_C.x, EFGR.09_R, EFGR.21_R, CKD.09_R, CKD.21_R, CKD.09.plus_R, CKD.21.plus_R, EGFR_cat.09_R, EGFR_cat.21_R,YEARS_BASELINE_RESURVEY) %>% 
  filter(if_all(c(PERSON_YEARS:Total_C.x), ~ !is.na(.x)))

#Population flowchart
#Total 
nrow(mcps) 
#Unmeasured or implausible creatinine
nrow(mcps)-nrow(mcps %>% filter(!is.na(creat)))
nrow(mcps)-nrow(mcps %>% filter(creat > 0)) 
# Total 2
nrow(mcps %>% filter(creat > 0)) 
# Unmeasured creatinine
nrow(mcps %>% filter(creat > 0))-nrow(mcps %>% filter(creat > 0) %>% filter(!is.na(D000) | !is.na(PERSON_YEARS))) 
# Total 3
t1<-nrow(mcps %>% filter(creat > 0) %>% filter(!is.na(D000) | !is.na(PERSON_YEARS)));t1
# Missing data
t1-nrow(mcps1)
#Complete data
nrow(mcps1) 

#Rename datasets
mcps_fin <- mcps1 
remove(mcps1)

#### Missing data analysis ####
#Characteristics of those reclassified with and without CKD
mcps$missing<-ifelse(mcps$PATID %in% mcps_fin$PATID,0,1)

mcps %>% 
  select(missing, EFGR.09,EFGR.21,AGE, MALE, COYOACAN, EDU_LEVEL, D000, cv_mort,kd_mort,diab, BASE_HYPERTENSION, BASE_CKD, cvd,BASE_HBA1C) %>% 
  tbl_summary(by = missing,
              statistic = list(c(AGE) ~ "{mean} (±{sd})",
                               c(EFGR.09, EFGR.21, BASE_HBA1C) ~ "{median} ({p25}, {p75})"),
              label = list(AGE ~ "Age",
                           MALE ~ "Men",
                           EFGR.09 ~ "eGFR with CKD-EPI 2009",
                           EFGR.09 ~ "eGFR with CKD-EPI 2021",
                           D000 ~ "All-cause mortality",
                           cv_mort ~ "Cardiovascular mortality",
                           kd_mort ~ "Kidney disease mortality",
                           diab ~ "Dibetes",
                           BASE_HYPERTENSION ~ "Hypertension",
                           BASE_CKD ~ "CKD",
                           cvd ~ "Cardiovascular disease",
                           BASE_HBA1C ~ "HbA1C")) %>% 
  add_p(pvalue_fun = label_style_pvalue(digits = 2)) %>% 
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))%>%
  add_overall() %>%
  as_gt() %>% 
  gt::gtsave("Tables/missing_tab.docx")

####eGFR categories and CKD prevalence#### 
#eGFR categories and CKD based on eGFR
prev_09 <- mcps_fin %>% group_by(EGFR_cat.09) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(sum = sum(n),
         pct = round((n / sum), 4),
         equation = "2009 equation",
         label_09 = paste0(n, " (", sprintf(pct*100, fmt = "%#.2f"), "%)")) %>% 
  rename(egfr = EGFR_cat.09) %>% 
  arrange(desc(n))

prev_ckd_09 <- mcps_fin %>% group_by(CKD.09) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(sum = sum(n),
         pct = round((n / sum), 4),
         equation = "2009 equation",
         label_09 = paste0(n, " (", sprintf(pct*100, fmt = "%#.2f"), "%)")) %>% 
  rename(ckd = CKD.09) %>% 
  filter(ckd == "CKD")

prev_21 <- mcps_fin %>% group_by(EGFR_cat.21) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(sum = sum(n),
         pct = round((n / sum), 4),
         equation = "2021 equation",
         label_21 = paste0(n, " (", sprintf(pct*100, fmt = "%#.2f"), "%)")) %>% 
  rename(egfr = EGFR_cat.21) %>% 
  arrange(desc(n))

prev_ckd_21 <- mcps_fin %>% group_by(CKD.21) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(sum = sum(n),
         pct = round((n / sum), 4),
         label_21 = paste0(n, " (", sprintf(pct*100, fmt = "%#.2f"), "%)")) %>% 
  rename(ckd = CKD.21) %>% 
  filter(ckd == "CKD")

egfr_prev <- cbind(prev_09, prev_21)
egfr_prev %>% select(1, label_09, label_21) %>% 
  gt() %>% 
  gtsave(filename = "Tables/tabla_egfr_mcps.docx")

cbind(prev_ckd_09, prev_ckd_21) %>% 
  select(1, label_09, label_21) %>% 
  gt() %>% 
  gtsave(filename = "Tables/tabla_ckd_mcps.docx")
  
####Previous diagnosis and CKD prevalence####
#Previous dignosis of CKD
ckd_prev <- mcps_fin %>% group_by(BASE_CKD) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(sum = sum(n),
         pct = round((n / sum), 4),
         label = paste0(n, " (", sprintf(pct*100, fmt = "%#.2f"), "%)"),
         group = "Previous diagnosis") %>% 
  filter(BASE_CKD == 1) %>% 
  select(5, 6)

#CKD eGFR 2009 + previous diagnosis 
ckd_plus_09 <- mcps_fin %>% group_by(CKD.09.plus) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(sum = sum(n),
         pct = round((n / sum), 4),
         label = paste0(n, " (", sprintf(pct*100, fmt = "%#.2f"), "%)"),
         group = "Previous diagnosis + eGFR 2009") %>% 
  filter(CKD.09.plus == "CKD") %>% 
  select(5, 6)

#CKD eGFR 2021 + previous diagnosis
ckd_plus_21 <- mcps_fin %>% group_by(CKD.21.plus) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(sum = sum(n),
         pct = round((n / sum), 4),
         label = paste0(n, " (", sprintf(pct*100, fmt = "%#.2f"), "%)"),
         group = "Previous diagnosis + eGFR 2021") %>% 
  filter(CKD.21.plus == "CKD") %>% 
  select(5, 6)

#Total table
rbind(ckd_prev, ckd_plus_09, ckd_plus_21) %>% 
  gt() %>% 
  gtsave(filename = "Tables/tabla_ckd_prev_mcps.docx")

####eGFR distribution between two equations####
summary(mcps_fin$EFGR.09)
summary(mcps_fin$EFGR.21)

#Distribution between two equations
fig0_1 <- mcps_fin %>% ggplot(aes()) + 
  geom_density(aes(x = EFGR.09, color = "CKD-EPI 2009"), show.legend = FALSE) +
  geom_density(aes(x = EFGR.21, color = "CKD-EPI 2021"), show.legend = FALSE) +
  labs(x = "eGFR",
       y = "Density",
       color = "") +
  geom_vline(xintercept = c(15,30,45,60,90), color = "gray", linetype = "dashed") +
  scale_x_continuous(breaks = c(15,30,45,60, 90)) +
  scale_color_manual(values = c("#3172B1", "#AE4532")) +
  theme_classic()

fig0_2 <- mcps_fin %>% filter(CKD.09 == 1 | CKD.21 == 1) %>% ggplot(aes()) + 
  geom_density(aes(x = EFGR.09, color = "CKD-EPI 2009")) +
  geom_density(aes(x = EFGR.21, color = "CKD-EPI 2021")) +
  labs(title = "eGFR distribution among participants with CKD",
       x = "eGFR",
       y = "Proportion",
       color = "CKD-EPI equation") +
  scale_color_manual(values = c("indianred2", "#009ACD")) +
  theme_classic()

####Reclassification in every eGFR category####
tab_pre <- mcps_fin %>% group_by(fct_rev(EGFR_cat.09), fct_rev(EGFR_cat.21)) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(sum = sum(n),
         pct = (n / sum)*100) %>% 
  data.frame()

fig0_3 <- tab_pre %>% filter(pct < 45) %>% 
  ggplot(aes(x = fct_rev.EGFR_cat.09., y = pct, fill = "Upward reclassification")) +
  geom_col() +
  labs(title = "Reclassification according to CKD-EPI 2021",
       x = "CKD-EPI 2009 eGFR",
       y = "Percentage reclassification",
       fill = "") +
  ylim(0, 45) +
  theme_classic()

table <- tribble(
  ~"G1",                      ~"G2",                        ~"G3a",                       ~"G3b",                       ~"G4",                        ~"G5",           ~"Total", 
  tab_pre[1,3],                 0,                            0,                            0,                            0,                            0,               tab_pre[1,3],
  tab_pre[2,3],                 tab_pre[3,3],                 0,                            0,                            0,                            0,               tab_pre[2,3] + tab_pre[3,3],
  0,                            tab_pre[4,3],                 tab_pre[5,3],                 0,                            0,                            0,               tab_pre[4,3] + tab_pre[5,3],
  0,                            0,                            tab_pre[6,3],                 tab_pre[7,3],                 0,                            0,               tab_pre[6,3] + tab_pre[7,3],
  0,                            0,                            0,                            tab_pre[8,3],                 tab_pre[9,3],                 0,               tab_pre[8,3] + tab_pre[9,3],
  0,                            0,                            0,                            0,                            tab_pre[10,3],                tab_pre[11,3],   tab_pre[10,3] + tab_pre[11,3],
  tab_pre[1,3] + tab_pre[2,3],  tab_pre[3,3] + tab_pre[4,3],  tab_pre[5,3] + tab_pre[6,3],  tab_pre[7,3] + tab_pre[8,3],  tab_pre[9,3] + tab_pre[10,3], tab_pre[11,3],   0
)

####Reclassification for CKD category####
tab_pre_2 <- mcps_fin %>% group_by(CKD.09, CKD.21) %>% 
  dplyr::summarize(n = n()) %>% 
  mutate(sum = sum(n),
         pct = (n / sum)*100) %>% 
  data.frame()


fig0_4 <- tab_pre_2 %>% filter(CKD.09 == 1 & CKD.21 == 0) %>% 
  ggplot(aes(x = CKD.09, y = pct, fill = "Upward reclassification")) +
  geom_col(width = 0.5) +
  labs(title = "CKD reclassification according to CKD-EPI 2021",
       x = "CKD-EPI 2009 CKD",
       y = "Percentage reclassification",
       fill = "") +
  ylim(0, 45) +
  theme_classic()

#Characteristics of those reclassified with and without CKD
mcps_fin %>% mutate(ckd_reclass = case_when(CKD.09 == "No CKD" & CKD.21 == "No CKD" ~ "No CKD (no reclassification)",
                                            CKD.09 == "CKD" & CKD.21 == "No CKD" ~ "Upward reclassification",
                                            CKD.09 == "CKD" & CKD.21 == "CKD" ~ "With CKD (no reclassification)")) %>% 
  select(ckd_reclass, AGE, MALE, EFGR.09, EFGR.21, D000, BASE_DIABETES, BASE_HYPERTENSION, BASE_CKD, BASE_HBA1C) %>% 
  tbl_summary(by = ckd_reclass,
              statistic = list(c(AGE) ~ "{mean} (±{sd})",
                               c(EFGR.09, EFGR.21, BASE_HBA1C) ~ "{median} ({p25}, {p75})"),
              label = list(AGE ~ "Age",
                           MALE ~ "Men",
                           EFGR.09 ~ "eGFR 2009",
                           EFGR.21 ~ "eGFR 2021",
                           D000 ~ "Events",
                           BASE_DIABETES ~ "Dibetes",
                           BASE_HYPERTENSION ~ "Hypertension",
                           BASE_CKD ~ "CKD",
                           BASE_HBA1C ~ "HbA1C"),
              missing = "no") %>% 
  add_p() %>% 
  as_gt() %>% 
  gt::gtsave("tab.docx")

#Characteristics of those reclassified with and without CKD
mcps_fin %>% 
  select(reclass_egfr_ckd, AGE, MALE, EFGR.09, EFGR.21, D000, BASE_DIABETES, BASE_HYPERTENSION, BASE_CKD, BASE_HBA1C) %>% 
  tbl_summary(by = reclass_egfr_ckd,
              statistic = list(c(AGE) ~ "{mean} (±{sd})",
                               c(EFGR.09, EFGR.21, BASE_HBA1C) ~ "{median} ({p25}, {p75})"),
              label = list(AGE ~ "Age",
                           MALE ~ "Men",
                           EFGR.09 ~ "eGFR 2009",
                           EFGR.21 ~ "eGFR 2021",
                           D000 ~ "Events",
                           BASE_DIABETES ~ "Dibetes",
                           BASE_HYPERTENSION ~ "Hypertension",
                           BASE_CKD ~ "CKD",
                           BASE_HBA1C ~ "HbA1C"),
              missing = "no") %>% 
  add_p() %>% 
  as_gt() %>% 
  gt::gtsave("tab_2.docx")

####Non-linearity of eGFR####
##eGFR estimated with 2009 CKD-EPI equation
#Non-linearity check with restricted cubic splines
m_fin <- datadist(mcps_fin)
options(datadist = "m_fin")

fit_egfr09_rcs <- cph(Surv(PERSON_YEARS, D000) ~ rcs(EFGR.09), data = mcps_fin,
                  x = TRUE, y = TRUE, surv = TRUE) 
plot(rms::Predict(fit_egfr09_rcs))

#Non-linearity check with penalized smoothing splines
fit_egfr09_psp <- coxph(Surv(PERSON_YEARS, D000) ~ pspline(EFGR.09, df = 4), data=mcps_fin) 

termplot(fit_egfr09_psp, term = "pspline(EFGR.09, df = 4)", se = TRUE, col.term = 1, col.se = 1)

##eGFR estimated with 2021 CKD-EPI equation
#Non-linearity check with restricted cubic splines
fit_egfr21_rcs <- cph(Surv(PERSON_YEARS, D000) ~ rcs(EFGR.21), data = mcps_fin,
                      x = TRUE, y = TRUE, surv = TRUE) 
plot(rms::Predict(fit_egfr21_rcs))

#Non-linearity check with penalized smoothing splines
fit_egfr21_psp <- coxph(Surv(PERSON_YEARS, D000) ~ pspline(EFGR.21, df = 4), data=mcps_fin) 

termplot(fit_egfr21_psp, term = "pspline(EFGR.21, df = 4)", se = TRUE, col.term = 1, col.se = 1)

#### --------- All-cause mortality analyses ---------####
#####Administrative censoring####
mcps_pre <- survSplit(Surv(PERSON_YEARS, D000) ~ ., data = mcps_fin, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

mcps_pre <- survSplit(Surv(PERSON_YEARS, D000) ~ ., data = mcps_fin, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)

#####Baseline Cox proportional hazards models####
#Five-year mortality risk
base5 <- coxph(Surv(PERSON_YEARS, D000) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
               x = TRUE, y = TRUE, data = mcps5)
cox.zph(base5)
cph_base5 <- cph(Surv(PERSON_YEARS, D000) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x,
                 x = TRUE, y = TRUE, surv = TRUE, data = mcps5)

#Ten-year mortality risk
base10 <- coxph(Surv(PERSON_YEARS, D000) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                x = TRUE, y = TRUE, data = mcps10)
cox.zph(base10)
cph_base10 <- cph(Surv(PERSON_YEARS, D000) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x,
                 x = TRUE, y = TRUE, surv = TRUE, data = mcps10)

#Number of events
events5 <- mcps5 %>% filter(D000 == 1) %>% count(MALE, D000)
events10 <- mcps10 %>% filter(D000 == 1) %>% count(MALE, D000)

#####5-year predictive mortality risk - Addition of eGFR estimated by 2009 CKD-EPI equation####
egfr09_5 <- update(base5, . ~ . + rcs(EFGR.09))
cox.zph(egfr09_5) %>% plot()

t1<-egfr09_5 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_5_09.docx")

anova(base5, egfr09_5)

#Model fitted using rms::cph()
m5 <- datadist(mcps5)
options(datadist = "m5")

cph_egfr09_5 <- cph(Surv(PERSON_YEARS, D000) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.09), data = mcps5,
                x = TRUE, y = TRUE, surv = TRUE)
?rcs
specs(cph_egfr09_5) #specs() function gives knots location

#Graph
Predict(cph_egfr09_5) %>% 
  ggplot()

#####5-year predictive mortality risk - Addition of eGFR estimated by 2021 CKD-EPI equation####
egfr21_5 <- update(base5, . ~ . + rcs(EFGR.21))
cox.zph(egfr21_5) %>% plot()
t1<-egfr21_5 %>%
  tbl_regression(exp=T) %>% 
  add_p(pvalue_fun = label_style_pvalue(digits = 2)) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_5_21.docx")

anova(base5, egfr21_5)

#Model fitted using rms::cph()
cph_egfr21_5 <- cph(Surv(PERSON_YEARS, D000) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.21), data = mcps5,
                x = TRUE, y = TRUE, surv = TRUE)

specs(cph_egfr21_5) #specs() function gives knots location

#Graph
Predict(cph_egfr21_5) %>% 
  ggplot()

#####5-year assessment of discrimination measures####
#Time range
c_base5 <- concordance(base5, keepstrata = TRUE)
c_egfr09_5 <- concordance(egfr09_5, keepstrata = TRUE)
c_egfr21_5 <- concordance(egfr21_5, keepstrata = TRUE)

#Fixed point for eGFR 2009
times <- seq(0.5, 5, by = 0.5)
data09 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base5, ymax = times[i])
  c_con <- concordance(egfr09_5, ymax = times[i])
  #Data
  data09[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_09_5 <- do.call(rbind, data09)

#Fixed point for eGFR 2021
times <- seq(0.5, 5, by = 0.5)
data21 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base5, ymax = times[i])
  c_con <- concordance(egfr21_5, ymax = times[i])
  #Data
  data21[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_21_5 <- do.call(rbind, data21)
cstat_5_final <- rbind(cstat_09_5, cstat_21_5)
cstat_5_final$equation <- c(rep("2009 equation", 20), rep("2021 equation", 20))
cstat_5_final$year <- c(rep("5-year", 40))
cstat_5_final$lower <- cstat_5_final$coef - 1.96 * cstat_5_final$se
cstat_5_final$upper <- cstat_5_final$coef + 1.96 * cstat_5_final$se

#####5-year Brier and scaled Brier score####
score_mcps5 <- riskRegression::Score(list("base" = base5,
                                          "egfr_09" = egfr09_5,
                                          "egfr_21" = egfr21_5), 
                                     formula = Surv(PERSON_YEARS, D000) ~ 1, 
                                     data = mcps5, 
                                     conf.int = TRUE, 
                                     times = 4.99,
                                     cens.model = "km", 
                                     metrics = "brier",
                                     summary = "ipa")

score_5 <- score_mcps5[["Brier"]][["score"]]
score_5 <- score_5 %>% select(1,3,5,6,7)

#####10-year predictive mortality risk - Addition of eGFR estimated by 2009 CKD-EPI equation####
egfr09_10 <- update(base10, . ~ . + rcs(EFGR.09))
cox.zph(egfr09_10) %>% plot()
t1<-egfr09_10 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_10_09.docx")
anova(base10, egfr09_10)

#Model fitted using rms::cph()
m10 <- datadist(mcps10)
options(datadist = "m10")

cph_egfr09_10 <- cph(Surv(PERSON_YEARS, D000) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.09), data = mcps10,
                    x = TRUE, y = TRUE, surv = TRUE)

specs(cph_egfr09_10) #specs() function gives knots location

#Graph
Predict(cph_egfr09_10) %>% 
  ggplot()

#####10-year predictive mortality risk - Addition of eGFR estimated by 2021 CKD-EPI equation####
egfr21_10 <- update(base10, . ~ . + rcs(EFGR.21))
cox.zph(egfr21_10) %>% plot()
t1<-egfr21_10 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_10_21.docx")

anova(base10, egfr09_10)

anova(base10, egfr21_10)

#Model fitted using rms::cph()
cph_egfr21_10 <- cph(Surv(PERSON_YEARS, D000) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.21), data = mcps10,
                    x = TRUE, y = TRUE, surv = TRUE)

specs(cph_egfr21_10) #specs() function gives knots location

#Graph
Predict(cph_egfr21_10) %>% 
  ggplot()

#####10-year assessment of discrimination measures####
#Time range
c_base10 <- concordance(base10, keepstrata = TRUE)
c_egfr09_10 <- concordance(egfr09_10, keepstrata = TRUE)
c_egfr21_10 <- concordance(egfr21_10, keepstrata = TRUE)

#Fixed point for eGFR 2009
times <- seq(6, 10, by = 0.5)
data09 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base10, ymax = times[i])
  c_con <- concordance(egfr09_10, ymax = times[i])
  #Data
  data09[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_09_10 <- do.call(rbind, data09)

#Fixed point for eGFR 2021
times <- seq(6, 10, by = 0.5)
data21 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base10, ymax = times[i])
  c_con <- concordance(egfr21_10, ymax = times[i])
  #Data
  data21[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_21_10 <- do.call(rbind, data21)
cstat_10_final <- rbind(cstat_09_10, cstat_21_10)
cstat_10_final$equation <- c(rep("2009 equation", 18), rep("2021 equation", 18))
cstat_10_final$year <- c(rep("10-year", 36))
cstat_10_final$lower <- cstat_10_final$coef - 1.96 * cstat_10_final$se
cstat_10_final$upper <- cstat_10_final$coef + 1.96 * cstat_10_final$se


#####10-year Brier and scaled Brier score####
score_mcps10 <- riskRegression::Score(list("base" = base10,
                                          "egfr_09" = egfr09_10,
                                          "egfr_21" = egfr21_10), 
                                     formula = Surv(PERSON_YEARS, D000) ~ 1, 
                                     data = mcps10, 
                                     conf.int = TRUE, 
                                     times = 10,
                                     cens.model = "km", 
                                     metrics = "brier",
                                     summary = "ipa")

score_10 <- score_mcps10[["Brier"]][["score"]]
score_10 <- score_10 %>% select(1,3,5,6,7)

#####Discrimination assessment over whole follow-up period####
cstat_10_final <- rbind(cstat_09_10, cstat_21_10)
cstat_10_final$equation <- c(rep("2009 equation", 18), rep("2021 equation", 18))
cstat_10_final$year <- c(rep("10-year", 36))
cstat_10_final$lower <- cstat_10_final$coef - 1.96 * cstat_10_final$se
cstat_10_final$upper <- cstat_10_final$coef + 1.96 * cstat_10_final$se

cstat_final <- rbind(cstat_5_final, cstat_10_final)

#C-statistic difference during the first 5 years of follow-up
cstat_diff_09 <- cstat_final %>% 
  filter(type != "Base model", time <= 5, equation == "2009 equation") %>% 
  rename(coef_09 = coef)

cstat_diff_21 <- cstat_final %>% 
  filter(type != "Base model", time <= 5, equation == "2021 equation") %>% 
  rename(coef_21 = coef)

cstat_diff <- cbind(cstat_diff_09, cstat_diff_21)
cstat_diff$diff <- cstat_diff$coef_09 - cstat_diff$coef_21

####Discrimination and performance table####
tab_dis_perf <- tibble("Year" = c("Baseline", "2009", "2021"),
                       "Discrimination" = c(paste0(round(cstat_final[75,2], 3), " (95%CI ", 
                                                   round(cstat_final[75,7], 3), "-",
                                                   round(cstat_final[75,8], 3), ")"),
                                            paste0(round(cstat_final[58,2], 3), " (95%CI ", 
                                                   round(cstat_final[58,7], 3), "-",
                                                   round(cstat_final[58,8], 3), ")"),
                                            paste0(round(cstat_final[76,2], 3), " (95%CI ", 
                                                   round(cstat_final[76,7], 3), "-",
                                                   round(cstat_final[76,8], 3), ")")),
                       "Brier score (5 year)" = c(paste0(round(score_5[2,2], 3), " (95%CI ", 
                                                         round(score_5[2,3], 3), "-",
                                                         round(score_5[2,4], 3), ")"),
                                                  paste0(round(score_5[3,2], 3), " (95%CI ", 
                                                         round(score_5[3,3], 3), "-",
                                                         round(score_5[3,4], 3), ")"),
                                                  paste0(round(score_5[4,2], 3), " (95%CI ", 
                                                         round(score_5[4,3], 3), "-",
                                                         round(score_5[4,4], 3), ")")),
                       "Scaled Brier score (5 year)" = c(round(score_5[2,5], 4)*100, 
                                                         round(score_5[3,5], 4)*100, 
                                                         round(score_5[4,5], 4)*100),
                       "Brier score (10 year)" = c(paste0(round(score_10[2,2], 3), " (95%CI ", 
                                                          round(score_10[2,3], 3), "-",
                                                          round(score_10[2,4], 3), ")"),
                                                   paste0(round(score_10[3,2], 3), " (95%CI ", 
                                                          round(score_10[3,3], 3), "-",
                                                          round(score_10[3,4], 3), ")"),
                                                   paste0(round(score_10[4,2], 3), " (95%CI ", 
                                                          round(score_10[4,3], 3), "-",
                                                          round(score_10[4,4], 3), ")")),
                       "Scaled Brier score (10 year)" = c(round(score_10[2,5], 4)*100, 
                                                          round(score_10[3,5], 4)*100, 
                                                          round(score_10[4,5], 4)*100))

tab_dis_perf %>% 
  gt() %>% 
  gtsave(filename = "Tables/tabla_discriminacion_performance_overall.docx")

#### Risk comparison ####
#2009 Equation
partial_09 <- termplot(egfr09_10, se = TRUE, plot = FALSE)
egfr_term_09 <- partial_09$EFGR.09[1:130610,]
center_09 <- -0.06258325 #Y en aprox 95 eGFR
ytemp_09 <- egfr_term_09$y + outer(egfr_term_09$se, c(0, -1.96, 1.96), '*')

#2021 Equation
partial_21 <- termplot(egfr21_10, se = TRUE, plot = FALSE)
egfr_term_21 <- partial_21$EFGR.21[1:130610,]
center_21 <- 0.007670158 #Y en aprox 95 eGFR
ytemp_21 <- egfr_term_21$y + outer(egfr_term_21$se, c(0, -1.96, 1.96), '*')

hr_pre09 <- as.data.frame(exp(ytemp_09 - center_09))
hr_pre21 <- as.data.frame(exp(ytemp_21 - center_21))
df <- data.frame(x_09 = egfr_term_09$x,
                 y_09 = hr_pre09$V1,
                 lower_09 = hr_pre09$V2,
                 upper_09 = hr_pre09$V3,
                 x_21 = egfr_term_21$x,
                 y_21 = hr_pre21$V1,
                 lower_21 = hr_pre21$V2,
                 upper_21 = hr_pre21$V3)

#Script para figura en ckd_prev.R


#### Chi-squared and log-likelihood ####
summary(base5)$logtest
summary(egfr09_5)$logtest
summary(egfr21_5)$logtest

anova(base5,egfr09_5)
anova(base5,egfr21_5)
anova(egfr09_5,egfr21_5)

summary(base10)$logtest
summary(egfr09_10)$logtest
summary(egfr21_10)$logtest

anova(base10,egfr09_10)
anova(base10,egfr21_10)
anova(egfr09_10,egfr21_10)

#### --------- Cardiovascular mortality analyses ---------####
#####Administrative censoring####
fg_data <- finegray(Surv(PERSON_YEARS, cv_mort2) ~ ., data = mcps_fin)

mcps_pre <- survSplit(Surv(fgstart, fgstop, fgstatus) ~ ., data = fg_data, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

mcps_pre <- survSplit(Surv(fgstart, fgstop, fgstatus) ~ ., data = fg_data, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)

#####Baseline Cox proportional hazards models####
#Five-year mortality risk
base5 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
               x = TRUE, y = TRUE, weight=fgwt,data = mcps5)

cph_base5 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x,
                 x = TRUE, y = TRUE, surv = TRUE, weight=fgwt,data = mcps5)

#Ten-year mortality risk
base10 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                x = TRUE, y = TRUE, weight=fgwt,data = mcps10)

cph_base10 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x,
                  x = TRUE, y = TRUE, surv = TRUE, weight=fgwt,data = mcps10)

#####5-year predictive mortality risk - Addition of eGFR estimated by 2009 CKD-EPI equation####
egfr09_5 <- update(base5, . ~ . + rcs(EFGR.09))
cox.zph(egfr09_5) %>% plot()

t1<-egfr09_5 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_05_09_cv.docx")
anova(base5, egfr09_5)

#Model fitted using rms::cph()
m5 <- datadist(mcps5)
options(datadist = "m5")

cph_egfr09_5 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.09), data = mcps5,
                    x = TRUE, y = TRUE, surv = TRUE,weight=fgwt)

specs(cph_egfr09_5) #specs() function gives knots location

#Graph
Predict(cph_egfr09_5) %>% 
  ggplot()

#####5-year predictive mortality risk - Addition of eGFR estimated by 2021 CKD-EPI equation####
egfr21_5 <- update(base5, . ~ . + rcs(EFGR.21))
cox.zph(egfr21_5) %>% plot()

t1<-egfr21_5 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_05_21_cv.docx")

anova(base5, egfr21_5)

#Model fitted using rms::cph()
cph_egfr21_5 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.21), data = mcps5,
                    x = TRUE, y = TRUE, surv = TRUE,weight=fgwt)

specs(cph_egfr21_5) #specs() function gives knots location

#Graph
Predict(cph_egfr21_5) %>% 
  ggplot()

#####5-year assessment of discrimination measures####
#Time range
c_base5 <- concordance(base5, keepstrata = TRUE)
c_egfr09_5 <- concordance(egfr09_5, keepstrata = TRUE)
c_egfr21_5 <- concordance(egfr21_5, keepstrata = TRUE)

#Fixed point for eGFR 2009
times <- seq(0.5, 5, by = 0.5)
data09 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base5, ymax = times[i])
  c_con <- concordance(egfr09_5, ymax = times[i])
  #Data
  data09[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_09_5 <- do.call(rbind, data09)

#Fixed point for eGFR 2021
times <- seq(0.5, 5, by = 0.5)
data21 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base5, ymax = times[i])
  c_con <- concordance(egfr21_5, ymax = times[i])
  #Data
  data21[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_21_5 <- do.call(rbind, data21)
cstat_5_final <- rbind(cstat_09_5, cstat_21_5)
cstat_5_final$equation <- c(rep("2009 equation", 20), rep("2021 equation", 20))
cstat_5_final$year <- c(rep("5-year", 40))
cstat_5_final$lower <- cstat_5_final$coef - 1.96 * cstat_5_final$se
cstat_5_final$upper <- cstat_5_final$coef + 1.96 * cstat_5_final$se

#####5-year Brier and scaled Brier score####
mcps_pre <- survSplit(Surv(PERSON_YEARS, cv_mort) ~ ., data = mcps_fin, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

base5 <- coxph(Surv(PERSON_YEARS, cv_mort) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
               x = TRUE, y = TRUE, data = mcps5)
egfr09_5 <- update(base5, . ~ . + rcs(EFGR.09))
egfr21_5 <- update(base5, . ~ . + rcs(EFGR.21))

score_mcps5 <- riskRegression::Score(list("base" = base5,
                                          "egfr_09" = egfr09_5,
                                          "egfr_21" = egfr21_5), 
                                     formula = Surv(PERSON_YEARS, cv_mort) ~ 1, 
                                     data = mcps5, 
                                     conf.int = TRUE, 
                                     times = 4.99,
                                     cens.model = "km", 
                                     metrics = "brier",
                                     summary = "ipa")

score_5 <- score_mcps5[["Brier"]][["score"]]
score_5 <- score_5 %>% select(1,3,5,6,7)

#####10-year predictive mortality risk - Addition of eGFR estimated by 2009 CKD-EPI equation####
egfr09_10 <- update(base10, . ~ . + rcs(EFGR.09))
cox.zph(egfr09_10) %>% plot()

t1<-egfr09_10 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_10_09_cv.docx")


anova(base10, egfr09_10)

#Model fitted using rms::cph()
m10 <- datadist(mcps10)
options(datadist = "m10")

cph_egfr09_10 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.09), data = mcps10,
                     x = TRUE, y = TRUE, surv = TRUE,weight=fgwt)

specs(cph_egfr09_10) #specs() function gives knots location

#Graph
Predict(cph_egfr09_10) %>% 
  ggplot()

#####10-year predictive mortality risk - Addition of eGFR estimated by 2021 CKD-EPI equation####
egfr21_10 <- update(base10, . ~ . + rcs(EFGR.21))
cox.zph(egfr21_10) %>% plot()

t1<-egfr21_10 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_10_21_cv.docx")

anova(base10, egfr21_10)

#Model fitted using rms::cph()
cph_egfr21_10 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.21), data = mcps10,
                     x = TRUE, y = TRUE, surv = TRUE,weight=fgwt)

specs(cph_egfr21_10) #specs() function gives knots location

#Graph
Predict(cph_egfr21_10) %>% 
  ggplot()

#####10-year assessment of discrimination measures####
#Time range
c_base10 <- concordance(base10, keepstrata = TRUE)
c_egfr09_10 <- concordance(egfr09_10, keepstrata = TRUE)
c_egfr21_10 <- concordance(egfr21_10, keepstrata = TRUE)

#Fixed point for eGFR 2009
times <- seq(6, 10, by = 0.5)
data09 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base10, ymax = times[i])
  c_con <- concordance(egfr09_10, ymax = times[i])
  #Data
  data09[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_09_10 <- do.call(rbind, data09)

#Fixed point for eGFR 2021
times <- seq(6, 10, by = 0.5)
data21 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base10, ymax = times[i])
  c_con <- concordance(egfr21_10, ymax = times[i])
  #Data
  data21[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_21_10 <- do.call(rbind, data21)
cstat_10_final <- rbind(cstat_09_10, cstat_21_10)
cstat_10_final$equation <- c(rep("2009 equation", 18), rep("2021 equation", 18))
cstat_10_final$year <- c(rep("10-year", 36))
cstat_10_final$lower <- cstat_10_final$coef - 1.96 * cstat_10_final$se
cstat_10_final$upper <- cstat_10_final$coef + 1.96 * cstat_10_final$se

cstat_final <- rbind(cstat_5_final, cstat_10_final)

#####10-year Brier and scaled Brier score####
mcps_pre <- survSplit(Surv(PERSON_YEARS, cv_mort) ~ ., data = mcps_fin, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)
base10 <- coxph(Surv(PERSON_YEARS, cv_mort) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                x = TRUE, y = TRUE, data = mcps10)
egfr09_10 <- update(base10, . ~ . + rcs(EFGR.09))
egfr21_10 <- update(base10, . ~ . + rcs(EFGR.21))

score_mcps10 <- riskRegression::Score(list("base" = base10,
                                           "egfr_09" = egfr09_10,
                                           "egfr_21" = egfr21_10), 
                                      formula = Surv(PERSON_YEARS, cv_mort) ~ 1, 
                                      data = mcps10, 
                                      conf.int = TRUE, 
                                      times = 9.99,
                                      cens.model = "km", 
                                      metrics = "brier",
                                      summary = "ipa")

score_10 <- score_mcps10[["Brier"]][["score"]]
score_10 <- score_10 %>% select(1,3,5,6,7)

####Discrimination and performance table####
tab_dis_perf <- tibble("Year" = c("Baseline", "2009", "2021"),
                       "Discrimination" = c(paste0(round(cstat_final[75,2], 3), " (95%CI ", 
                                                   round(cstat_final[75,7], 3), "-",
                                                   round(cstat_final[75,8], 3), ")"),
                                            paste0(round(cstat_final[58,2], 3), " (95%CI ", 
                                                   round(cstat_final[58,7], 3), "-",
                                                   round(cstat_final[58,8], 3), ")"),
                                            paste0(round(cstat_final[76,2], 3), " (95%CI ", 
                                                   round(cstat_final[76,7], 3), "-",
                                                   round(cstat_final[76,8], 3), ")")),
                       "Brier score (5 year)" = c(paste0(round(score_5[2,2], 3), " (95%CI ", 
                                                         round(score_5[2,3], 3), "-",
                                                         round(score_5[2,4], 3), ")"),
                                                  paste0(round(score_5[3,2], 3), " (95%CI ", 
                                                         round(score_5[3,3], 3), "-",
                                                         round(score_5[3,4], 3), ")"),
                                                  paste0(round(score_5[4,2], 3), " (95%CI ", 
                                                         round(score_5[4,3], 3), "-",
                                                         round(score_5[4,4], 3), ")")),
                       "Scaled Brier score (5 year)" = c(round(score_5[2,5], 4)*100, 
                                                         round(score_5[3,5], 4)*100, 
                                                         round(score_5[4,5], 4)*100),
                       "Brier score (10 year)" = c(paste0(round(score_10[2,2], 3), " (95%CI ", 
                                                          round(score_10[2,3], 3), "-",
                                                          round(score_10[2,4], 3), ")"),
                                                   paste0(round(score_10[3,2], 3), " (95%CI ", 
                                                          round(score_10[3,3], 3), "-",
                                                          round(score_10[3,4], 3), ")"),
                                                   paste0(round(score_10[4,2], 3), " (95%CI ", 
                                                          round(score_10[4,3], 3), "-",
                                                          round(score_10[4,4], 3), ")")),
                       "Scaled Brier score (10 year)" = c(round(score_10[2,5], 4)*100, 
                                                          round(score_10[3,5], 4)*100, 
                                                          round(score_10[4,5], 4)*100))

tab_dis_perf %>% 
  gt() %>% 
  gtsave(filename = "tabla_discriminacion_performance_cvd.docx")

table(mcps_fin$reclass_egfr) %>% prop.table()

#### Chi-squared and log-likelihood ####
summary(base5)$logtest
summary(egfr09_5)$logtest
summary(egfr21_5)$logtest

anova(base5,egfr09_5)
anova(base5,egfr21_5)
anova(egfr09_5,egfr21_5)

summary(base10)$logtest
summary(egfr09_10)$logtest
summary(egfr21_10)$logtest

anova(base10,egfr09_10)
anova(base10,egfr21_10)
anova(egfr09_10,egfr21_10)

#### --------- Kidney-related mortality mortality analyses ---------####
#####Administrative censoring####
fg_data <- finegray(Surv(PERSON_YEARS, kd_mort2) ~ ., data = mcps_fin)

mcps_pre <- survSplit(Surv(fgstart, fgstop, fgstatus) ~ ., data = fg_data, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

mcps_pre <- survSplit(Surv(fgstart, fgstop, fgstatus) ~ ., data = fg_data, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)

events10 <- mcps10 %>% filter(kd_mort == 1) %>% count(MALE, kd_mort)

#####Baseline Cox proportional hazards models####
#Five-year mortality risk
base5 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
               x = TRUE, y = TRUE, weight=fgwt,data = mcps5)

cph_base5 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x,
                 x = TRUE, y = TRUE, surv = TRUE, weight=fgwt,data = mcps5)

#Ten-year mortality risk
base10 <- coxph(Surv(fgstart, fgstop, fgstatus) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                x = TRUE, y = TRUE, weight=fgwt,data = mcps10)

cph_base10 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x,
                  x = TRUE, y = TRUE, surv = TRUE, weight=fgwt,data = mcps10)


#####5-year predictive mortality risk - Addition of eGFR estimated by 2009 CKD-EPI equation####
egfr09_5 <- update(base5, . ~ . + rcs(EFGR.09))
cox.zph(egfr09_5) %>% plot()

t1<-egfr09_5 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_05_09_kd.docx")
anova(base5, egfr09_5)

#Model fitted using rms::cph()
m5 <- datadist(mcps5)
options(datadist = "m5")

cph_egfr09_5 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.09), data = mcps5,
                    x = TRUE, y = TRUE, surv = TRUE,weight=fgwt)

specs(cph_egfr09_5) #specs() function gives knots location

#Graph
Predict(cph_egfr09_5) %>% 
  ggplot()

#####5-year predictive mortality risk - Addition of eGFR estimated by 2021 CKD-EPI equation####
egfr21_5 <- update(base5, . ~ . + rcs(EFGR.21))
cox.zph(egfr21_5) %>% plot()

t1<-egfr21_5 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_05_21_kd.docx")

anova(base5, egfr21_5)

#Model fitted using rms::cph()
cph_egfr21_5 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.21), data = mcps5,
                    x = TRUE, y = TRUE, surv = TRUE,weight=fgwt)

specs(cph_egfr21_5) #specs() function gives knots location

#Graph
Predict(cph_egfr21_5) %>% 
  ggplot()

#####5-year assessment of discrimination measures####
#Time range
c_base5 <- concordance(base5, keepstrata = TRUE)
c_egfr09_5 <- concordance(egfr09_5, keepstrata = TRUE)
c_egfr21_5 <- concordance(egfr21_5, keepstrata = TRUE)

#Fixed point for eGFR 2009
times <- seq(0.5, 5, by = 0.5)
data09 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base5, ymax = times[i])
  c_con <- concordance(egfr09_5, ymax = times[i])
  #Data
  data09[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_09_5 <- do.call(rbind, data09)

#Fixed point for eGFR 2021
times <- seq(0.5, 5, by = 0.5)
data21 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base5, ymax = times[i])
  c_con <- concordance(egfr21_5, ymax = times[i])
  #Data
  data21[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_21_5 <- do.call(rbind, data21)
cstat_5_final <- rbind(cstat_09_5, cstat_21_5)
cstat_5_final$equation <- c(rep("2009 equation", 20), rep("2021 equation", 20))
cstat_5_final$year <- c(rep("5-year", 40))
cstat_5_final$lower <- cstat_5_final$coef - 1.96 * cstat_5_final$se
cstat_5_final$upper <- cstat_5_final$coef + 1.96 * cstat_5_final$se

#####5-year Brier and scaled Brier score####
mcps_pre <- survSplit(Surv(PERSON_YEARS, kd_mort) ~ ., data = mcps_fin, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

base5 <- coxph(Surv(PERSON_YEARS, kd_mort) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
               x = TRUE, y = TRUE, data = mcps5)
egfr09_5 <- update(base5, . ~ . + rcs(EFGR.09))
egfr21_5 <- update(base5, . ~ . + rcs(EFGR.21))

score_mcps5 <- riskRegression::Score(list("base" = base5,
                                          "egfr_09" = egfr09_5,
                                          "egfr_21" = egfr21_5), 
                                     formula = Surv(PERSON_YEARS, cv_mort) ~ 1, 
                                     data = mcps5, 
                                     conf.int = TRUE, 
                                     times = 4.99,
                                     cens.model = "km", 
                                     metrics = "brier",
                                     summary = "ipa")

score_5 <- score_mcps5[["Brier"]][["score"]]
score_5 <- score_5 %>% select(1,3,5,6,7)

#####10-year predictive mortality risk - Addition of eGFR estimated by 2009 CKD-EPI equation####
egfr09_10 <- update(base10, . ~ . + rcs(EFGR.09))
cox.zph(egfr09_10) %>% plot()

t1<-egfr09_10 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_10_09_kd.docx")


anova(base10, egfr09_10)

#Model fitted using rms::cph()
m10 <- datadist(mcps10)
options(datadist = "m10")

cph_egfr09_10 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.09), data = mcps10,
                     x = TRUE, y = TRUE, surv = TRUE,weight=fgwt)

specs(cph_egfr09_10) #specs() function gives knots location

#Graph
Predict(cph_egfr09_10) %>% 
  ggplot()

#####10-year predictive mortality risk - Addition of eGFR estimated by 2021 CKD-EPI equation####
egfr21_10 <- update(base10, . ~ . + rcs(EFGR.21))
cox.zph(egfr21_10) %>% plot()

t1<-egfr21_10 %>%
  tbl_regression(exp=T) %>% 
  as_flex_table() 

# Guardar como .docx
read_docx() %>%
  body_add_flextable(t1) %>%
  print(target = "Tables/model_10_21_kd.docx")

anova(base10, egfr21_10)

#Model fitted using rms::cph()
cph_egfr21_10 <- cph(Surv(fgstart, fgstop, fgstatus) ~ strat(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x + rcs(EFGR.21), data = mcps10,
                     x = TRUE, y = TRUE, surv = TRUE,weight=fgwt)

specs(cph_egfr21_10) #specs() function gives knots location

#Graph
Predict(cph_egfr21_10) %>% 
  ggplot()

#####10-year assessment of discrimination measures####
#Time range
c_base10 <- concordance(base10, keepstrata = TRUE)
c_egfr09_10 <- concordance(egfr09_10, keepstrata = TRUE)
c_egfr21_10 <- concordance(egfr21_10, keepstrata = TRUE)

#Fixed point for eGFR 2009
times <- seq(6, 10, by = 0.5)
data09 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base10, ymax = times[i])
  c_con <- concordance(egfr09_10, ymax = times[i])
  #Data
  data09[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_09_10 <- do.call(rbind, data09)

#Fixed point for eGFR 2021
times <- seq(6, 10, by = 0.5)
data21 <- list()

for(i in 1:length(times)) {
  #Harrel C
  c_sin <- concordance(base10, ymax = times[i])
  c_con <- concordance(egfr21_10, ymax = times[i])
  #Data
  data21[[i]] <- data.frame("time" = c(rep(times[i], 2)),
                            "coef" = c(coef(c_sin), coef(c_con)),
                            "se" = c(sqrt(vcov(c_sin)), sqrt(vcov(c_con))),
                            "type" = c("Base model", "eGFR model"))
}

cstat_21_10 <- do.call(rbind, data21)
cstat_10_final <- rbind(cstat_09_10, cstat_21_10)
cstat_10_final$equation <- c(rep("2009 equation", 18), rep("2021 equation", 18))
cstat_10_final$year <- c(rep("10-year", 36))
cstat_10_final$lower <- cstat_10_final$coef - 1.96 * cstat_10_final$se
cstat_10_final$upper <- cstat_10_final$coef + 1.96 * cstat_10_final$se

cstat_final <- rbind(cstat_5_final, cstat_10_final)

#####10-year Brier and scaled Brier score####
mcps_pre <- survSplit(Surv(PERSON_YEARS, kd_mort) ~ ., data = mcps_fin, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)
base10 <- coxph(Surv(PERSON_YEARS, kd_mort) ~ strata(age_cat) + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                x = TRUE, y = TRUE, data = mcps10)
egfr09_10 <- update(base10, . ~ . + rcs(EFGR.09))
egfr21_10 <- update(base10, . ~ . + rcs(EFGR.21))

score_mcps10 <- riskRegression::Score(list("base" = base10,
                                           "egfr_09" = egfr09_10,
                                           "egfr_21" = egfr21_10), 
                                      formula = Surv(PERSON_YEARS, cv_mort) ~ 1, 
                                      data = mcps10, 
                                      conf.int = TRUE, 
                                      times = 9.99,
                                      cens.model = "km", 
                                      metrics = "brier",
                                      summary = "ipa")

score_10 <- score_mcps10[["Brier"]][["score"]]
score_10 <- score_10 %>% select(1,3,5,6,7)

####Discrimination and performance table####
tab_dis_perf <- tibble("Year" = c("Baseline", "2009", "2021"),
                       "Discrimination" = c(paste0(round(cstat_final[75,2], 3), " (95%CI ", 
                                                   round(cstat_final[75,7], 3), "-",
                                                   round(cstat_final[75,8], 3), ")"),
                                            paste0(round(cstat_final[58,2], 3), " (95%CI ", 
                                                   round(cstat_final[58,7], 3), "-",
                                                   round(cstat_final[58,8], 3), ")"),
                                            paste0(round(cstat_final[76,2], 3), " (95%CI ", 
                                                   round(cstat_final[76,7], 3), "-",
                                                   round(cstat_final[76,8], 3), ")")),
                       "Brier score (5 year)" = c(paste0(round(score_5[2,2], 3), " (95%CI ", 
                                                         round(score_5[2,3], 3), "-",
                                                         round(score_5[2,4], 3), ")"),
                                                  paste0(round(score_5[3,2], 3), " (95%CI ", 
                                                         round(score_5[3,3], 3), "-",
                                                         round(score_5[3,4], 3), ")"),
                                                  paste0(round(score_5[4,2], 3), " (95%CI ", 
                                                         round(score_5[4,3], 3), "-",
                                                         round(score_5[4,4], 3), ")")),
                       "Scaled Brier score (5 year)" = c(round(score_5[2,5], 4)*100, 
                                                         round(score_5[3,5], 4)*100, 
                                                         round(score_5[4,5], 4)*100),
                       "Brier score (10 year)" = c(paste0(round(score_10[2,2], 3), " (95%CI ", 
                                                          round(score_10[2,3], 3), "-",
                                                          round(score_10[2,4], 3), ")"),
                                                   paste0(round(score_10[3,2], 3), " (95%CI ", 
                                                          round(score_10[3,3], 3), "-",
                                                          round(score_10[3,4], 3), ")"),
                                                   paste0(round(score_10[4,2], 3), " (95%CI ", 
                                                          round(score_10[4,3], 3), "-",
                                                          round(score_10[4,4], 3), ")")),
                       "Scaled Brier score (10 year)" = c(round(score_10[2,5], 4)*100, 
                                                          round(score_10[3,5], 4)*100, 
                                                          round(score_10[4,5], 4)*100))

tab_dis_perf %>% 
  gt() %>% 
  gtsave(filename = "Tables/tabla_discriminacion_performance_kd.docx")


#### Chi-squared and log-likelihood ####
summary(base5)$logtest
summary(egfr09_5)$logtest
summary(egfr21_5)$logtest

anova(base5,egfr09_5)
anova(base5,egfr21_5)
anova(egfr09_5,egfr21_5)

summary(base10)$logtest
summary(egfr09_10)$logtest
summary(egfr21_10)$logtest

anova(base10,egfr09_10)
anova(base10,egfr21_10)
anova(egfr09_10,egfr21_10)

#### --------- Reclassification relevance  ####

##All cause mortality for eGFR categories
mcps_pre <- survSplit(Surv(PERSON_YEARS, D000) ~ ., data = mcps_fin, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

mcps_pre <- survSplit(Surv(PERSON_YEARS, D000) ~ ., data = mcps_fin, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)

reclassified_adj_5 <- coxph(Surv(PERSON_YEARS, D000) ~ reclass_egfr+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                            x = TRUE, y = TRUE, data = mcps5 %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

reclassified_adj_5_ckd <- coxph(Surv(PERSON_YEARS, D000) ~ reclass_ckd+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                            x = TRUE, y = TRUE, data = mcps5 %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

reclassified_adj_10 <- coxph(Surv(PERSON_YEARS, D000) ~ reclass_egfr+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                            x = TRUE, y = TRUE, data = mcps10 %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

reclassified_adj_10_ckd <- coxph(Surv(PERSON_YEARS, D000) ~ reclass_ckd+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                            x = TRUE, y = TRUE, data = mcps10 %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

### CV mortality for eGFR categories

fg_data <- finegray(Surv(PERSON_YEARS, cv_mort2) ~ ., data = mcps_fin)

mcps_pre <- survSplit(Surv(fgstart, fgstop, fgstatus) ~ ., data = fg_data, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

mcps_pre <- survSplit(Surv(fgstart, fgstop, fgstatus) ~ ., data = fg_data, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)

reclassified_adj_5_cv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ reclass_egfr+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                            x = TRUE, y = TRUE, weight=fgwt, data = mcps5 %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

reclassified_adj_5_ckd_cv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ reclass_ckd+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                                x = TRUE, y = TRUE, weight=fgwt, data = mcps5 %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

reclassified_adj_10_cv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ reclass_egfr+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                             x = TRUE, y = TRUE, weight=fgwt, data = mcps10 %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

reclassified_adj_10_ckd_cv <- coxph(Surv(fgstart, fgstop, fgstatus) ~ reclass_ckd+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                                 x = TRUE, y = TRUE, weight=fgwt, data = mcps10 %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

### Kidney-related mortality for eGFR categories

fg_data <- finegray(Surv(PERSON_YEARS, kd_mort2) ~ ., data = mcps_fin)

mcps_pre <- survSplit(Surv(fgstart, fgstop, fgstatus) ~ ., data = fg_data, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

mcps_pre <- survSplit(Surv(fgstart, fgstop, fgstatus) ~ ., data = fg_data, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)

reclassified_adj_5_kd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ reclass_egfr+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                               x = TRUE, y = TRUE,  weight=fgwt,data = mcps5 %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

reclassified_adj_5_ckd_kd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ reclass_ckd+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                                   x = TRUE, y = TRUE,  weight=fgwt,data = mcps5 %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

reclassified_adj_10_kd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ reclass_egfr+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                                x = TRUE, y = TRUE,  weight=fgwt,data = mcps10 %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

reclassified_adj_10_ckd_kd <- coxph(Surv(fgstart, fgstop, fgstatus) ~ reclass_ckd+strata(age_cat)+MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x, 
                                    x = TRUE, y = TRUE,  weight=fgwt,data = mcps10 %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

### All-cause mortality ###
models <- list(reclassified_adj_5,reclassified_adj_10)
fit_model0 <- function(model, group) {
  model %>% tbl_regression(exponentiate = TRUE, include = c("reclass_egfr"), label="reclass_egfr"~"Reclassification with CKD-EPI 2021") %>%
    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
results <- lapply(models,fit_model0)
row1 <- tbl_merge(results, tab_spanner = c("5-year", "10-year"))

### All-cause mortality ###
models <- list(reclassified_adj_5_ckd,reclassified_adj_10_ckd)
fit_model0 <- function(model, group) {
  model %>% tbl_regression(exponentiate = TRUE, include = c("reclass_ckd"), label="reclass_ckd"~"Reclassification with CKD-EPI 2021") %>%
    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
results <- lapply(models,fit_model0)
row1_ckd <- tbl_merge(results, tab_spanner = c("5-year", "10-year"))

### CV mortality ###
models <- list(reclassified_adj_5_cv,reclassified_adj_10_cv)
fit_model0 <- function(model, group) {
  model %>% tbl_regression(exponentiate = TRUE, include = c("reclass_egfr"), label="reclass_egfr"~"Reclassification with CKD-EPI 2021") %>%
    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
results <- lapply(models,fit_model0)
row2 <- tbl_merge(results, tab_spanner = c("5-year", "10-year"))

### CV mortality ###
models <- list(reclassified_adj_5_ckd_cv,reclassified_adj_10_ckd_cv)
fit_model0 <- function(model, group) {
  model %>% tbl_regression(exponentiate = TRUE, include = c("reclass_ckd"), label="reclass_ckd"~"Reclassification with CKD-EPI 2021") %>%
    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
results <- lapply(models,fit_model0)
row2_ckd <- tbl_merge(results, tab_spanner = c("5-year", "10-year"))

### KD mortality ###
models <- list(reclassified_adj_5_kd,reclassified_adj_10_kd)
fit_model0 <- function(model, group) {
  model %>% tbl_regression(exponentiate = TRUE, include = c("reclass_egfr"), label="reclass_egfr"~"Reclassification with CKD-EPI 2021") %>%
    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
results <- lapply(models,fit_model0)
row3 <- tbl_merge(results, tab_spanner = c("5-year", "10-year"))

### KD mortality ###
models <- list(reclassified_adj_5_ckd_kd,reclassified_adj_10_ckd_kd)
fit_model0 <- function(model, group) {
  model %>% tbl_regression(exponentiate = TRUE, include = c("reclass_ckd"), label="reclass_ckd"~"Reclassification with CKD-EPI 2021") %>%
    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
results <- lapply(models,fit_model0)
row3_ckd <- tbl_merge(results, tab_spanner = c("5-year", "10-year"))

#Save table
tab2 <-tbl_stack(list(row1, row1_ckd,row2,row2_ckd,row3,row3_ckd), group_header = c("All-cause death","All-cause death CKD","CVD death","CVD death CKD","KD death","KD death CKD"))%>% as_flex_table() %>% align(align = "center",part = "all") %>% autofit()
doc <- read_docx() %>% body_add_flextable(value = tab2, split = TRUE) %>%  body_end_section_landscape() %>%
  print(target = "Tables/tabla2_5.docx"); remove(row1, row2)


#### Event rates per reclassification category #####
mcps_pre <- survSplit(Surv(PERSON_YEARS, D000) ~ ., data = mcps_fin, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

mcps_pre <- survSplit(Surv(PERSON_YEARS, D000) ~ ., data = mcps_fin, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)

#Age and sex standardized rates
rates1.M <- mcps5 %>% select(PERSON_YEARS,D000,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Male") %>% group_by(D000) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,D000,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Female") %>% group_by(D000) %>% un.ag.st()
r1 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps5 %>% select(PERSON_YEARS,D000,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Male") %>% group_by(D000) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,D000,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Female") %>% group_by(D000) %>% un.ag.st()
r2 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps5 %>% select(PERSON_YEARS,D000,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Male") %>% group_by(D000) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,D000,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Female") %>% group_by(D000) %>% un.ag.st()
r3 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps5 %>% select(PERSON_YEARS,D000,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Male") %>% group_by(D000) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,D000,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Female") %>% group_by(D000) %>% un.ag.st()
r4 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

round(rbind(r1,r2,r3,r4),2)

#Age and sex standardized rates
rates1.M <- mcps10 %>% select(PERSON_YEARS,D000,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Male") %>% group_by(D000) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,D000,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Female") %>% group_by(D000) %>% un.ag.st()
r1 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps10 %>% select(PERSON_YEARS,D000,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Male") %>% group_by(D000) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,D000,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Female") %>% group_by(D000) %>% un.ag.st()
r2 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps10 %>% select(PERSON_YEARS,D000,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Male") %>% group_by(D000) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,D000,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Female") %>% group_by(D000) %>% un.ag.st()
r3 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps10 %>% select(PERSON_YEARS,D000,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Male") %>% group_by(D000) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,D000,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Female") %>% group_by(D000) %>% un.ag.st()
r4 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

round(rbind(r1,r2,r3,r4),2)

### CV mortality
mcps_pre <- survSplit(Surv(PERSON_YEARS, cv_mort) ~ ., data = mcps_fin, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

mcps_pre <- survSplit(Surv(PERSON_YEARS, cv_mort) ~ ., data = mcps_fin, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)

#Age and sex standardized rates
rates1.M <- mcps5 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Male") %>% group_by(cv_mort) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Female") %>% group_by(cv_mort) %>% un.ag.st()
r1 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps5 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Male") %>% group_by(cv_mort) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Female") %>% group_by(cv_mort) %>% un.ag.st()
r2 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps5 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Male") %>% group_by(cv_mort) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Female") %>% group_by(cv_mort) %>% un.ag.st()
r3 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps5 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Male") %>% group_by(cv_mort) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Female") %>% group_by(cv_mort) %>% un.ag.st()
r4 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

round(rbind(r1,r2,r3,r4),2)

#Age and sex standardized rates
rates1.M <- mcps10 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Male") %>% group_by(cv_mort) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Female") %>% group_by(cv_mort) %>% un.ag.st()
r1 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps10 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Male") %>% group_by(cv_mort) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Female") %>% group_by(cv_mort) %>% un.ag.st()
r2 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps10 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Male") %>% group_by(cv_mort) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Female") %>% group_by(cv_mort) %>% un.ag.st()
r3 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps10 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Male") %>% group_by(cv_mort) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,cv_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Female") %>% group_by(cv_mort) %>% un.ag.st()
r4 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

round(rbind(r1,r2,r3,r4),2)

### Kidney-related mortality

mcps_pre <- survSplit(Surv(PERSON_YEARS, kd_mort) ~ ., data = mcps_fin, cut = c(5), episode = "epoch")
mcps5 <- subset(mcps_pre, epoch == 1)

mcps_pre <- survSplit(Surv(PERSON_YEARS, kd_mort) ~ ., data = mcps_fin, cut = c(10), episode = "epoch")
mcps10 <- subset(mcps_pre, epoch == 1)

#Age and sex standardized rates
rates1.M <- mcps5 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Male") %>% group_by(kd_mort) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Female") %>% group_by(kd_mort) %>% un.ag.st()
r1 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps5 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Male") %>% group_by(kd_mort) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Female") %>% group_by(kd_mort) %>% un.ag.st()
r2 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps5 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Male") %>% group_by(kd_mort) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Female") %>% group_by(kd_mort) %>% un.ag.st()
r3 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps5 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Male") %>% group_by(kd_mort) %>% un.ag.st()
rates1.F <- mcps5 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Female") %>% group_by(kd_mort) %>% un.ag.st()
r4 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

round(rbind(r1,r2,r3,r4),2)

#Age and sex standardized rates
rates1.M <- mcps10 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Male") %>% group_by(kd_mort) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="No reclassification",MALE=="Female") %>% group_by(kd_mort) %>% un.ag.st()
r1 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps10 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Male") %>% group_by(kd_mort) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="eGFR reclassification",MALE=="Female") %>% group_by(kd_mort) %>% un.ag.st()
r2 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps10 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Male") %>% group_by(kd_mort) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD reclassification",MALE=="Female") %>% group_by(kd_mort) %>% un.ag.st()
r3 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

rates1.M <- mcps10 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2, reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Male") %>% group_by(kd_mort) %>% un.ag.st()
rates1.F <- mcps10 %>% select(PERSON_YEARS,kd_mort,MALE,age_cat2,reclass_egfr_ckd) %>% filter(reclass_egfr_ckd=="CKD (No reclassification)",MALE=="Female") %>% group_by(kd_mort) %>% un.ag.st()
r4 <- rbind(rates1.M, rates1.F) %>% un.ag.st.CI %>% round(2)

round(rbind(r1,r2,r3,r4),2)


#### --------- Resurvey data handling ####
mcps_r<-mcps %>% 
  filter(creat_resurvey > 0) %>% 
  select(PERSON_YEARS, D000, D003, D008, D015, cv_mort, EFGR.09, EFGR.21, CKD.09, CKD.21, CKD.09.plus, CKD.21.plus, EGFR_cat.09, EGFR_cat.21, reclass_egfr,reclass_ckd,reclass_egfr_ckd,
         AGE, age_cat, MALE,R_Age,R_MALE, COYOACAN,R_COYOACAN, EDU_LEVEL, R_SBP1,diab, smoke, cvd, SBP1, BASE_DIABETES, BASE_HYPERTENSION, BASE_CKD, R_RENTP,R_DIALYSIS,
         BASE_HBA1C, Total_C.x,Total_C.y, reclass_egfr, reclass_ckd, EFGR.09_R, EFGR.21_R, CKD.09_R, CKD.21_R, CKD.09.plus_R, CKD.21.plus_R, EGFR_cat.09_R, EGFR_cat.21_R,YEARS_BASELINE_RESURVEY) %>% 
  filter(!is.na(YEARS_BASELINE_RESURVEY) &!is.na(EFGR.09_R) & !is.na(EFGR.21_R)) %>%
  mutate(CKD.09.plus_R=ifelse(is.na(CKD.09.plus_R),"No CKD","CKD"),CKD.21.plus_R=ifelse(is.na(CKD.21.plus_R),"No CKD","CKD"))%>%
  mutate(ckd21=ifelse(CKD.21.plus_R=="CKD",1,0), ckd09=ifelse(CKD.09.plus_R=="CKD",1,0),reclass=ifelse(reclass_ckd=="CKD reclassificatio",1,0)) %>%
  mutate(progression_21=ifelse((CKD.21.plus=="No CKD" & CKD.21.plus_R =="CKD")==T,1,0),progression_09=ifelse((CKD.09.plus=="No CKD" & CKD.09.plus_R =="CKD")==T,1,0)) %>% 
  mutate(reclass_egfr_R = if_else(EGFR_cat.09_R != EGFR_cat.21_R, "Upward reclassification", "No reclassification"),
         reclass_ckd_R = if_else(CKD.09_R == CKD.21_R,"No reclassification","CKD reclassification"),
         reclass_egfr_ckd_R = case_when(reclass_egfr_R=="Upward reclassification" & reclass_ckd_R=="No reclassification"~"eGFR reclassification",
                                        reclass_egfr_R=="Upward reclassification" & reclass_ckd_R=="CKD reclassification"~"CKD reclassification",
                                        reclass_egfr_R=="No reclassification" & reclass_ckd_R=="No reclassification"~"No reclassification",
                                        reclass_egfr_R=="No reclassification" & reclass_ckd_R=="CKD reclassification"~"CKD reclassification")) %>%
  mutate(reclass_ckd_R=relevel(factor(reclass_ckd_R, ordered = F), ref="No reclassification"), reclass_egfr_ckd_R=relevel(factor(reclass_egfr_ckd_R, ordered = F), ref="No reclassification"),
         years2=PERSON_YEARS-YEARS_BASELINE_RESURVEY, age_cat_R=cut(R_Age, breaks = c(1, 35, 45, 55, 65, 75, Inf), right = FALSE))
levels(mcps_r$reclass_egfr_ckd_R)<-c(levels(mcps_r$reclass_egfr_ckd_R), "CKD (No reclassification)")
mcps_r$reclass_egfr_ckd_R[mcps_r$CKD.09_R == "CKD" & mcps_r$CKD.21_R =="CKD"]<-"CKD (No reclassification)"
table(mcps_r$reclass_egfr_ckd_R)

mcps_r$progression_egfr09<-0
mcps_r$progression_egfr09[mcps_r$EGFR_cat.09=="G1" & ((mcps_r$EGFR_cat.09_R %in% c("G2","G3a","G3b","G4","G4")))]<-1
mcps_r$progression_egfr09[mcps_r$EGFR_cat.09=="G2" & ((mcps_r$EGFR_cat.09_R %in% c("G3a","G3b","G4","G5")))]<-1
mcps_r$progression_egfr09[mcps_r$EGFR_cat.09=="G3a" & ((mcps_r$EGFR_cat.09_R %in% c("G3b","G4","G5")))]<-1
mcps_r$progression_egfr09[mcps_r$EGFR_cat.09=="G3b" & ((mcps_r$EGFR_cat.09_R %in% c("G4","G5")))]<-1
mcps_r$progression_egfr09[mcps_r$EGFR_cat.09=="G4" & ((mcps_r$EGFR_cat.09_R %in% c("G5")))]<-1
table(mcps_r$progression_egfr09)

mcps_r$progression_egfr21<-0
mcps_r$progression_egfr21[mcps_r$EGFR_cat.21=="G1" & ((mcps_r$EGFR_cat.21_R %in% c("G2","G3a","G3b","G4","G5")))]<-1
mcps_r$progression_egfr21[mcps_r$EGFR_cat.21=="G2" & ((mcps_r$EGFR_cat.21_R %in% c("G3a","G3b","G4","G5")))]<-1
mcps_r$progression_egfr21[mcps_r$EGFR_cat.21=="G3a" & ((mcps_r$EGFR_cat.21_R %in% c("G3b","G4","G5")))]<-1
mcps_r$progression_egfr21[mcps_r$EGFR_cat.21=="G3b" & ((mcps_r$EGFR_cat.21_R %in% c("G4","G5")))]<-1
mcps_r$progression_egfr21[mcps_r$EGFR_cat.21=="G4" & ((mcps_r$EGFR_cat.21_R %in% c("G5")))]<-1
table(mcps_r$progression_egfr21)

mcps_r$esrd09<-ifelse(mcps_r$EGFR_cat.09_R=="G5",1,0)
mcps_r$esrd21<-ifelse(mcps_r$EGFR_cat.21_R=="G5",1,0)

#### --------- Reclassification of CKD over time ####

m1<-coxph(Surv(YEARS_BASELINE_RESURVEY,progression_egfr09)~reclass_egfr+strata(age_cat) + MALE + R_COYOACAN + EDU_LEVEL + smoke + diab + cvd + R_SBP1 + Total_C.y, data=mcps_r %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))

m1_ckd<-coxph(Surv(YEARS_BASELINE_RESURVEY,progression_egfr09)~reclass_ckd+strata(age_cat) + MALE + R_COYOACAN + EDU_LEVEL + smoke + diab + cvd + R_SBP1 + Total_C.y, data=mcps_r %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")))
summary(m1_ckd)

m2<-coxph(Surv(YEARS_BASELINE_RESURVEY,progression_egfr21)~reclass_egfr_ckd+strata(age_cat) + MALE + R_COYOACAN + EDU_LEVEL + smoke + diab + cvd + R_SBP1 + Total_C.y, data=mcps_r)
summary(m2)

m1a<-coxph(Surv(YEARS_BASELINE_RESURVEY,ckd09)~reclass_egfr_ckd+strata(age_cat) + MALE + R_COYOACAN + EDU_LEVEL + smoke + diab + cvd + R_SBP1 + Total_C.y, data=mcps_r)
summary(m1)

m2b<-coxph(Surv(YEARS_BASELINE_RESURVEY,ckd21)~reclass_egfr_ckd+strata(age_cat) + MALE + R_COYOACAN + EDU_LEVEL + smoke + diab + cvd + R_SBP1 + Total_C.y, data=mcps_r)
summary(m2b)


### eGFR progression ###
models <- list(m1, m2)
fit_model0 <- function(model, group) {
  model %>% tbl_regression(exponentiate = TRUE, include = c("reclass_egfr_ckd"), label="reclass_egfr_ckd"~"Reclassification with CKD-EPI 2021") %>%
    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
results <- lapply(models,fit_model0)
row1 <- tbl_merge(results, tab_spanner = c("Progression using CKD-EPI 2009", "Progression using CKD-EPI 2021"))

### Progression to CKD ###
models <- list(m1a, m2b)
fit_model0 <- function(model, group) {
  model %>% tbl_regression(exponentiate = TRUE, include = c("reclass_egfr_ckd"), label="reclass_egfr_ckd"~"Reclassification with CKD-EPI 2021") %>%
    bold_labels() %>% italicize_levels() %>% modify_table_body(dplyr::select, -p.value)}
results <- lapply(models,fit_model0)
row2 <- tbl_merge(results, tab_spanner = c("Progression using CKD-EPI 2009", "Progression using CKD-EPI 2021"))


#Save table
tab2 <-tbl_stack(list(row1, row2), group_header = c("eGFR progression","CKD progression"))%>% as_flex_table() %>% align(align = "center",part = "all") %>% autofit()
doc <- read_docx() %>% body_add_flextable(value = tab2, split = TRUE) %>%  body_end_section_landscape() %>%
  print(target = "Tables/tabla2_0.docx"); remove(row1, row2)


#### --------- Adjusted Kaplan-Meier curves ####
mcps_fin$reclass_egfr<-factor(mcps_fin$reclass_egfr, labels = c("None", "Upward"))
mcps_fin$reclass_egfr_cat<-ifelse(mcps_fin$reclass_egfr=="Upward",1,0)
treatment_model <- glm(reclass_egfr_cat ~ AGE + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x,
                                  data = mcps_fin %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")), family="binomial")
s_iptw <- adjustedsurv(data=mcps_fin %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")),
                      variable = "reclass_egfr",
                      ev_time = "PERSON_YEARS",
                      event = "D000",
                      cause=1,
                      method = "iptw_km",
                      treatment_model=treatment_model,
                      conf_int = TRUE,
                      bootstrap = FALSE)

fig4a <- plot(s_iptw, risk_table = TRUE, risk_table_stratify_color = TRUE,
             risk_table_stratify = TRUE, risk_table_size = 4,risk_table_digits = 0, x_n_breaks = 8,
             risk_table_title_size= 10, ylab = "Adjusted all-cause mortality",cif = TRUE,
             gg_theme= theme_bw(), risk_table_theme = theme_pubclean(),
             legend.position="top",censoring_ind = "none",max_t = 10,
             xlab="Time in Years", additional_layers = list(scale_y_continuous(labels = scales::percent), labs(col="Reclassification")))


s_iptw2 <- adjustedsurv(data=mcps_fin %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")),
                       variable = "reclass_egfr",
                       ev_time = "PERSON_YEARS",
                       event = "cv_mort",
                       cause=1,
                       method = "iptw_km",
                       treatment_model=treatment_model,
                       conf_int = TRUE,
                       bootstrap = FALSE)

fig4b <- plot(s_iptw2, risk_table = TRUE, risk_table_stratify_color = TRUE,
              risk_table_stratify = TRUE, risk_table_size = 4,risk_table_digits = 0, x_n_breaks = 8,
              risk_table_title_size= 10, ylab = "Adjusted cardiovascular mortality",cif = TRUE,
              gg_theme= theme_bw(), risk_table_theme = theme_pubclean(),
              legend.position="top",censoring_ind = "none",max_t = 10,
              xlab="Time in Years", additional_layers = list(scale_y_continuous(labels = scales::percent), labs(col="Reclassification")))

s_iptw3 <- adjustedsurv(data=mcps_fin %>% filter(!reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")),
                        variable = "reclass_egfr",
                        ev_time = "PERSON_YEARS",
                        event = "kd_mort",
                        cause=1,
                        method = "iptw_km",
                        treatment_model=treatment_model,
                        conf_int = TRUE,
                        bootstrap = FALSE)

fig4c <- plot(s_iptw3, risk_table = TRUE, risk_table_stratify_color = TRUE,
              risk_table_stratify = TRUE, risk_table_size = 4,risk_table_digits = 0, x_n_breaks = 8,
              risk_table_title_size= 10, ylab = "Adjusted kidney-related mortality",cif = TRUE,
              gg_theme= theme_bw(), risk_table_theme = theme_pubclean(),
              legend.position="top",censoring_ind = "none",max_t = 10,
              xlab="Time in Years", additional_layers = list(scale_y_continuous(labels = scales::percent), labs(col="Reclassification")))
  
treatment_model2 <- glm(reclass_egfr_cat ~ AGE + MALE + COYOACAN + EDU_LEVEL + smoke + diab + cvd + SBP1 + Total_C.x,
                       data = mcps_fin %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")), family="binomial")
s_iptw <- adjustedsurv(data=mcps_fin %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")),
                       variable = "reclass_egfr",
                       ev_time = "PERSON_YEARS",
                       event = "D000",
                       cause=1,
                       method = "iptw_km",
                       treatment_model=treatment_model2,
                       conf_int = TRUE,
                       bootstrap = FALSE)

fig4d <- plot(s_iptw, risk_table = TRUE, risk_table_stratify_color = TRUE,
              risk_table_stratify = TRUE, risk_table_size = 4,risk_table_digits = 0, x_n_breaks = 8,
              risk_table_title_size= 10, ylab = "Adjusted all-cause mortality",cif = TRUE,
              gg_theme= theme_bw(), risk_table_theme = theme_pubclean(),
              legend.position="top",censoring_ind = "none",max_t = 10,
              xlab="Time in Years", additional_layers = list(scale_y_continuous(labels = scales::percent), labs(col="Reclassification")))


s_iptw2 <- adjustedsurv(data=mcps_fin %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")),
                        variable = "reclass_egfr",
                        ev_time = "PERSON_YEARS",
                        event = "cv_mort",
                        cause=1,
                        method = "iptw_km",
                        treatment_model=treatment_model,
                        conf_int = TRUE,
                        bootstrap = FALSE)

fig4e <- plot(s_iptw2, risk_table = TRUE, risk_table_stratify_color = TRUE,
              risk_table_stratify = TRUE, risk_table_size = 4,risk_table_digits = 0, x_n_breaks = 8,
              risk_table_title_size= 10, ylab = "Adjusted cardiovascular mortality",cif = TRUE,
              gg_theme= theme_bw(), risk_table_theme = theme_pubclean(),
              legend.position="top",censoring_ind = "none",max_t = 10,
              xlab="Time in Years", additional_layers = list(scale_y_continuous(labels = scales::percent), labs(col="Reclassification")))

s_iptw3 <- adjustedsurv(data=mcps_fin %>% filter(reclass_egfr_ckd %in% c("CKD (No reclassification)","CKD reclassification")),
                        variable = "reclass_egfr",
                        ev_time = "PERSON_YEARS",
                        event = "kd_mort",
                        cause=1,
                        method = "iptw_km",
                        treatment_model=treatment_model,
                        conf_int = TRUE,
                        bootstrap = FALSE)

fig4f <- plot(s_iptw3, risk_table = TRUE, risk_table_stratify_color = TRUE,
              risk_table_stratify = TRUE, risk_table_size = 4,risk_table_digits = 0, x_n_breaks = 8,
              risk_table_title_size= 10, ylab = "Adjusted kidney-related mortality",cif = TRUE,
              gg_theme= theme_bw(), risk_table_theme = theme_pubclean(),
              legend.position="top",censoring_ind = "none",max_t = 10,
              xlab="Time in Years", additional_layers = list(scale_y_continuous(labels = scales::percent), labs(col="Reclassification")))


fig4<-ggarrange(fig4a, fig4b,fig4c, fig4d,fig4e,fig4f, common.legend = T, labels = LETTERS[1:6])

ggsave(fig4,filename = "Figure4.pdf", 
       width = 60, 
       height = 30,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE) 

#### --------- Unadjusted Kaplan-Meier curves ####
km1 <-survfit(Surv(PERSON_YEARS, D000) ~ 1, data = mcps_fin)
f1 <- ggsurvplot(km1, fontsize = 4.5*0.75, data = mcps_fin, conf.int = T, risk.table = T, xlab="Time (years)", censor=F,
                 ylab="All-cause mortality (%)", title="Rates of all-cause mortality in MCPS",fun = "event",cumevents = T,
                 break.y.by= c(0.1), xlim=c(0,10),break.x.by= c(2), ggtheme =  ggpubr::theme_pubclean() + 
                   (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15*0.8),
                          text = element_text(hjust=0.5, family = "sans", size=15*0.75),
                          legend.text = element_text(hjust=0.5, size=11))))
f1$cumevents <- f1$cumevents + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5,  margin = margin(r=20*0.7), face="bold", size=12*0.75, angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size=15*0.75)) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
f1$plot <- f1$plot + labs(fill=NULL, color=NULL) + theme(legend.position="none") + scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,0.1))
f1$cumevents <- f1$cumevents  + xlab("") + theme(panel.grid = element_blank(), title = element_blank(), text = element_text(hjust=0.5, family = "sans", size=8)) +
  scale_x_continuous(labels = NULL,breaks = NULL,sec.axis = dup_axis(name = element_blank(),
                                                                     breaks = c(seq(0,20, by=2)), labels = c(seq(0,20, by=2))))
fig1a <- ggarrange(f1$plot, f1$cumevents, heights = c(2, 0.4), ncol = 1, nrow = 2)

km2 <-survfit(Surv(PERSON_YEARS, cv_mort) ~ 1, data = mcps_fin)
f2 <- ggsurvplot(km2, fontsize = 4.5*0.75, data = mcps_fin, conf.int = T, risk.table = T, xlab="Time (years)", censor=F,
                 ylab="Cardiovascular disease mortality (%)", title="Rates of cardiovascular disease mortality in MCPS",fun = "event",cumevents = T,
                 break.y.by= c(0.1), xlim=c(0,10),break.x.by= c(2), ggtheme =  ggpubr::theme_pubclean() + 
                   (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15*0.8),
                          text = element_text(hjust=0.5, family = "sans", size=15*0.75),
                          legend.text = element_text(hjust=0.5, size=11))))
f2$cumevents <- f2$cumevents + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5,  margin = margin(r=20*0.7), face="bold", size=12*0.75, angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size=15*0.75)) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
f2$plot <- f2$plot + labs(fill=NULL, color=NULL) + theme(legend.position="none") + scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,0.04))
f2$cumevents <- f2$cumevents  + xlab("") + theme(panel.grid = element_blank(), title = element_blank(), text = element_text(hjust=0.5, family = "sans", size=8)) +
  scale_x_continuous(labels = NULL,breaks = NULL,sec.axis = dup_axis(name = element_blank(),
                                                                     breaks = c(seq(0,20, by=2)), labels = c(seq(0,20, by=2))))
fig1b <- ggarrange(f2$plot, f2$cumevents, heights = c(2, 0.4), ncol = 1, nrow = 2)

km3 <-survfit(Surv(PERSON_YEARS, kd_mort) ~ 1, data = mcps_fin)
f3 <- ggsurvplot(km3, fontsize = 4.5*0.75, data = mcps_fin, conf.int = T, risk.table = T, xlab="Time (years)", censor=F,
                 ylab="Kidney disease mortality (%)", title="Rates of kidney disease mortality in MCPS",fun = "event",cumevents = T,
                 break.y.by= c(0.1), xlim=c(0,10),break.x.by= c(2), ggtheme =  ggpubr::theme_pubclean() + 
                   (theme(plot.title = element_text(hjust = 0.5, face = "bold", size=15*0.8),
                          text = element_text(hjust=0.5, family = "sans", size=15*0.75),
                          legend.text = element_text(hjust=0.5, size=11))))
f3$cumevents <- f3$cumevents + theme(axis.text.y = element_blank()) +
  theme(axis.text.y = element_text(hjust=0.5, vjust=0.5,  margin = margin(r=20*0.7), face="bold", size=12*0.75, angle=0)) +
  theme(plot.title = element_text(hjust=0.5,vjust=0,margin = margin(b = 0.5), face="bold", colour="black", size=15*0.75)) +
  theme(axis.ticks.y = element_blank()) + theme(axis.title = element_blank())
f3$plot <- f3$plot + labs(fill=NULL, color=NULL) + theme(legend.position="none") + scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.015))
f3$cumevents <- f3$cumevents  + xlab("") + theme(panel.grid = element_blank(), title = element_blank(), text = element_text(hjust=0.5, family = "sans", size=8)) +
  scale_x_continuous(labels = NULL,breaks = NULL,sec.axis = dup_axis(name = element_blank(),
                                                                     breaks = c(seq(0,20, by=2)), labels = c(seq(0,20, by=2))))
fig1c <- ggarrange(f3$plot, f3$cumevents, heights = c(2, 0.4), ncol = 1, nrow = 2)


fig1_a<-ggarrange(fig1a, fig1b, labels = LETTERS[1:2], nrow=1, ncol=2)
fig1_b<-ggarrange(NULL, fig1c, NULL,labels = c("","C",""), nrow=1, ncol=3, widths = c(0.25,0.5,0.25))
fig1<-ggarrange(fig1_a,"",fig1_b, nrow=3, ncol=1,heights = c(0.49,0.02,0.49))
ggsave(fig1,filename = "SuppFig1.jpg", 
       width = 25, 
       height = 20,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE) 
