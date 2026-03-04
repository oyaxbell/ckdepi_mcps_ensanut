## CHRONIC KIDNEY DISEASE PREVALENCE IN MEXICO
#Analysis: Paulina Sanchez-Castro & Carlos A. Fermin-Martinez
#Latest version of analysis: February 2023

#### Custom function####
wt_prop <- function(x, design) {
  f <- as.formula(paste0("~", x))
  prop <- svymean(f, design, na.rm = TRUE)
  ci <- confint(prop)
  table <- cbind(prop, ci)
  df <- data.frame(table)
}

wt_total <- function(x, design, year) {
  f <- as.formula(paste0("~", x))
  total <- svytotal(f, design, na.rm = TRUE)
  ci <- confint(total)
  table <- cbind(total, ci)
  df <- data.frame(table)
  df <- df %>% 
    rename(lower = X2.5..,
           upper = X97.5..) %>% 
    mutate(year = year)
  df
}

#### Load data----- ####
#setwd("C:/Users/pauli/OneDrive/Documentos")
#setwd("/Users/carlosfermin/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/Mexico City Cohort Study/Proyectos ENSANUT")
setwd("~/Mi unidad (obello@facmed.unam.mx)/Datasets/ENSANUT")

pacman::p_load(
  readr, haven, tidyverse, ggpubr, flextable,  gtsummary, forestmodel, nephro, survey, gt, officer)

pacman::p_load(conflicted)
conflicted::conflict_prefer(name = "filter", winner = "dplyr")
conflicted::conflict_prefer(name = "select", winner = "dplyr")

ensanut_2016 <- read_csv("ensanut2016_fin.csv")
ensanut_2018 <- read_csv("ensanut2018_fin.csv")
ensanut_2020 <- read_csv("ensanut2020_fin.csv")
ensanut_2021 <- read_csv("ensanut2021_fin.csv")
ensanut_2022 <- read_csv("ensanut2022_fin.csv")
ensanut_2023 <- read_csv("ensanut2023_fin.csv")

setwd("~/Mi unidad (obello@facmed.unam.mx)/Mexico City Cohort Study/Proyectos/2_Proyectos ENSANUT/CKD_EPI MX")

#### Data cleaning- ####
ensanut_2016$Sex <- ifelse(ensanut_2016$Sex==1, 0, 1)
ensanut_2018$Sex <- ifelse(ensanut_2018$Sex==1, 0, 1)
ensanut_2020$Sex <- ifelse(ensanut_2020$Sex==1, 0, 1)
ensanut_2021$Sex <- ifelse(ensanut_2021$Sex==1, 0, 1)
ensanut_2022$Sex <- ifelse(ensanut_2022$Sex==1, 0, 1)
ensanut_2023$Sex <- ifelse(ensanut_2023$Sex==1, 0, 1)

ensanut_2023$Year <- 2023

ensanut_2023$Creatinine1 <-  as.numeric(gsub(",", ".", ensanut_2023$Creatinine))

ensanut_2023$HBA1C_1 <- if_else(ensanut_2023$HBA1C < 10, ensanut_2023$HBA1C*10, ensanut_2023$HBA1C)
ensanut_2023$HBA1C_1 <- ensanut_2023$HBA1C_1 / 10

ensanut_2016$diabetes <- with(
  ensanut_2016, ifelse((Glucose>=126|HBA1C>=6.5|HX_T2D==1),1,0))
ensanut_2018$diabetes <- with(
  ensanut_2018, ifelse((Glucose>=126|HBA1C>=6.5|HX_T2D==1),1,0))
ensanut_2020$HX_T2D[is.na(ensanut_2020$HX_T2D)] <- 0
ensanut_2020$diabetes <- with(
  ensanut_2020, ifelse((Glucose>=126|HBA1C>=6.5|HX_T2D==1),1,0))
ensanut_2021$diabetes <- with(
  ensanut_2021, ifelse((Glucose>=126|HBA1C>=6.5|HX_T2D==1),1,0))
ensanut_2022$diabetes <- with(
  ensanut_2022, ifelse((Glucose>=126|HBA1C>=6.5|HX_T2D==1),1,0))
ensanut_2023$diabetes <- with(
  ensanut_2023, ifelse((Glucose>=126|HBA1C_1>=6.5|HX_T2D==1),1,0))

ensanut_2016$hypertension <- with(
  ensanut_2016, ifelse((SBP>=140|DBP>=90|HX_HBP==1),1,0))
ensanut_2018$hypertension <- with(
  ensanut_2018, ifelse((SBP>=140|DBP>=90|HX_HBP==1),1,0))
ensanut_2020$HX_HBP[is.na(ensanut_2020$HX_HBP)] <- 0
ensanut_2020$hypertension <- with(
  ensanut_2020, ifelse((SBP>=140|DBP>=90|HX_HBP==1),1,0))
ensanut_2021$hypertension <- with(
  ensanut_2021, ifelse((SBP>=140|DBP>=90|HX_HBP==1),1,0))
ensanut_2022$hypertension <- with(
  ensanut_2022, ifelse((SBP>=140|DBP>=90|HX_HBP==1),1,0))
ensanut_2023$hypertension <- with(
  ensanut_2023, ifelse((SBP>=140|DBP>=90|HX_HBP==1),1,0))

ensanut_2016$ethnicity <- rep(0, nrow(ensanut_2016))
ensanut_2018$ethnicity <- rep(0, nrow(ensanut_2018))
ensanut_2020$ethnicity <- rep(0, nrow(ensanut_2020))
ensanut_2021$ethnicity <- rep(0, nrow(ensanut_2021))
ensanut_2022$ethnicity <- rep(0, nrow(ensanut_2022))
ensanut_2023$ethnicity <- rep(0, nrow(ensanut_2023))

##-- Estimated glomerular filtration rate (eGFR) --##
#CKD-EPI 2009
ensanut_2016$egfr_12 <- with(ensanut_2016, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut_2018$egfr_12 <- with(ensanut_2018, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut_2020$egfr_12 <- with(ensanut_2020, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut_2021$egfr_12 <- with(ensanut_2021, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut_2022$egfr_12 <- with(ensanut_2022, CKDEpi.creat(
  creatinine=Creatinine, sex=Sex, age=Age, ethnicity=ethnicity))
ensanut_2023$egfr_12 <- with(ensanut_2023, CKDEpi.creat(
  creatinine=Creatinine1, sex=Sex, age=Age, ethnicity=ethnicity))

#CKD-EPI 2021
ensanut_2016$egfr_21 <- with(ensanut_2016, CKDEpi2021.creat(
  creatinine=Creatinine, sex=Sex, age=Age))
ensanut_2018$egfr_21 <- with(ensanut_2018, CKDEpi2021.creat(
  creatinine=Creatinine, sex=Sex, age=Age))
ensanut_2020$egfr_21 <- with(ensanut_2020, CKDEpi2021.creat(
  creatinine=Creatinine, sex=Sex, age=Age))
ensanut_2021$egfr_21 <- with(ensanut_2021, CKDEpi2021.creat(
  creatinine=Creatinine, sex=Sex, age=Age))
ensanut_2022$egfr_21 <- with(ensanut_2022, CKDEpi2021.creat(
  creatinine=Creatinine, sex=Sex, age=Age))
ensanut_2023$egfr_21 <- with(ensanut_2023, CKDEpi2021.creat(
  creatinine=Creatinine1, sex=Sex, age=Age))

ensanut_2022$egfr_12[is.infinite(ensanut_2022$egfr_12)] <- NA
ensanut_2022$egfr_21[is.infinite(ensanut_2022$egfr_21)] <- NA

##-- eGFR categories --##
#CKD-EPI 2009
ensanut_2016$egfr_12_cat <- cut(ensanut_2016$egfr_12, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2018$egfr_12_cat <- cut(ensanut_2018$egfr_12, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2020$egfr_12_cat <- cut(ensanut_2020$egfr_12, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2021$egfr_12_cat <- cut(ensanut_2021$egfr_12, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2022$egfr_12_cat <- cut(ensanut_2022$egfr_12, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2023$egfr_12_cat <- cut(ensanut_2023$egfr_12, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))

#CKD-EPI 2021
ensanut_2016$egfr_21_cat <- cut(ensanut_2016$egfr_21, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2018$egfr_21_cat <- cut(ensanut_2018$egfr_21, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2020$egfr_21_cat <- cut(ensanut_2020$egfr_21, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2021$egfr_21_cat <- cut(ensanut_2021$egfr_21, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2022$egfr_21_cat <- cut(ensanut_2022$egfr_21, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))
ensanut_2023$egfr_21_cat <- cut(ensanut_2023$egfr_21, breaks = c(-Inf,15,30,45,60,90,Inf),
                                right = FALSE, labels = c("G5","G4","G3b","G3a","G2","G1"))

##-- Chronic Kidney Disease (CKD) assessed by eGFR --##
#CKD-EPI 2009
ensanut_2016$CKD_12<-ifelse(ensanut_2016$egfr_12<60,1,0)
ensanut_2018$CKD_12<-ifelse(ensanut_2018$egfr_12<60,1,0)
ensanut_2020$CKD_12<-ifelse(ensanut_2020$egfr_12<60,1,0)
ensanut_2021$CKD_12<-ifelse(ensanut_2021$egfr_12<60,1,0)
ensanut_2022$CKD_12<-ifelse(ensanut_2022$egfr_12<60,1,0)
ensanut_2023$CKD_12<-ifelse(ensanut_2023$egfr_12<60,1,0)

#CKD-EPI 2012
ensanut_2016$CKD_21<-ifelse(ensanut_2016$egfr_21<60,1,0)
ensanut_2018$CKD_21<-ifelse(ensanut_2018$egfr_21<60,1,0)
ensanut_2020$CKD_21<-ifelse(ensanut_2020$egfr_21<60,1,0)
ensanut_2021$CKD_21<-ifelse(ensanut_2021$egfr_21<60,1,0)
ensanut_2022$CKD_21<-ifelse(ensanut_2022$egfr_21<60,1,0)
ensanut_2023$CKD_21<-ifelse(ensanut_2023$egfr_21<60,1,0)

##-- CKD: assessed by eGFR + baseline self-report --##
#CKD-EPI 2009
ensanut_2016$CKD_12.B<-(ensanut_2016 %>% select(HX_CKD, CKD_12) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)
ensanut_2018$CKD_12.B<-(ensanut_2018 %>% select(HX_CKD, CKD_12) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)
ensanut_2020$CKD_12.B <- ensanut_2020$CKD_12
ensanut_2021$CKD_12.B<-(ensanut_2021 %>% select(HX_CKD, CKD_12) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)
ensanut_2022$CKD_12.B<-(ensanut_2022 %>% select(HX_CKD, CKD_12) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)
ensanut_2023$CKD_12.B<-(ensanut_2023 %>% select(HX_CKD, CKD_12) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)
#CKD-EPI 2021
ensanut_2016$CKD_21.B<-(ensanut_2016 %>% select(HX_CKD, CKD_21) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)
ensanut_2018$CKD_21.B<-(ensanut_2018 %>% select(HX_CKD, CKD_21) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)
ensanut_2020$CKD_21.B <- ensanut_2020$CKD_21
ensanut_2021$CKD_21.B<-(ensanut_2021 %>% select(HX_CKD, CKD_21) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)
ensanut_2022$CKD_21.B<-(ensanut_2022 %>% select(HX_CKD, CKD_21) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)
ensanut_2023$CKD_21.B<-(ensanut_2023 %>% select(HX_CKD, CKD_21) %>%
                          apply(1, sum, na.rm=F)==0) %>% ifelse(0,1)

#Recode datasets
ensanut_2016 <- ensanut_2016 %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = factor(Sex, 0:1, c("Women","Men")),
  "Age_cat" = cut(Age, c(-Inf,60,Inf)) %>% ordered(
    labels=c("<60 years","≥60 years")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf)) %>% ordered(
    labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")), 
  "DM2" = ordered(diabetes, 0:1, c("No diabetes", "Diabetes")), 
  "HBP" = ordered(hypertension,0:1, c("No hypertension", "Hypertension")))

ensanut_2018 <- ensanut_2018 %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = factor(Sex, 0:1, c("Women","Men")),
  "Age_cat" = cut(Age, c(-Inf,60,Inf)) %>% ordered(
    labels=c("<60 years","≥60 years")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf)) %>% ordered(
    labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")), 
  "DM2" = ordered(diabetes, 0:1, c("No diabetes", "Diabetes")), 
  "HBP" = ordered(hypertension,0:1, c("No hypertension", "Hypertension")))

ensanut_2020 <- ensanut_2020 %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = factor(Sex, 0:1, c("Women","Men")),
  "Age_cat" = cut(Age, c(-Inf,60,Inf)) %>% ordered(
    labels=c("<60 years","≥60 years")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf)) %>% ordered(
    labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")), 
  "DM2" = ordered(diabetes, 0:1, c("No diabetes", "Diabetes")), 
  "HBP" = ordered(hypertension,0:1, c("No hypertension", "Hypertension")))

ensanut_2021 <- ensanut_2021 %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = factor(Sex, 0:1, c("Women","Men")),
  "Age_cat" = cut(Age, c(-Inf,60,Inf)) %>% ordered(
    labels=c("<60 years","≥60 years")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf)) %>% ordered(
    labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")), 
  "DM2" = ordered(diabetes, 0:1, c("No diabetes", "Diabetes")), 
  "HBP" = ordered(hypertension,0:1, c("No hypertension", "Hypertension")))

ensanut_2022 <- ensanut_2022 %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = factor(Sex, 0:1, c("Women","Men")),
  "Age_cat" = cut(Age, c(-Inf,60,Inf)) %>% ordered(
    labels=c("<60 years","≥60 years")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf)) %>% ordered(
    labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")), 
  "DM2" = ordered(diabetes, 0:1, c("No diabetes", "Diabetes")), 
  "HBP" = ordered(hypertension,0:1, c("No hypertension", "Hypertension")))

ensanut_2023 <- ensanut_2023 %>% mutate(
  "Smoking" = ordered(Smoking, 0:2, c("Never", "Former", "Current")),
  "Sex2" = factor(Sex, 0:1, c("Women","Men")),
  "Age_cat" = cut(Age, c(-Inf,60,Inf)) %>% ordered(
    labels=c("<60 years","≥60 years")),
  "BMI_cat" = cut(BMI, c(-Inf,30,Inf)) %>% ordered(
    labels=c("No obesity","Obesity")),
  "Indigenous2" = ordered(Indigenous, 0:1, c("Non-indigenous","Indigenous")), 
  "DM2" = ordered(diabetes, 0:1, c("No diabetes", "Diabetes")), 
  "HBP" = ordered(hypertension,0:1, c("No hypertension", "Hypertension")))

#### Reclassification in 2023####
ensanut_2023 <- ensanut_2023 %>% 
  mutate(reclass_egfr = if_else(egfr_12_cat != egfr_21_cat, "Upward reclassification", "No reclassification"),
         reclass_ckd = if_else(CKD_12 != CKD_21, "CKD reclassification", "No reclassification"))

#### Survey health- ####
#Adjust options for single PSU
options(survey.adjust.domain.lonely = TRUE, survey.lonely.psu="adjust")

#---------------------------------- 2016 ----------------------------------#
#Filter 1: age ≥20 years
pob_2016<-nrow(ensanut_2016)
ensanut_2016 <- ensanut_2016 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2016_2 <- ensanut_2016 %>% filter(!is.na(ponde_f.x.x), !is.na(est_var))
#Sampling design
ensanut_2016_survey_h <- svydesign(data=ensanut_2016_2, ids=~code_upm.x, #Data & PSU
                                  weights=~ponde_f.x.x, #Venous blood weights
                                  strata=~est_var, nest=TRUE) #Strata 1st stage
#CKD prevalence
svymean(~HX_CKD, ensanut_2016_survey_h, na.rm = T)*100

#Flowchart
pob_2016
nrow(ensanut_2016)
nrow(ensanut_2016 %>% filter(!is.na(Creatinine)))

#---------------------------------- 2018 ----------------------------------#
#Filter 1: age ≥20 years
pob_2018<-nrow(ensanut_2018)
ensanut_2018 <- ensanut_2018 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2018_2 <- ensanut_2018 %>% filter(!is.na(pond_sal), !is.na(est_sal))
#Sampling design
ensanut_2018_survey_h <- svydesign(data=ensanut_2018_2, ids=~UPM_sal, #Data & PSU
                                  weights=~pond_sal, #Venous blood weights
                                  strata=~est_sal, nest=TRUE) #Strata 1st stage
#CKD prevalence
svymean(~HX_CKD, ensanut_2018_survey_h, na.rm = T)*100

#Flowchart
pob_2018
nrow(ensanut_2018)
nrow(ensanut_2018 %>% filter(!is.na(Creatinine)))

#---------------------------------- 2021 ----------------------------------#
#Filter 1: age ≥20 years
pob_2021<-nrow(ensanut_2021)
ensanut_2021 <- ensanut_2021 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2021_2 <- ensanut_2021 %>% filter(!is.na(pond_sal), !is.na(est_final))
#Sampling design
ensanut_2021_survey_h <- svydesign(data=ensanut_2021_2, ids=~upm.x, #Data & PSU
                                  weights=~pond_sal, #Venous blood weights
                                  strata=~est_final, nest=TRUE) #Strata 1st stage
#CKD prevalence
svymean(~HX_CKD, ensanut_2021_survey_h, na.rm = T)*100

#---------------------------------- 2022 ----------------------------------#
#Filter 1: age ≥20 years
pob_2022<-nrow(ensanut_2022)
ensanut_2022 <- ensanut_2022 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2022_2 <- ensanut_2022 %>% filter(!is.na(ponde_f.y), !is.na(est_sel.y))
#Sampling design
ensanut_2022_survey_h <- svydesign(data=ensanut_2022_2, ids=~upm.y, #Data & PSU
                                  weights=~ponde_f.y, #Venous blood weights
                                  strata=~est_sel.y, nest=TRUE) #Strata 1st stage
#CKD prevalence
svymean(~HX_CKD, ensanut_2022_survey_h, na.rm = T)*100

#Flowchart
pob_2022
nrow(ensanut_2022)
nrow(ensanut_2022 %>% filter(!is.na(Creatinine)))
#---------------------------------- 2023 ----------------------------------#
#Filter 1: age ≥20 years
ensanut_2023 <- ensanut_2023 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2023_2 <- ensanut_2023 %>% filter(!is.na(ponde_f.y), !is.na(est_sel.y))
#Sampling design
ensanut_2023_survey_h <- svydesign(data=ensanut_2023_2, ids=~upm.y, #Data & PSU
                                   weights=~ponde_f.y, #Venous blood weights
                                   strata=~est_sel.y, nest=TRUE) #Strata 1st stage

#CKD prevalence
svymean(~HX_CKD, ensanut_2023_survey_h, na.rm = T)*100

#### Survey venous- ####
#Adjust options for single PSU
options(survey.adjust.domain.lonely = TRUE, survey.lonely.psu="adjust")

#---------------------------------- 2016 ----------------------------------#
#Filter 1: age ≥20 years
ensanut_2016 <- ensanut_2016 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2016_2 <- ensanut_2016 %>% filter(!is.na(ponde_f_vv), !is.na(est_var))
#Sampling design
ensanut_2016_survey2 <- svydesign(data=ensanut_2016_2, ids=~code_upm.x, #Data & PSU
                               weights=~ponde_f_vv, #Venous blood weights
                               strata=~est_var, nest=TRUE) #Strata 1st stage
#Filter 3: complete eGFR data
ensanut_2016_survey <- subset(ensanut_2016_survey2,!is.na(egfr_12))
#CKD prevalence
svymean(~CKD_12+CKD_21, ensanut_2016_survey, na.rm = T)*100

#---------------------------------- 2018 ----------------------------------#
#Filter 1: age ≥20 years
ensanut_2018 <- ensanut_2018 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2018_2 <- ensanut_2018 %>% filter(!is.na(pond_san), !is.na(est_san))
#Sampling design
ensanut_2018_survey2 <- svydesign(data=ensanut_2018_2, ids=~UPM_san, #Data & PSU
                               weights=~pond_san, #Venous blood weights
                               strata=~est_san, nest=TRUE) #Strata 1st stage
#Filter 3: complete eGFR data
ensanut_2018_survey <- subset(ensanut_2018_survey2, !is.na(egfr_12))
#CKD prevalence
svymean(~CKD_12+CKD_21, ensanut_2018_survey, na.rm = T)*100

#---------------------------------- 2020 ----------------------------------#
#Filter 1: age ≥20 years
ensanut_2020 <- ensanut_2020 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2020_2 <- ensanut_2020 %>% filter(!is.na(pond_lab), !is.na(est_final))
#Sampling design
ensanut_2020_survey2 <- svydesign(data=ensanut_2020_2, ids=~Upm.x, #Data & PSU
                               weights=~pond_lab, #Venous blood weights
                               strata=~est_final, nest=TRUE) #Strata 1st stage
#Filter 3: complete eGFR data
ensanut_2020_survey <- subset(ensanut_2020_survey2, !is.na(egfr_12))
#CKD prevalence
svymean(~CKD_12+CKD_21, ensanut_2020_survey, na.rm = T)*100

#---------------------------------- 2021 ----------------------------------#
#Filter 1: age ≥20 years
ensanut_2021 <- ensanut_2021 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2021_2 <- ensanut_2021 %>% filter(!is.na(pond_lab), !is.na(est_final))
#Sampling design
ensanut_2021_survey2 <- svydesign(data=ensanut_2021_2, ids=~upm.x, #Data & PSU
                               weights=~pond_lab, #Venous blood weights
                               strata=~est_final, nest=TRUE) #Strata 1st stage
#Filter 3: complete eGFR data
ensanut_2021_survey <- subset(ensanut_2021_survey2, !is.na(egfr_12))
#CKD prevalence
svymean(~CKD_12+CKD_21, ensanut_2021_survey, na.rm = T)*100

#---------------------------------- 2022 ----------------------------------#
#Filter 1: age ≥20 years
ensanut_2022 <- ensanut_2022 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2022_2 <- ensanut_2022 %>% filter(!is.na(ponde_v), !is.na(est_sel.y))
#Sampling design
ensanut_2022_survey2 <- svydesign(data=ensanut_2022_2, ids=~upm.y, #Data & PSU
                               weights=~ponde_v, #Venous blood weights
                               strata=~est_sel.y, nest=TRUE) #Strata 1st stage
#Filter 3: complete eGFR data
ensanut_2022_survey <- subset(ensanut_2022_survey2, !is.na(egfr_12))
#CKD prevalence
svymean(~CKD_12+CKD_21, ensanut_2022_survey, na.rm = T)*100

#---------------------------------- 2023 ----------------------------------#
#Filter 1: age ≥20 years
ensanut_2023 <- ensanut_2023 %>% filter(Age>=20)
#Filter 2: complete data for weights and strata
ensanut_2023_2 <- ensanut_2023 %>% filter(!is.na(ponde_suero_2), !is.na(est_sel.y.y))
#Sampling design
ensanut_2023_survey2 <- svydesign(data=ensanut_2023_2, ids=~upm.y.y, #Data & PSU
                                  weights=~ponde_suero_2, #Venous blood weights
                                  strata=~est_sel.y.y, nest=TRUE) #Strata 1st stage
#Filter 3: complete eGFR data
ensanut_2023_survey <- subset(ensanut_2023_survey2, !is.na(egfr_12))
#CKD prevalence
svymean(~CKD_12+CKD_21, ensanut_2023_survey, na.rm = T)*100


#Additional prevalence
#2016
svymean(~HX_CKD, ensanut_2016_survey_h, na.rm = T)*100; svymean(
  ~CKD_12.B+CKD_21.B, ensanut_2016_survey, na.rm = T)*100
#2018
svymean(~HX_CKD, ensanut_2018_survey_h, na.rm = T)*100; svymean(
  ~CKD_12.B+CKD_21.B, ensanut_2018_survey, na.rm = T)*100
#2021
svymean(~HX_CKD, ensanut_2021_survey_h, na.rm = T)*100; svymean(
  ~CKD_12.B+CKD_21.B, ensanut_2021_survey, na.rm = T)*100
#2022
svymean(~HX_CKD, ensanut_2022_survey_h, na.rm = T)*100; svymean(
  ~CKD_12.B+CKD_21.B, ensanut_2022_survey, na.rm = T)*100

#### Flowchart ####

#2016
pob_2021
nrow(ensanut_2021)
nrow(ensanut_2021 %>% filter(!is.na(Creatinine)))

#2018
pob_2021
nrow(ensanut_2021)
nrow(ensanut_2021 %>% filter(!is.na(Creatinine)))

#2021
pob_2021
nrow(ensanut_2021)
nrow(ensanut_2021 %>% filter(!is.na(Creatinine)))

#2022
pob_2021
nrow(ensanut_2021)
nrow(ensanut_2021 %>% filter(!is.na(Creatinine)))

#2023
pob_2021
nrow(ensanut_2021)
nrow(ensanut_2021 %>% filter(!is.na(Creatinine)))

####---------------------------- eGFR CATEGORIES PREVALENCE: CKD-EPI 2021 ------------------#### ----####
#### Prevalence by year ---- ####
####ENSANUT 2016
egfr_21_prev_2016 <- wt_prop("egfr_21_cat", ensanut_2016_survey)

####ENSANUT 2018
egfr_21_prev_2018 <- wt_prop("egfr_21_cat", ensanut_2018_survey)

####ENSANUT 2020
egfr_21_prev_2020 <- wt_prop("egfr_21_cat", ensanut_2020_survey)

####ENSANUT 2021
egfr_21_prev_2021 <- wt_prop("egfr_21_cat", ensanut_2021_survey)

####ENSANUT 2022
egfr_21_prev_2022 <- wt_prop("egfr_21_cat", ensanut_2022_survey)

####ENSANUT 2023
egfr_21_prev_2023 <- wt_prop("egfr_21_cat", ensanut_2023_survey)

#Total prevalence
egfr_21_prev_total <- rbind(egfr_21_prev_2016, egfr_21_prev_2018, egfr_21_prev_2020, egfr_21_prev_2021, egfr_21_prev_2022, egfr_21_prev_2023)
egfr_21_prev_total <- egfr_21_prev_total %>% 
  rename(lower = X2.5..,
         upper = X97.5..) %>% 
  mutate(cat = rep(c("G5","G4","G3b","G3a","G2","G1"), times = 6),
         cat = factor(cat, levels = c("G1", "G2", "G3a", "G3b", "G4", "G5")),
         year = rep(c(2016, 2018, 2020, 2021, 2022, 2023), each = 6),
         eq = rep(c("2021 equation"), times = 36),
         lower = if_else(lower < 0, 0, lower)) %>% 
  group_by(year) %>% 
  arrange(cat, .by_group = TRUE)

#### eGFR prevalence table -####
#eGFR Table
egfr_tab_21 <- egfr_21_prev_total %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         label = paste0(sprintf(prop*100, fmt = "%#.2f"), "%" , " (", sprintf(lower*100, fmt = "%#.2f"), "-", sprintf(upper*100, fmt = "%#.2f"), ")")) %>% 
  select(cat, year, label) %>% 
  pivot_wider(names_from = year, values_from = label) %>% 
  gt() %>% 
  gtsave(filename = "tabla_egfr_2021.docx")

####---------------------------- eGFR CATEGORIES PREVALENCE: CKD-EPI 2009 ------------------#### ----####
#### Prevalence by year ---- ####
####ENSANUT 2016
egfr_09_prev_2016 <- wt_prop("egfr_12_cat", ensanut_2016_survey)

####ENSANUT 2018
egfr_09_prev_2018 <- wt_prop("egfr_12_cat", ensanut_2018_survey)

####ENSANUT 2020
egfr_09_prev_2020 <- wt_prop("egfr_12_cat", ensanut_2020_survey)

####ENSANUT 2021
egfr_09_prev_2021 <- wt_prop("egfr_12_cat", ensanut_2021_survey)

####ENSANUT 2022
egfr_09_prev_2022 <- wt_prop("egfr_12_cat", ensanut_2022_survey)

####ENSANUT 2023
egfr_09_prev_2023 <- wt_prop("egfr_12_cat", ensanut_2023_survey)

#Total prevalence
egfr_09_prev_total <- rbind(egfr_09_prev_2016, egfr_09_prev_2018, egfr_09_prev_2020, egfr_09_prev_2021, egfr_09_prev_2022, egfr_09_prev_2023)
egfr_09_prev_total <- egfr_09_prev_total %>% 
  rename(lower = X2.5..,
         upper = X97.5..) %>% 
  mutate(cat = rep(c("G5","G4","G3b","G3a","G2","G1"), times = 6),
         cat = factor(cat, levels = c("G1", "G2", "G3a", "G3b", "G4", "G5")),
         year = rep(c(2016, 2018, 2020, 2021, 2022, 2023), each = 6),
         eq = rep(c("2009 equation"), times = 36),
         lower = if_else(lower < 0, 0, lower)) %>% 
  group_by(year) %>% 
  arrange(cat, .by_group = TRUE)

#### eGFR prevalence table -####
#eGFR Table
egfr_tab_09 <- egfr_09_prev_total %>% 
  mutate(lower = if_else(lower < 0, 0, lower),
         label = paste0(sprintf(prop*100, fmt = "%#.2f"), "%" , " (", sprintf(lower*100, fmt = "%#.2f"), "-", sprintf(upper*100, fmt = "%#.2f"), ")")) %>% 
  select(cat, year, label) %>% 
  pivot_wider(names_from = year, values_from = label) %>% 
  gt() %>% 
  gtsave(filename = "tabla_egfr_2009.docx")

####---------------------------- CKD PREVALENCE: CKD-EPI 2021 ------------------#### ----####
#### Prevalence by year ----####
ckd_21_prev_2016 <- wt_prop("CKD_21", ensanut_2016_survey)

####ENSANUT 2018
ckd_21_prev_2018 <- wt_prop("CKD_21", ensanut_2018_survey)

####ENSANUT 2020
ckd_21_prev_2020 <- wt_prop("CKD_21", ensanut_2020_survey)

####ENSANUT 2021
ckd_21_prev_2021 <- wt_prop("CKD_21", ensanut_2021_survey)

####ENSANUT 2022
ckd_21_prev_2022 <- wt_prop("CKD_21", ensanut_2022_survey)

####ENSANUT 2023
ckd_21_prev_2023 <- wt_prop("CKD_21", ensanut_2023_survey)

ckd_21_prev_total <- rbind(ckd_21_prev_2016, ckd_21_prev_2018, ckd_21_prev_2020, ckd_21_prev_2021, ckd_21_prev_2022, ckd_21_prev_2023)
ckd_21_prev_total <- ckd_21_prev_total %>% 
  rename(lower = X2.5..,
         upper = X97.5..) %>% 
  mutate(year = c(2016, 2018, 2020, 2021, 2022, 2023),
         lower = if_else(lower < 0, 0, lower))

ckf_21_tab <- ckd_21_prev_total %>% 
  mutate(label = paste0(sprintf(prop*100, fmt = "%#.2f"), "%" , " (", sprintf(lower*100, fmt = "%#.2f"), "-", sprintf(upper*100, fmt = "%#.2f"), ")")) %>% 
  select(year, label) %>% 
  pivot_wider(names_from = year, values_from = label)

####---------------------------- CKD PREVALENCE: CKD-EPI 2009 ------------------#### ----####
#### Prevalence by year-- ####
###ENSANUT 2016
prev_year_2016<-svyby(~CKD_12, by=~Year, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_year_2018<-svyby(~CKD_12, by=~Year, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_year_2020<-svyby(~CKD_12, by=~Year, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_year_2021<-svyby(~CKD_12, by=~Year, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

#ENSANUT 2022
prev_year_2022<-svyby(~CKD_12, by=~Year, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

#ENSANUT 2023
prev_year_2023<-svyby(~CKD_12, by=~Year, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)


CKD12_prev_year <- rbind(prev_year_2016, prev_year_2018, prev_year_2020, prev_year_2021, prev_year_2022, prev_year_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL)

#### Prevalence by sex--- ####
###ENSANUT 2016
prev_sex_2016<-svyby(~CKD_12, by=~Year+Sex2, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_sex_2018<-svyby(~CKD_12, by=~Year+Sex2, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_sex_2020<-svyby(~CKD_12, by=~Year+Sex2, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_sex_2021<-svyby(~CKD_12, by=~Year+Sex2, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

#ENSANUT 2022
prev_sex_2022<-svyby(~CKD_12, by=~Year+Sex2, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

#ENSANUT 2023
prev_sex_2023<-svyby(~CKD_12, by=~Year+Sex2, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)


CKD12_prev_sex <- rbind(prev_sex_2016, prev_sex_2018, prev_sex_2020, prev_sex_2021, prev_sex_2022, prev_sex_2023) %>% 
  as.data.frame() %>% mutate(group=rep(c("Women","Men"),6)) %>% `rownames<-`(NULL)


#### Prevalence by age--- ####
###ENSANUT 2016
prev_age_2016<-svyby(~CKD_12, by=~Year+Age_cat, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_age_2018<-svyby(~CKD_12, by=~Year+Age_cat, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_age_2020<-svyby(~CKD_12, by=~Year+Age_cat, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_age_2021<-svyby(~CKD_12, by=~Year+Age_cat, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

#ENSANUT 2022
prev_age_2022<-svyby(~CKD_12, by=~Year+Age_cat, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

#ENSANUT 2023
prev_age_2023<-svyby(~CKD_12, by=~Year+Age_cat, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

CKD12_prev_age <- rbind(
  prev_age_2016, prev_age_2018, prev_age_2020, prev_age_2021, prev_age_2022, prev_age_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL) %>% mutate(
    group=rep(c(1,2),6) %>% ordered(labels=c("<60 years","≥60 years")))

#### Prevalence by T2D--- #### 
###ENSANUT 2016
prev_dm_2016<-svyby(~CKD_12, by=~Year+diabetes, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_dm_2018<-svyby(~CKD_12, by=~Year+diabetes, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_dm_2020<-svyby(~CKD_12, by=~Year+diabetes, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_dm_2021<-svyby(~CKD_12, by=~Year+diabetes, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2022
prev_dm_2022<-svyby(~CKD_12, by=~Year+diabetes, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2023
prev_dm_2023<-svyby(~CKD_12, by=~Year+diabetes, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

CKD12_prev_dm <- rbind(
  prev_dm_2016, prev_dm_2018, prev_dm_2020, prev_dm_2021, prev_dm_2022, prev_dm_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL) %>% mutate(
    group=rep(c(1,2),6) %>% ordered(labels=c("No diabetes","Diabetes")))



#### Prevalence by HBP--- #### 
###ENSANUT 2016
prev_hbp_2016<-svyby(~CKD_12, by=~Year+hypertension, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_hbp_2018<-svyby(~CKD_12, by=~Year+hypertension, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_hbp_2020<-svyby(~CKD_12, by=~Year+hypertension, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_hbp_2021<-svyby(~CKD_12, by=~Year+hypertension, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2022
prev_hbp_2022<-svyby(~CKD_12, by=~Year+hypertension, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2023
prev_hbp_2023<-svyby(~CKD_12, by=~Year+hypertension, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_12*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2009") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)


CKD12_prev_hbp <- rbind(
  prev_hbp_2016, prev_hbp_2018, prev_hbp_2020, prev_hbp_2021,prev_hbp_2022, prev_hbp_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL) %>% mutate(
    group=rep(c(1,2),6) %>% ordered(labels=c("No hypertension","Hypertension")))


####---------------------------- CKD PREVALENCE: CKD-EPI 2021 ------------------#### ----####
#### Prevalence by year-- ####
###ENSANUT 2016
prev_year_2016_2<-svyby(~CKD_21, by=~Year, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_year_2018_2<-svyby(~CKD_21, by=~Year, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_year_2020_2<-svyby(~CKD_21, by=~Year, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_year_2021_2<-svyby(~CKD_21, by=~Year, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2022
prev_year_2022_2<-svyby(~CKD_21, by=~Year, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2023
prev_year_2023_2<-svyby(~CKD_21, by=~Year, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)


CKD21_prev_year <- rbind(prev_year_2016_2, prev_year_2018_2, prev_year_2020_2, prev_year_2021_2, prev_year_2022_2, prev_year_2023_2) %>%
  as.data.frame() %>% `rownames<-`(NULL)


#### Prevalence by sex--- ####
###ENSANUT 2016
prev_sex_2016_2<-svyby(~CKD_21, by=~Year+Sex2, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_sex_2018_2<-svyby(~CKD_21, by=~Year+Sex2, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_sex_2020_2<-svyby(~CKD_21, by=~Year+Sex2, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_sex_2021_2<-svyby(~CKD_21, by=~Year+Sex2, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2022
prev_sex_2022_2<-svyby(~CKD_21, by=~Year+Sex2, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2023
prev_sex_2023_2<-svyby(~CKD_21, by=~Year+Sex2, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)


CKD21_prev_sex <- rbind(prev_sex_2016_2, prev_sex_2018_2, prev_sex_2020_2, prev_sex_2021_2, prev_sex_2022_2, prev_sex_2023_2) %>% as.data.frame() %>%
  mutate(group=rep(c("Women","Men"),6)) %>% `rownames<-`(NULL)


#### Prevalence by age--- ####
###ENSANUT 2016
prev_age_2016_2<-svyby(~CKD_21, by=~Year+Age_cat, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_age_2018_2<-svyby(~CKD_21, by=~Year+Age_cat, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_age_2020_2<-svyby(~CKD_21, by=~Year+Age_cat, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_age_2021_2<-svyby(~CKD_21, by=~Year+Age_cat, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2022
prev_age_2022_2<-svyby(~CKD_21, by=~Year+Age_cat, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2023
prev_age_2023_2<-svyby(~CKD_21, by=~Year+Age_cat, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)


CKD21_prev_age <- rbind(
  prev_age_2016_2, prev_age_2018_2, prev_age_2020_2, prev_age_2021_2,
  prev_age_2022_2, prev_age_2023_2) %>% as.data.frame() %>% `rownames<-`(NULL) %>% mutate(
    group=rep(c(1,2),6) %>% ordered(labels=c("<60 years","≥60 years")))

#### Prevalence by T2D--- #### 
###ENSANUT 2016
prev_dm_2016_2<-svyby(~CKD_21, by=~Year+diabetes, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_dm_2018_2<-svyby(~CKD_21, by=~Year+diabetes, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_dm_2020_2<-svyby(~CKD_21, by=~Year+diabetes, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_dm_2021_2<-svyby(~CKD_21, by=~Year+diabetes, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2022
prev_dm_2022_2<-svyby(~CKD_21, by=~Year+diabetes, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2023
prev_dm_2023_2<-svyby(~CKD_21, by=~Year+diabetes, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)


CKD21_prev_dm<- rbind(
  prev_dm_2016_2, prev_dm_2018_2, prev_dm_2020_2, prev_dm_2021_2,
  prev_dm_2022_2, prev_dm_2023_2) %>% as.data.frame() %>% `rownames<-`(NULL) %>% mutate(
    group=rep(c(1,2),6) %>% ordered(labels=c("No diabetes","Diabetes")))

#### Prevalence by HBP--- #### 
###ENSANUT 2016
prev_hbp_2016_2<-svyby(~CKD_21, by=~Year+hypertension, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_hbp_2018_2<-svyby(~CKD_21, by=~Year+hypertension, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_hbp_2020_2<-svyby(~CKD_21, by=~Year+hypertension, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_hbp_2021_2<-svyby(~CKD_21, by=~Year+hypertension, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2022
prev_hbp_2022_2<-svyby(~CKD_21, by=~Year+hypertension, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2023
prev_hbp_2023_2<-svyby(~CKD_21, by=~Year+hypertension, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)


CKD21_prev_hbp <- rbind(
  prev_hbp_2016_2, prev_hbp_2018_2, prev_hbp_2020_2, prev_hbp_2021_2,
  prev_hbp_2022_2, prev_hbp_2023_2) %>% as.data.frame() %>% `rownames<-`(NULL) %>% mutate(
    group=rep(c(1,2),6) %>% ordered(labels=c("No hypertension","Hypertension")))

#### Prevalence by Obesity--- #### 
###ENSANUT 2016
prev_obes_2016_2<-svyby(~CKD_21, by=~Year+BMI_cat, design=ensanut_2016_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2018
prev_obes_2018_2<-svyby(~CKD_21, by=~Year+BMI_cat, design=ensanut_2018_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2020
prev_obes_2020_2<-svyby(~CKD_21, by=~Year+BMI_cat, design=ensanut_2020_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2021
prev_obes_2021_2<-svyby(~CKD_21, by=~Year+BMI_cat, design=ensanut_2021_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2022
prev_obes_2022_2<-svyby(~CKD_21, by=~Year+BMI_cat, design=ensanut_2022_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)

###ENSANUT 2023
prev_obes_2023_2<-svyby(~CKD_21, by=~Year+BMI_cat, design=ensanut_2023_survey, svymean, na.rm=T) %>%  
  mutate(prop=round(CKD_21*100, digits=1), IC95=round((se*1.96)*100, digits=2)) %>%
  mutate(lIC95=prop-IC95,uIC95=prop+IC95, cluster="CKD-EPI 2021") %>% dplyr::select(prop,IC95,lIC95,uIC95,cluster,Year)


CKD21_prev_obes <- rbind(
  prev_obes_2016_2, prev_obes_2018_2, prev_obes_2020_2, prev_obes_2021_2,
  prev_obes_2022_2, prev_obes_2023_2) %>% as.data.frame() %>% `rownames<-`(NULL) %>% mutate(
    group=rep(c(1,2),6) %>% ordered(labels=c("No obesity","Obesity")))


####----------- Total number of individuals with CKD -----------####
#### CKD-EPI 2009 ------- ####
total_2016_eq09 <- wt_total("CKD_12", ensanut_2016_survey, 2016)
total_2018_eq09 <- wt_total("CKD_12", ensanut_2018_survey, 2018)
total_2020_eq09 <- wt_total("CKD_12", ensanut_2020_survey, 2020)
total_2021_eq09 <- wt_total("CKD_12", ensanut_2021_survey, 2021)
total_2022_eq09 <- wt_total("CKD_12", ensanut_2022_survey, 2022)
total_2023_eq09 <- wt_total("CKD_12", ensanut_2023_survey, 2023)

total_eq09 <- rbind(total_2016_eq09, total_2018_eq09, total_2020_eq09, total_2021_eq09, total_2022_eq09, total_2023_eq09)
total_eq09$eq <- "CKD-EPI 2009"

#### CKD-EPI 2021 ------- ####
total_2016_eq21 <- wt_total("CKD_21", ensanut_2016_survey, 2016)
total_2018_eq21 <- wt_total("CKD_21", ensanut_2018_survey, 2018)
total_2020_eq21 <- wt_total("CKD_21", ensanut_2020_survey, 2020)
total_2021_eq21 <- wt_total("CKD_21", ensanut_2021_survey, 2021)
total_2022_eq21 <- wt_total("CKD_21", ensanut_2022_survey, 2022)
total_2023_eq21 <- wt_total("CKD_21", ensanut_2023_survey, 2023)

total_eq21 <- rbind(total_2016_eq21, total_2018_eq21, total_2020_eq21, total_2021_eq21, total_2022_eq21, total_2023_eq21)
total_eq21$eq <- "CKD-EPI 2021"


#### CKD total number table -- ####
total_ckd_number <- rbind(total_eq09, total_eq21) 
total_labels <- total_ckd_number %>% 
  mutate(label = paste0(sprintf(total, fmt = "%#.2f"), " (", sprintf(lower, fmt = "%#.2f"), "-", sprintf(upper, fmt = "%#.2f"), ")")) %>% 
  select(year, eq, label) %>% 
  pivot_wider(names_from = eq,
              values_from = label) %>% 
  relocate(`CKD-EPI 2021`, .before =  `CKD-EPI 2009`)

total_diff <- total_ckd_number %>% 
  select(year, eq, total) %>% 
  pivot_wider(names_from = eq,
              values_from = total) %>% 
  mutate(diff = `CKD-EPI 2021` - `CKD-EPI 2009`)

cbind(total_labels, total_diff$diff) %>% 
  gt() %>% 
  gtsave(filename = "tabla_ckd_diff.docx")

####----------- eGFR and CKD Reclassification ------------ #####
#### Descriptive characteristics ---- ####
ensanut_2023_survey %>% 
  tbl_svysummary(include = c(edad, Age_cat, Creatinine1, egfr_12, egfr_21, Sex2, DM2, HBP),
                 by = reclass_egfr,
                 statistic = list(all_categorical() ~ "{p}%",
                                  all_continuous() ~ "{mean}"),
                 label = list(edad ~ "Mean age",
                              Age_cat ~ "Age category",
                              Creatinine1 ~ "Creatinine",
                              Sex2 ~ "Sex",
                              DM2 ~ "Diabetes",
                              HBP ~ "Hypertension",
                              egfr_12 ~ "eGFRcr 2009",
                              egfr_21 ~ "eGFRcr 2021"),
                 missing = "no") %>% 
  add_ci(pattern = "{stat} ({ci})") %>% 
  as_gt() %>% 
  gt::gtsave(filename = "table_reclass_egfr.docx")

ensanut_2023_survey %>% 
  tbl_svysummary(include = c(edad, Age_cat, Creatinine1, egfr_12, egfr_21, Sex2, DM2, HBP),
                 by = reclass_ckd,
                 statistic = list(all_categorical() ~ "{p}%",
                                  all_continuous() ~ "{mean}"),
                 label = list(edad ~ "Mean age",
                              Age_cat ~ "Age category",
                              Creatinine1 ~ "Creatinine",
                              Sex2 ~ "Sex",
                              DM2 ~ "Diabetes",
                              HBP ~ "Hypertension",
                              egfr_12 ~ "eGFRcr 2009",
                              egfr_21 ~ "eGFRcr 2021"),
                 missing = "no") %>% 
  add_ci(pattern = "{stat} ({ci})") %>% 
  as_gt() %>% 
  gt::gtsave(filename = "table_reclass_ckd.docx")


#### eGFR distribution -------------- #####
median_egfr_09 <- svyquantile(~ egfr_12, ensanut_2023_survey, quantiles = c(0.25, 0.5, 0.75))
median_egfr_21 <- svyquantile(~ egfr_21, ensanut_2023_survey, quantiles = c(0.25, 0.5, 0.75))

egfr_09 <- svysmooth(~ egfr_12, ensanut_2023_survey)
egfr_21 <- svysmooth(~ egfr_21, ensanut_2023_survey)

df1 <- data.frame(x_09 = egfr_09$egfr_12$x,
                  y_09 = egfr_09$egfr_12$y,
                  x_21 = egfr_21$egfr_21$x,
                  y_21 = egfr_21$egfr_21$y)

fig_dens <- df1 %>% ggplot(aes()) +
  geom_line(aes(x = x_09, y = y_09, color = "2009 CKD-EPI equation")) +
  geom_line(aes(x = x_21, y = y_21, color = "2021 CKD-EPI equation")) +
  geom_vline(xintercept = c(15,30,45,60,90), color = "gray", linetype = "dashed") +
  scale_x_continuous(breaks = c(15,30,45,60, 90)) +
  scale_color_manual(values = c("2009 CKD-EPI equation" = "#3172B1", "2021 CKD-EPI equation" = "#AE4532")) +
  labs(x = "eGFR",
       y = "Density",
       color = "") +
  theme_classic()
  
fig_dens_ckd <- df1 %>% ggplot(aes()) +
  geom_line(aes(x = x_09, y = y_09, color = "2009 CKD-EPI equation")) +
  geom_line(aes(x = x_21, y = y_21, color = "2021 CKD-EPI equation")) +
  geom_vline(xintercept = c(15,30,45,60,90), color = "gray", linetype = "dashed") +
  scale_x_continuous(breaks = c(15,30,45,60), limits = c(0,60)) +
  scale_y_continuous(limits = c(0, 0.004)) +
  scale_color_manual(values = c("2009 CKD-EPI equation" = "#3172B1", "2021 CKD-EPI equation" = "#AE4532")) +
  labs(x = "eGFR",
       y = "Density",
       color = "") +
  theme_classic()

fig_dens_final <- ggarrange(fig_dens, fig_dens_ckd, nrow = 1, ncol = 2, 
                            labels = "AUTO", legend = "bottom", common.legend = TRUE)

ggsave(fig_dens_final, file="fig_dens.png", bg="transparent",
       width = 20, height = 8, units = c("cm"), dpi = 600, limitsize = FALSE)

####---------------------------- TOTAL NUMBER + PREVALENCE OF CKD ---------------#### ----####
#### T: Prior CKD diagnosis------ ####
###ENSANUT 2016
total_2016 <- svytotal(~HX_CKD, ensanut_2016_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2016)

###ENSANUT 2018
total_2018 <- svytotal(~HX_CKD, ensanut_2018_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2018)

###ENSANUT 2021
total_2021 <- svytotal(~HX_CKD, ensanut_2021_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2021)

###ENSANUT 2022
total_2022 <- svytotal(~HX_CKD, ensanut_2022_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2022)

###ENSANUT 2023
total_2023 <- svytotal(~HX_CKD, ensanut_2023_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2023)


TOTALCKD.1 <- rbind(total_2016, total_2018, total_2021, total_2022, total_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL)

tab_ckd_total1 <- TOTALCKD.1 %>% 
  mutate(label = paste0(sprintf(tot, fmt = "%#.2f") , " (", sprintf(lower, fmt = "%#.2f"), "-", sprintf(upper, fmt = "%#.2f"), ")")) %>% 
  select(Year, label) %>% 
  pivot_wider(names_from = Year, values_from = label) %>%
  mutate(type = "Prior CKD diagnosis", .before = `2016`)

#### T: Prior DX + CKD-EPI 2009-- ####
###ENSANUT 2016
total_2016 <- svytotal(~CKD_12.B, ensanut_2016_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2009)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2016)

###ENSANUT 2018
total_2018 <- svytotal(~CKD_12.B, ensanut_2018_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2009)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2018)

###ENSANUT 2021
total_2021 <- svytotal(~CKD_12.B, ensanut_2021_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2009)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2021)

###ENSANUT 2022
total_2022 <- svytotal(~CKD_12.B, ensanut_2022_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2009)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2022)

###ENSANUT 2023
total_2023 <- svytotal(~CKD_12.B, ensanut_2023_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2009)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2023)


TOTALCKD.2 <- rbind(total_2016, total_2018, total_2021, total_2022, total_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL)

tab_ckd_total2 <- TOTALCKD.2 %>% 
  mutate(label = paste0(sprintf(tot, fmt = "%#.2f") , " (", sprintf(lower, fmt = "%#.2f"), "-", sprintf(upper, fmt = "%#.2f"), ")")) %>% 
  select(Year, label) %>% 
  pivot_wider(names_from = Year, values_from = label) %>%
  mutate(type = "Prior DX + CKD-EPI 2009", .before = `2016`)

#### T: Prior DX + CKD-EPI 2021-- ####
###ENSANUT 2016
total_2016 <- svytotal(~CKD_21.B, ensanut_2016_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2021)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2016)

###ENSANUT 2018
total_2018 <- svytotal(~CKD_21.B, ensanut_2018_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2021)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2018)

###ENSANUT 2021
total_2021 <- svytotal(~CKD_21.B, ensanut_2021_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2021)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2021)

###ENSANUT 2022
total_2022 <- svytotal(~CKD_21.B, ensanut_2022_survey) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2021)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2022)

###ENSANUT 2023
total_2023 <- svytotal(~CKD_21.B, ensanut_2023_survey) %>%
  as_tibble() %>% `names<-`(c("tot","se")) %>% mutate(
    tot=round(tot, digits=1), "IC95"=(se*1.96)) %>%
  mutate("lower"=tot-IC95, "upper"=tot+IC95,
         cluster="+ decreased eGFR\n(CKD-EPI 2021)") %>%
  dplyr::select(tot,IC95,lower,upper,cluster) %>% mutate("Year"=2023)

TOTALCKD.3 <- rbind(total_2016, total_2018, total_2021, total_2022, total_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL)

tab_ckd_total3 <- TOTALCKD.3 %>% 
  mutate(label = paste0(sprintf(tot, fmt = "%#.2f") , " (", sprintf(lower, fmt = "%#.2f"), "-", sprintf(upper, fmt = "%#.2f"), ")")) %>% 
  select(Year, label) %>% 
  pivot_wider(names_from = Year, values_from = label) %>%
  mutate(type = "Prior DX + CKD-EPI 2021", .before = `2016`)

tab_ckd_total <- rbind(tab_ckd_total1, tab_ckd_total2, tab_ckd_total3)
tab_ckd_total %>%  gt() %>% 
  gtsave(filename = "tabla_ckd_diag_total.docx")


#### P: Prior CKD diagnosis------ ####
###ENSANUT 2016
prev2_2016 <- svymean(~HX_CKD, ensanut_2016_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2016) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2018
prev2_2018 <- svymean(~HX_CKD, ensanut_2018_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2018) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2021
prev2_2021 <- svymean(~HX_CKD, ensanut_2021_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2021) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2022
prev2_2022 <- svymean(~HX_CKD, ensanut_2022_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2022) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2023
prev2_2023 <- svymean(~HX_CKD, ensanut_2023_survey_h, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="Prior diagnosis") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2023) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

PREVCKD.1 <- rbind(prev2_2016, prev2_2018, prev2_2021, prev2_2022, prev2_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL)

tab_ckd_diag1 <- PREVCKD.1 %>% 
  mutate(label = paste0(sprintf(prev, fmt = "%#.2f"), "%" , " (", sprintf(lower, fmt = "%#.2f"), "-", sprintf(upper, fmt = "%#.2f"), ")")) %>% 
  select(Year, label) %>% 
  pivot_wider(names_from = Year, values_from = label) %>%
  mutate(type = "Prior CKD diagnosis", .before = `2016`)

#### P: Prior DX + CKD-EPI 2009-- ####
###ENSANUT 2016
prev2_2016 <- svymean(~CKD_12.B, ensanut_2016_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate(
    "IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2012") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2016) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2018
prev2_2018 <- svymean(~CKD_12.B, ensanut_2018_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2012") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2018) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2021
prev2_2021 <- svymean(~CKD_12.B, ensanut_2021_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2012") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2021) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2022
prev2_2022 <- svymean(~CKD_12.B, ensanut_2022_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2012") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2022) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2023
prev2_2023 <- svymean(~CKD_12.B, ensanut_2023_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2012") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2023) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

PREVCKD.2 <- rbind(prev2_2016, prev2_2018, prev2_2021, prev2_2022, prev2_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL)

tab_ckd_diag2 <- PREVCKD.2 %>% 
  mutate(label = paste0(sprintf(prev, fmt = "%#.2f"), "%" , " (", sprintf(lower, fmt = "%#.2f"), "-", sprintf(upper, fmt = "%#.2f"), ")")) %>% 
  select(Year, label) %>% 
  pivot_wider(names_from = Year, values_from = label) %>% 
  mutate(type = "Prior DX + CKD-EPI 2009", .before = `2016`)
  

#### P: Prior DX + CKD-EPI 2021-- ####
###ENSANUT 2016
prev2_2016 <- svymean(~CKD_21.B, ensanut_2016_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2021") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2016) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2018
prev2_2018 <- svymean(~CKD_21.B, ensanut_2018_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2021") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2018) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2021
prev2_2021 <- svymean(~CKD_21.B, ensanut_2021_survey, na.rm=T) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2021") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2021) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2022
prev2_2022 <- svymean(~CKD_21.B, ensanut_2022_survey) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2021") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2022) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

###ENSANUT 2023
prev2_2023 <- svymean(~CKD_21.B, ensanut_2023_survey) %>%
  as_tibble() %>% `names<-`(c("prev","se")) %>% mutate("IC95"=(se*1.96)) %>%
  mutate("lower"=prev-IC95, "upper"=prev+IC95, cluster="CKD-EPI 2021") %>%
  dplyr::select(prev,IC95,lower,upper,cluster) %>% mutate("Year"=2023) %>%
  mutate(prev=prev*100,IC95=IC95*100,lower=lower*100,upper=upper*100)

PREVCKD.3 <- rbind(prev2_2016, prev2_2018, prev2_2021, prev2_2022, prev2_2023) %>%
  as.data.frame() %>% `rownames<-`(NULL)

tab_ckd_diag3 <- PREVCKD.3 %>% 
  mutate(label = paste0(sprintf(prev, fmt = "%#.2f"), "%" , " (", sprintf(lower, fmt = "%#.2f"), "-", sprintf(upper, fmt = "%#.2f"), ")")) %>% 
  select(Year, label) %>% 
  pivot_wider(names_from = Year, values_from = label) %>% 
  mutate(type = "Prior DX + CKD-EPI 2021", .before = `2016`)

tab_ckd_diag <- rbind(tab_ckd_diag1, tab_ckd_diag2, tab_ckd_diag3)
tab_ckd_diag %>%  gt() %>% 
  gtsave(filename = "tabla_ckd_diag_prev.docx")

####---------------------------- FIGURES ---------------------------------------#### ----####
#### Figure 1a: eGFR prevalence ---------- ####
fig1a <- egfr_21_prev_total %>% 
  filter(cat == "G1") %>% 
  ggplot(aes(x = factor(year), y = prop, ymin = lower, ymax = upper, color = cat)) +
  geom_pointrange() +
  geom_line(aes(group = cat)) +
  labs(x = "Year",
       y = "Prevalence",
       color = "") +
  scale_y_continuous(labels = scales::percent, limits = c(0.5, 1)) +
  scale_color_manual(values = c("#34733F")) +
  theme_classic()

fig1b <- egfr_21_prev_total %>% 
  filter(cat == "G2") %>% 
  ggplot(aes(x = factor(year), y = prop, ymin = lower, ymax = upper, color = cat)) +
  geom_pointrange() +
  geom_line(aes(group = cat)) +
  labs(x = "Year",
       y = "Prevalence",
       color = "") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.20)) +
  scale_color_manual(values = c("#4EAC5B")) +
  theme_classic()

fig1c <- egfr_21_prev_total %>% 
  filter(cat != "G1" & cat != "G2") %>% 
  ggplot(aes(x = factor(year), y = prop, ymin = lower, ymax = upper, color = cat)) +
  geom_pointrange(alpha = 0.8) +
  geom_line(aes(group = cat), alpha = 0.8) +
  labs(x = "Year",
       y = "Prevalence",
       color = "") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.03)) +
  scale_color_manual(values = c("#F8D353", "#E99B57", "#EB503A", "#AF3023")) +
  theme_classic()

fig1d <- egfr_21_prev_total %>% 
  ggplot(aes(x = factor(year), y = prop, ymin = lower, ymax = upper, color = cat)) +
  geom_pointrange() +
  labs(color = "eGFR Categories") +
  scale_color_manual(values = c("#34733F", "#4EAC5B","#F8D353", "#E99B57", "#EB503A", "#AF3023")) +
  theme_classic()

fig_1a <- ggarrange(fig1a, fig1b, fig1c, nrow = 1, ncol = 3, labels = "AUTO", legend = "bottom",
                   legend.grob = get_legend(fig1d, position = "bottom"))

#### Figure 1b: CKD prevalence ----------- ####
#Overall CKD prevalence
fig2a <- CKD21_prev_year %>% 
  mutate(ckd = "CKD",
         label = paste0(prop, "%")) %>% 
  ggplot(aes(x = factor(Year), y = prop/100, ymin = lIC95/100, ymax = uIC95/100, color = ckd)) +
  geom_pointrange() +
  geom_line(aes(group = ckd)) +
  labs(x = "Year",
       y = "Prevalence",
       color = "") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.05)) +
  scale_color_manual(values = c("black")) +
  theme_classic()

#CKD prevalence by sex
fig2b <- CKD21_prev_sex %>% 
  ggplot(aes(x = factor(Year), y = prop/100, ymin = lIC95/100, ymax = uIC95/100, color = group)) +
  geom_pointrange() +
  geom_line(aes(group = group)) +
  labs(x = "Year",
       y = "Prevalence",
       color = "") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.2)) +
  scale_color_manual(values = c("#99292A", "#5388B8")) +
  theme_classic()

#CKD prevalence by age
fig2c <- CKD21_prev_age %>% 
  mutate(lIC95 = if_else(lIC95 < 0, 0, lIC95)) %>% 
  ggplot(aes(x = factor(Year), y = prop/100, ymin = lIC95/100, ymax = uIC95/100, color = group)) +
  geom_pointrange() +
  geom_line(aes(group = group)) +
  labs(x = "Year",
       y = "Prevalence",
       color = "") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.2)) +
  scale_color_manual(values = c("#99292A", "#5388B8")) +
  theme_classic()

#CKD prevalence by diabetes
fig2d <- CKD21_prev_dm %>% 
  ggplot(aes(x = factor(Year), y = prop/100, ymin = lIC95/100, ymax = uIC95/100, color = group)) +
  geom_pointrange() +
  geom_line(aes(group = group)) +
  labs(x = "Year",
       y = "Prevalence",
       color = "") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.2)) +
  scale_color_manual(values = c("#99292A", "#5388B8")) +
  theme_classic()

#CKD prevalence by hypertension
fig2e <- CKD21_prev_hbp %>% 
  ggplot(aes(x = factor(Year), y = prop/100, ymin = lIC95/100, ymax = uIC95/100, color = group)) +
  geom_pointrange() +
  geom_line(aes(group = group)) +
  labs(x = "Year",
       y = "Prevalence",
       color = "") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.2)) +
  scale_color_manual(values = c("#99292A", "#5388B8")) +
  theme_classic()

#CKD prevalence by obesity
fig2f <- CKD21_prev_obes %>% 
  ggplot(aes(x = factor(Year), y = prop/100, ymin = lIC95/100, ymax = uIC95/100, color = group)) +
  geom_pointrange() +
  geom_line(aes(group = group)) +
  labs(x = "Year",
       y = "Prevalence",
       color = "") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.2)) +
  scale_color_manual(values = c("#99292A", "#5388B8")) +
  theme_classic()

fig_1b <- ggarrange(fig2a, fig2b, fig2c, fig2d, fig2e, fig2f, nrow = 2, ncol = 3, 
                   labels = c(LETTERS[4:9]), legend = "bottom")

fig_1 <- ggarrange(fig_1a, fig_1b, nrow = 2, ncol = 1, heights = c(0.9,1.7))

ggsave(fig_1, file="Figure1.pdf", bg="transparent",
       width = 26, height = 28, units = c("cm"), dpi = 600, limitsize = FALSE)


#### Figure 2: Comparison of CKD prevalence between equations####
#Overall
CKD12_prev_year <- rbind(prev_year_2016, prev_year_2018, prev_year_2020, 
                         prev_year_2021, prev_year_2022, prev_year_2023)
CKD21_prev_year <- rbind(prev_year_2016_2, prev_year_2018_2, prev_year_2020_2, 
                         prev_year_2021_2, prev_year_2022_2, prev_year_2023_2)
ckd_overall <- rbind(CKD12_prev_year, CKD21_prev_year)

fig2a <- ckd_overall %>% 
  mutate(prop = prop/100, lower = lIC95/100, upper = uIC95/100) %>% 
  ggplot(aes(x = factor(Year), y = prop, ymin = lower, ymax = upper, color = cluster)) +
  geom_pointrange(alpha = 0.8, show.legend = FALSE) +
  geom_line(aes(group = cluster), alpha = 0.8) +
  labs(x = "Year",
       y = "Prevalence",
       color = "Equation") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.07)) +
  scale_color_manual(values = c("#3172B1", "#AE4532")) +
  theme_classic()

#Sex
ckd_sex_total <- rbind(CKD12_prev_sex, CKD21_prev_sex)

fig2b <- ckd_sex_total %>%
  mutate(prop = prop/100, lower = lIC95/100, upper = uIC95/100) %>% 
  ggplot(aes(x = factor(Year), y = prop, ymin = lower, ymax = upper, color = cluster)) +
  geom_pointrange(alpha = 0.8, show.legend = FALSE) +
  geom_line(aes(group = cluster), alpha = 0.8) +
  labs(x = "Year",
       y = "Prevalence",
       color = "Equation") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.07)) +
  scale_color_manual(values = c("#3172B1", "#AE4532")) +
  facet_wrap(~ group) +
  theme_classic()

#Age
ckd_age_total <- rbind(CKD12_prev_age, CKD21_prev_age)

fig2c <- ckd_age_total %>%
  mutate(prop = prop/100, lower = lIC95/100, upper = uIC95/100,
         lower = if_else(lower < 0, 0, lower)) %>% 
  ggplot(aes(x = factor(Year), y = prop, ymin = lower, ymax = upper, color = cluster)) +
  geom_pointrange(alpha = 0.8, show.legend = FALSE) +
  geom_line(aes(group = cluster), alpha = 0.8) +
  labs(x = "Year",
       y = "Prevalence",
       color = "Equation") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  scale_color_manual(values = c("#3172B1", "#AE4532")) +
  facet_wrap(~ group) +
  theme_classic()

#Diabetes
ckd_diab_total <- rbind(CKD12_prev_dm, CKD21_prev_dm)

fig2d <- ckd_diab_total %>%
  mutate(prop = prop/100, lower = lIC95/100, upper = uIC95/100) %>% 
  ggplot(aes(x = factor(Year), y = prop, ymin = lower, ymax = upper, color = cluster)) +
  geom_pointrange(alpha = 0.8, show.legend = FALSE) +
  geom_line(aes(group = cluster), alpha = 0.8) +
  labs(x = "Year",
       y = "Prevalence",
       color = "Equation") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  scale_color_manual(values = c("#3172B1", "#AE4532")) +
  facet_wrap(~ group) +
  theme_classic()

#Diabetes
ckd_hbp_total <- rbind(CKD12_prev_hbp, CKD21_prev_hbp)

fig2e <- ckd_hbp_total %>%
  mutate(prop = prop/100, lower = lIC95/100, upper = uIC95/100) %>% 
  ggplot(aes(x = factor(Year), y = prop, ymin = lower, ymax = upper, color = cluster)) +
  geom_pointrange(alpha = 0.8, show.legend = FALSE) +
  geom_line(aes(group = cluster), alpha = 0.8) +
  labs(x = "Year",
       y = "Prevalence",
       color = "Equation") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.25)) +
  scale_color_manual(values = c("#3172B1", "#AE4532")) +
  facet_wrap(~ group) +
  theme_classic()

fig_2_pre <- ggarrange(fig_dens, fig2a, nrow = 1, ncol = 2, 
                      labels = c("A", "B"), legend = "none", common.legend = TRUE)

fig_2_pre_2 <- ggarrange(fig2b, fig2c, fig2d, fig2e, nrow = 2, ncol = 2, 
                   labels = c("C", "D", "E", "F"), 
                   common.legend = TRUE,
                   heights = c( 1, 1, 1, 1), legend = "bottom")
fig_2<-ggarrange(fig_2_pre,fig_2_pre_2, heights = c(0.33, 0.66), nrow=2,ncol=1,common.legend = T, legend="bottom")

ggsave(fig_2, file="Figure2.pdf", width = 24, height = 24, units = c("cm"), dpi = 600, limitsize = FALSE)

#### Figure 3: MCPS estimations ####
# Load all data in ckd_mcps.R #
#eGFR distribution
fig3a <- mcps_fin %>% ggplot(aes()) + 
  geom_density(aes(x = EFGR.09, color = "CKD-EPI 2009"), show.legend = FALSE) +
  geom_density(aes(x = EFGR.21, color = "CKD-EPI 2021"), show.legend = FALSE) +
  labs(x = "eGFR",
       y = "Density",
       color = "") +
  geom_vline(xintercept = c(15,30,45,60,90), color = "gray", linetype = "dashed") +
  scale_x_continuous(breaks = c(15,30,45,60, 90)) +
  scale_color_manual(values = c("#3172B1", "#AE4532")) +
  theme_classic()

#HR comparison
fig3b <- df %>% ggplot(aes()) +
  geom_line(aes(x = x_09, y = y_09, color = "2009 CKD-EPI equation")) +
  geom_ribbon(aes(x = x_09, y = y_09, ymin = lower_09, ymax = upper_09), 
              alpha = 0.15, show.legend = FALSE, color = NA, fill = "#3172B1") +
  geom_line(aes(x = x_21, y = y_21, color = "2021 CKD-EPI equation")) +
  geom_ribbon(aes(x = x_21, y = y_21, ymin = lower_21, ymax = upper_21), 
              alpha = 0.15, show.legend = FALSE, color = NA, fill = "#AE4532") +
  geom_vline(xintercept = c(15,30,45,60,90), color = "gray", linetype = "dashed") +
  geom_hline(yintercept = 1, color = "gray", linetype = "dashed") +
  scale_x_continuous(breaks = c(15,30,45,60, 90)) +
  scale_y_log10() +
  scale_color_manual(values = c("2009 CKD-EPI equation" = "#3172B1", "2021 CKD-EPI equation" = "#AE4532")) +
  labs(x = "eGFR",
       y = "HR (95% CI)",
       color = "") +
  theme_classic()

#Concordance figures
fig3c <- cstat_final %>% filter(year == "5-year") %>% 
  mutate(color_group = ifelse(equation == "2009 equation", 
                              factor(type), 
                              paste(factor(type), equation))) %>%
  ggplot(aes(x = time, y = coef)) +
  geom_pointrange(aes(ymax = upper, ymin = lower, color = color_group), show.legend = FALSE) +
  geom_line(aes(color = color_group), show.legend = FALSE) +
  labs(x = "Years",
       y = "Concordance",
       color = "") +
  ylim(0, 1) +
  facet_grid(rows = vars(equation)) +
  scale_color_manual(values = c("#616463", "#3172B1", "#616463", "#AE4532")) +
  theme_classic()

fig3d <- cstat_final %>% filter(year == "10-year") %>% 
  mutate(color_group = ifelse(equation == "2009 equation", 
                              factor(type), 
                              paste(factor(type), equation))) %>%
  ggplot(aes(x = time, y = coef)) +
  geom_pointrange(aes(ymax = upper, ymin = lower, color = color_group), show.legend = FALSE) +
  geom_line(aes(color = color_group), show.legend = FALSE) +
  labs(x = "Years",
       y = "Concordance",
       color = "") +
  ylim(0, 1) +
  facet_grid(rows = vars(equation)) +
  scale_color_manual(values = c("#616463", "#3172B1", "#616463", "#AE4532")) +
  theme_classic()

fig_3 <- ggarrange(fig3a, fig3b, fig3c, fig3d, nrow = 2, ncol = 2, labels = "AUTO", legend = "bottom",
                    legend.grob = get_legend(fig3b, position = "bottom"))

ggsave(fig_3, file="Figure3.pdf", bg="transparent",
       width = 22, height = 16, units = c("cm"), dpi = 600, limitsize = FALSE)


#### Supplementary Figure 1: Total by eGFR + prior diagnosis ####
#Prevalence by eGFR alone
fnc1 <- function(x,y,z){
  paste0(sprintf(x/1000,fmt="%#.1f"),"\n(",
         sprintf(y/1000,fmt="%#.1f"),", ",
         sprintf(z/1000,fmt="%#.1f"),")")}
fnc2 <- function(x,y,z){
  paste0(sprintf(x,fmt="%#.1f"),"%","\n(",
         sprintf(y,fmt="%#.1f"),", ",
         sprintf(z,fmt="%#.1f"),")")}

CKD12_prev_year <- rbind(prev_year_2016, prev_year_2018,
                         prev_year_2020, prev_year_2021, prev_year_2022)
CKD21_prev_year <- rbind(prev_year_2016_2, prev_year_2018_2,
                         prev_year_2020_2, prev_year_2021_2, prev_year_2022_2)
combined_data <- rbind(mutate(CKD12_prev_year, Formula = "CKD-EPI 2009"),
                       mutate(CKD21_prev_year, Formula = "CKD-EPI 2021")) %>%
  mutate("Year"=as.factor(Year), "label"=fnc2(prop, lIC95, uIC95),
         "y.pos"=ifelse(Formula == "CKD-EPI 2009", 0.2, -0.2))

val.min <- ((c(combined_data$lIC95) %>% min*100)%>% floor)/100
val.max <- ((c(combined_data$uIC95) %>% max*100)%>% ceiling)/100

#Prevalence and total by eGFR + prior diagnosis
combined_data2 <- rbind(TOTALCKD.1, TOTALCKD.2, TOTALCKD.3) %>%
  select(Year,cluster,tot,lower,upper) %>%
  rename(tot.l=lower,tot.u=upper, year=Year) %>% cbind(
    rbind(PREVCKD.1, PREVCKD.2, PREVCKD.3) %>% select(
      prev,lower,upper) %>% rename(prev.l=lower,prev.u=upper)) %>%
  mutate("year"=as.factor(year), "label1"=fnc1(tot, tot.l, tot.u),
         "label2"=fnc2(prev, prev.l, prev.u),
         "y.pos"=ifelse(cluster == "CKD-EPI 2009", 0.125, -0.125))

Fig1B <- combined_data2 %>% mutate(
  "tot"=tot/1000,"tot.l"=tot.l/1000,"tot.u"=tot.u/1000, "year"=ordered(
    year, levels=c("2023", "2022","2021","2018","2016")),"cluster"=ordered(
      cluster, levels=c("+ decreased eGFR\n(CKD-EPI 2021)",
                        "+ decreased eGFR\n(CKD-EPI 2009)",
                        "Prior diagnosis")),) %>% 
  ggplot(aes(x=year, y=tot, fill=cluster), xLabels=NA) + geom_col(
    color="black", linetype=1, position=position_dodge2(width=1)) + 
  geom_errorbar(aes(ymin = tot.l, ymax = tot.u), col="black",
                position = position_dodge2(width = 1, padding = 0.65)) +
  labs(fill="") + ggpubr::theme_pubclean() +
  xlab("ENSANUT cycle") + ylab ("Total number (in thousands)") +
  scale_fill_manual(values=c("#392220","#CA4E49","#D49797"),
                    guide=guide_legend(reverse = T)) + 
  ggtitle("Mexicans with CKD (by prior diagnosis or decreased eGFR)") +
  coord_flip() + geom_text(
    aes(label=label1, y=tot.u+100), size=4, col="black", fontface="bold",
            position=position_dodge2(width=1), hjust=0) + ylim(0,5000) +
  theme(legend.position = "bottom", plot.title = element_text(
    hjust=0.5, vjust=-2, size=17, face = "bold")); Fig1B

ggsave(Fig1B, file="fig_supp.png", bg="white",
       width=32, height=24, units=c("cm"), dpi=600, limitsize = FALSE)


####---------------------------- TABLES ----------------------------------------#### ----####
#### Table 1: Prevalence of CKD by year -------------------- ####
#Decreased eGFR
fnc3 <- function(x,y,z){
  paste0(sprintf(x,fmt="%#.1f"),"%"," (",
         sprintf(y,fmt="%#.1f"),"% - ",
         sprintf(z,fmt="%#.1f"),"%)")}
tab_prev_12 <- rbind(
  CKD12_prev_year %>% mutate(group="Overall"), CKD12_prev_sex,
  CKD12_prev_age, CKD12_prev_dm, CKD12_prev_hbp) %>% mutate(
    "prev_f"=fnc3(prop, lIC95, uIC95))
tab_prev_21 <- rbind(
  CKD21_prev_year %>% mutate(group="Overall"), CKD21_prev_sex,
  CKD21_prev_age, CKD21_prev_dm, CKD21_prev_hbp) %>% mutate(
    "prev_f"=fnc3(prop, lIC95, uIC95))
tab_prev_12[1:5,] %>% transmute(Year, "C09"=prev_f) %>%
  cbind(tab_prev_21[1:5,] %>% transmute("C21"=prev_f)) %>%
  `rownames<-`(NULL) %>% mutate("Prior"="-") %>%
  select(Year, Prior, C09, C21) -> tab1A

#Prior diagnosis or decreased eGFR
tab_prev_pr.B <- PREVCKD.1 %>% mutate(group="Overall") %>%
  mutate("prev_f"=fnc3(prev, lower, upper))
tab_prev_12.B <- PREVCKD.2 %>% mutate(group="Overall") %>%
  mutate("prev_f"=fnc3(prev, lower, upper))
tab_prev_21.B <- PREVCKD.3 %>% mutate(group="Overall") %>%
  mutate("prev_f"=fnc3(prev, lower, upper))
tab_prev_pr.B %>% transmute(Year, "Prior"=prev_f) %>%
  cbind(tab_prev_12.B %>% transmute("C09"=prev_f)) %>%
  cbind(tab_prev_21.B %>% transmute("C21"=prev_f)) %>% 
  `rownames<-`(NULL) %>% select(Year, Prior, C09, C21) -> tab1B

#Join and save
char1 <- "Decreased eGFR"
char2 <- "Prior diagnosis or decreased eGFR"

rbind(c(char1,"","",""), tab1A, c(char2,"","",""), tab1B) %>%  `names<-`(c(
  "ENSANUT year","Prior diagnosis","CKD-EPI 2009","CKD-EPI 2021")) %>% 
  flextable() %>% add_header_row(top=T, values=c(
    "ENSANUT year",rep("Prevalence (95% CI)",3))) %>% 
  merge_h(part = "header") %>% merge_v(part = "header") %>% 
  merge_at(i = 1, j = 1:4, part = "body") %>% 
  merge_at(i = 7, j = 1:4, part = "body") %>% bold(part="header") %>% 
  bold(i = c(1,7), j = 1:4, part = "body") %>%
  italic(i = c(1,7), j = 1:4, part = "body") %>% 
  align(align = "center", part = "all") %>% autofit() %>% save_as_docx(
    path="Table1.docx", pr_section=prop_section(
      page_size = page_size(orient = "landscape")))


#### Supp Table 1: Population characteristics--------------- ####
ensanut_all <- rbind(
  ensanut_2016 %>% transmute(
    "Year"=2016, "Age"=Age, "Sex"=Sex, "Smo"=Smoking, "Ind"=Indigenous,
    "Urb"=Area, "BMI"=BMI, "WC"=Waist, "Glu"=Glucose, "A1c"=HBA1C, "TC"=CT, TG,
    Creatinine, "HTN"=hypertension, "T2D"=diabetes, "CVD"=HX_CVD, "CKD"=HX_CKD),
  ensanut_2018 %>% transmute(
    "Year"=2018, "Age"=Age, "Sex"=Sex, "Smo"=Smoking, "Ind"=Indigenous,
    "Urb"=Area, "BMI"=BMI, "WC"=Waist, "Glu"=Glucose, "A1c"=HBA1C, "TC"=CT, TG,
    Creatinine, "HTN"=hypertension, "T2D"=diabetes, "CVD"=HX_CVD, "CKD"=HX_CKD),
  ensanut_2020 %>% transmute(
    "Year"=2020, "Age"=Age, "Sex"=Sex, "Smo"=Smoking, "Ind"=Indigenous,
    "Urb"=Area, "BMI"=BMI, "WC"=Waist, "Glu"=Glucose, "A1c"=HBA1C, "TC"=CT, TG,
    Creatinine, "HTN"=hypertension, "T2D"=diabetes, "CVD"=HX_CVD, "CKD"=HX_CKD),
  ensanut_2021 %>% transmute(
    "Year"=2021, "Age"=Age, "Sex"=Sex, "Smo"=Smoking, "Ind"=Indigenous,
    "Urb"=Area, "BMI"=BMI, "WC"=Waist, "Glu"=Glucose, "A1c"=HBA1C, "TC"=CT, TG,
    Creatinine, "HTN"=hypertension, "T2D"=diabetes, "CVD"=HX_CVD, "CKD"=HX_CKD),
  ensanut_2022 %>% transmute(
    "Year"=2022, "Age"=Age, "Sex"=Sex, "Smo"=Smoking, "Ind"=Indigenous,
    "Urb"=Area, "BMI"=BMI, "WC"=Waist, "Glu"=Glucose, "A1c"=HBA1C, "TC"=CT, TG,
    Creatinine, "HTN"=hypertension, "T2D"=diabetes, "CVD"=HX_CVD, "CKD"=HX_CKD),
  ensanut_2023 %>% transmute(
    "Year"=2023, "Age"=Age, "Sex"=Sex, "Smo"=Smoking, "Ind"=Indigenous,
    "Urb"=Area, "BMI"=BMI, "WC"=Waist, "Glu"=Glucose, "A1c"=HBA1C, "TC"=CT, TG,
    Creatinine, "HTN"=hypertension, "T2D"=diabetes, "CVD"=HX_CVD, "CKD"=HX_CKD))

lab.in <- ensanut_all[-1] %>% names; lab.out <- c(
  "Age (years)", "Sex (women)", "Smoking status", "Indigenous language",
  "Urban area", "Body mass index (kg/m^2)", "Waist circumference (cm)",
  "Fasting plasma glucose (mg/dL)", "HbA1c (%)", "Total cholesterol (mg/dL)",
  "Triglycerides (mg/dL)", "Creatinine (mg/L)", "Hypertension (%)",
  "Type 2 diabetes (%)","Cardiovascular disease (%)", "Prior CKD (%)"
  ); lab.all <- list(
    lab.in[1]~lab.out[1], lab.in[2]~lab.out[2], lab.in[3]~lab.out[3],
    lab.in[4]~lab.out[4], lab.in[5]~lab.out[5], lab.in[6]~lab.out[6],
    lab.in[7]~lab.out[7], lab.in[8]~lab.out[8], lab.in[9]~lab.out[9],
    lab.in[10]~lab.out[10], lab.in[11]~lab.out[11], lab.in[12]~lab.out[12],
    lab.in[13]~lab.out[13], lab.in[14]~lab.out[14], lab.in[15]~lab.out[15],
    lab.in[16]~lab.out[16])

gtsummary::tbl_summary(
  ensanut_all %>% filter(!is.na(Creatinine)), by=Year, missing_text="Missing",
  label = lab.all) %>% bold_labels() %>% modify_table_body(~.x %>% mutate(
    stat_3 = ifelse(stat_3%in%c(
      "0 (NA%)", "NA (NA, NA)","0 (0%)", "1 (100%)"), "-", stat_3))) %>% 
  as_flex_table() %>% align(align = "center", part = "all") %>% autofit() %>%
  save_as_docx(
    path="SuppTable1.docx", pr_section=prop_section(
      page_size = page_size(orient = "landscape")))

#### Supp Table 2: Prevalence of CKD by subgroup------------ ####
tab_prev_12[-(1:5),] %>% transmute(Year=factor(Year), group, "C09"=prev_f) %>%
  cbind(tab_prev_21[-(1:5),] %>% transmute("C21"=prev_f)) %>% mutate(
    "strata"=c(rep(1,12),rep(2,12), rep(3,12), rep(4,12)) %>% ordered(labels=c(
      "Sex","Age group","Diabetes\nstatus","Hypertension\nstatus"))) %>%
  select(strata, group, Year, C09, C21) %>% arrange(strata, group) %>%
  `rownames<-`(NULL) %>% flextable() %>% set_header_labels(values = c(
    "Subgroup","Subgroup","ENSANUT year","CKD-EPI 2009","CKD-EPI 2021")) %>%
  add_header_row(top=T, values=c("Subgroup","Subgroup", "ENSANUT year",
                                 rep("Prevalence (95% CI)",2))) %>% 
  merge_h(part = "header") %>% merge_v(part = "header") %>% 
  merge_v(part = "body") %>% bold(part="header") %>% bold(j=1, part="body") %>% 
  hline(i=c(5,10,15,20,25,30,35), border=fp_border(
    color="black", width=0.5), part="body") %>% 
  align(align = "center", part = "all") %>% autofit() %>% save_as_docx(
    path="SuppTable2.docx", pr_section=prop_section(
      page_size = page_size(orient = "landscape")))

#### Drafts ####
#### BMI
val.min <- c(CKD12_prev_bmi$lIC95, CKD21_prev_bmi$lIC95) %>% min %>% floor
val.max <- c(CKD12_prev_bmi$uIC95, CKD21_prev_bmi$uIC95) %>% max %>% ceiling
CKD12_prev_bmi %>% mutate("A"=c("CKD-EPI 2012")) %>%
  mutate(group=ordered(group, c("Non- obesity","Obesity"), c("No obesity","Obesity"))) %>% 
  ggplot(aes(x=as.numeric(Year), y=prop,group=factor(group),colour=factor(group))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("CKD-EPI 2009") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E49", "#392220")) + ylim(val.min,val.max) +
  theme(legend.position = "bottom", plot.title = element_text(face = "italic", hjust = 0.5, size = 14)) -> F2a
CKD21_prev_bmi %>% mutate("A"=c("CKD-EPI 2021")) %>%
  mutate(group=ordered(group, c("Non- obesity","Obesity"), c("No obesity","Obesity"))) %>% 
  ggplot(aes(x=as.numeric(Year), y=prop,group=factor(group),colour=factor(group))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("CKD-EPI 2021") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E49", "#392220")) + ylim(val.min,val.max) +
  theme(legend.position = "bottom", plot.title = element_text(face = "italic", hjust = 0.5, size = 14)) -> F2b
F2.bmi <- ggarrange(F2a, F2b, nrow=1, ncol=2) %>% annotate_figure(top = text_grob(
  "CKD prevalence by BMI category", color = "black", face = "bold", size = 18))
#### Ind
val.min <- c(CKD12_prev_ind$lIC95, CKD21_prev_ind$lIC95) %>% min %>% floor
val.max <- c(CKD12_prev_ind$uIC95, CKD21_prev_ind$uIC95) %>% max %>% ceiling
CKD12_prev_ind %>% mutate("A"=c("CKD-EPI 2012")) %>%
  mutate(group=ordered(group, c("Non-indigenous","Indigenous"))) %>% 
  ggplot(aes(x=as.numeric(Year), y=prop,group=factor(group),colour=factor(group))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("CKD-EPI 2009") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E49", "#392220")) + ylim(val.min,val.max) +
  theme(legend.position = "bottom", plot.title = element_text(face = "italic", hjust = 0.5, size = 14)) -> F2a
CKD21_prev_ind %>% mutate("A"=c("CKD-EPI 2021")) %>%
  mutate(group=ordered(group, c("Non-indigenous","Indigenous"))) %>% 
  ggplot(aes(x=as.numeric(Year), y=prop,group=factor(group),colour=factor(group))) + geom_line(size=1.5) +
  geom_point(size=2) + scale_x_continuous(breaks = c(2016, 2018, 2020, 2021, 2022)) +
  geom_errorbar(aes(ymin=prop+IC95, ymax=prop-IC95), width=.2, position=position_dodge(0.01))+
  ylab("Weighted prevalence (%)") + xlab("ENSANUT cycle")+ labs(colour=NULL) + ggtitle("CKD-EPI 2021") +
  theme_pubclean() + scale_color_manual(values=c("#CA4E49", "#392220")) + ylim(val.min,val.max) +
  theme(legend.position = "bottom", plot.title = element_text(face = "italic", hjust = 0.5, size = 14)) -> F2b
F2.ind <- ggarrange(F2a, F2b, nrow=1, ncol=2) %>% annotate_figure(top = text_grob(
  "CKD prevalence by Indigenous identity", color = "black", face = "bold", size = 18))




ensanut_2023$bp_diag[ensanut_2023$HX_HBP == 0 & (ensanut_2023$SBP < 130 & ensanut_2023$DBP < 80)] <- 0 #Sin hipertensión
ensanut_2023$bp_diag[ensanut_2023$HX_HBP == 1] <- 1 #Diagnosticado
ensanut_2023$bp_diag[ensanut_2023$HX_HBP == 0 & (ensanut_2023$SBP >= 130 | ensanut_2023$DBP >= 80)] <- 2 #No diagnosticado
ensanut_2023$bp_diag <- factor(ensanut_2023$bp_diag, labels = c("Without hypertension", "Diagnosed hypertension", "Undiagnosed hypertension"))

ensanut_2022$bp_diag[ensanut_2022$HX_HBP == 0 & (ensanut_2022$SBP < 130 & ensanut_2022$DBP < 80)] <- 0 #Sin hipertensión
ensanut_2022$bp_diag[ensanut_2022$HX_HBP == 1] <- 1 #Diagnosticado
ensanut_2022$bp_diag[ensanut_2022$HX_HBP == 0 & (ensanut_2022$SBP >= 130 | ensanut_2022$DBP >= 80)] <- 2 #No diagnosticado
ensanut_2022$bp_diag <- factor(ensanut_2022$bp_diag, labels = c("Without hypertension", "Diagnosed hypertension", "Undiagnosed hypertension"))

svy_2022 <- svyby(~bp_diag, ~diabetes, ensanut_2022_survey, svymean, na.rm = TRUE)
svy_2023 <- svyby(~bp_diag, ~diabetes, ensanut_2023_survey, svymean, na.rm = TRUE)


svy_total_2022 <- svymean(~bp_diag, ensanut_2022_survey, na.rm = TRUE)
svy_total_2023 <- svymean(~bp_diag, ensanut_2023_survey, na.rm = TRUE)

v <- svy_2022 %>% select(1:4) %>% 
  pivot_longer(cols = c(2:4),
               names_to = "bp_2022",
               values_to = "prop")

vv <- svy_2023 %>% select(1:4) %>% 
  pivot_longer(cols = c(2:4),
               names_to = "bp_2023",
               values_to = "prop") %>% 
  select(-1)

a <- cbind(v, vv)
