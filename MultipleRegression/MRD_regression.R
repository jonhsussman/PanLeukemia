#### Introduction ####
# This script analyzes correlation between leukemic cell type fraction and MRD.
# Multiple regression is used to analyze correlation while adjusting for
# biological covariates. Figures 3E and S6D

#### Load Packages ####
library(tidyverse)
library(readxl)
library(Matrix)
library(Hmisc)
library(ggplot2)
library(scales)
library(patchwork)

#### Load in and Format Data ####

## AML (TARGET) Data
NCI.cell.frac <- read_delim("./scRNA/Scripts/CibersortX/AML/Output/CIBERSORTx_Job1_Results_AML_NCI_min0.5_rep50_sampling0.5.csv",
                            delim=",")
FHCRC.cell.frac <- read_delim("./scRNA/Scripts/CibersortX/AML/Output/CIBERSORTx_Job1_Results_AML_FHCRC_min0.5_rep50_sampling0.5.csv",
                              delim=",")[-2,] #repeated entry
AML.cell.frac <- rbind(NCI.cell.frac, FHCRC.cell.frac)
# Generate patient ID column from sample names
AML.cell.frac$Patient_ID <- sapply(AML.cell.frac$Mixture, function(n) {paste(strsplit(n,"-",fixed=TRUE)[[1]][c(1:3)],collapse="-")})
AML.cell.frac <- AML.cell.frac[,c(1,dim(AML.cell.frac)[2],2:(dim(AML.cell.frac)[2]-1))]
# Keep only IDX samples
AML.cell.frac <- AML.cell.frac[sort(c(grep("09A",AML.cell.frac$Mixture),grep("03A",AML.cell.frac$Mixture),
                                      grep("09B",AML.cell.frac$Mixture),grep("03B",AML.cell.frac$Mixture))),]
# Clinical Data
AML.dir <- "./BulkRNA/Target/AML/"
clinical.disc <- read_xlsx(paste0(AML.dir,"TARGET_AML_ClinicalData_Discovery_20211201.xlsx"))
clinical.valid <- read_xlsx(paste0(AML.dir,"TARGET_AML_ClinicalData_Validation_20211201.xlsx"))
clinical.other <- read_xlsx(paste0(AML.dir,"TARGET_AML_ClinicalData_AML1031_20211201.xlsx"))
clinical.disc["source"] <- rep("disc", dim(clinical.disc)[1])
clinical.valid["source"] <- rep("valid", dim(clinical.valid)[1])
clinical.other["source"] <- rep("other", dim(clinical.other)[1])
# Combine into one data frame. NOTE: This uses discovery data as default source
# for repeated clinical data (listed first), alternate therapy may be missing
AML.outcome <- rbind(clinical.disc,clinical.valid,clinical.other) %>%
  distinct(`TARGET USI`, .keep_all=TRUE)
colnames(AML.outcome)[1] <- "Patient_ID"
# Set MRD names from clinical data for later
AML.MRD1.name <- "MRD % at end of course 1"
AML.MRD2.name <- "MRD % at end of course 2"


## B-ALL (TARGET) Data
# Deconvoluted leukemic cell fractions
BALL.frac.file <- "./scRNA/Scripts/CibersortX/BALL/Output/CIBERSORTx_Job1_Results_BALL_withMyeloid_BCCA_min0.5_rep50_sampling0.5.csv"
BALL.cell.frac <- read_delim(BALL.frac.file, delim=",")
# Generate patient ID column from sample names
BALL.cell.frac$Patient_ID <- sapply(BALL.cell.frac$Mixture, function(n) {paste(strsplit(n,"-",fixed=TRUE)[[1]][c(1:3)],collapse="-")})
BALL.cell.frac <- BALL.cell.frac[,c(1,dim(BALL.cell.frac)[2],2:(dim(BALL.cell.frac)[2]-1))]
# Check sample matrices for IDX samples
BALL.dir <- "./BulkRNA/Target/Phase2/"
sample.disc <- read_xlsx(paste0(BALL.dir,"TARGET_ALL_SampleMatrix_Phase2_Discovery_20190606.xlsx"))
sample.valid <- read_xlsx(paste0(BALL.dir,"TARGET_ALL_SampleMatrix_Phase2_Validation_20190606.xlsx"))
IDX.samples <- c(na.omit(sample.disc$`Diagnostic Tumor RNA Sample ID`),na.omit(sample.valid$`Diagnostic Tumor RNA Sample ID`))
IDX.samples <- unlist(strsplit(IDX.samples,","))
# Keep only IDX samples
BALL.cell.frac <- BALL.cell.frac[BALL.cell.frac$Mixture %in% IDX.samples,]
# Clinical Data
clinical.disc <- read_xlsx(paste0(BALL.dir,"TARGET_ALL_ClinicalData_Phase_II_Discovery_20211118.xlsx"))
clinical.valid <- read_xlsx(paste0(BALL.dir,"TARGET_ALL_ClinicalData_Phase_II_Validation_20211118.xlsx"))
clinical.all <- read_xlsx(paste0(BALL.dir,"TARGET_ALL_ClinicalData_Dicentric_20211118.xlsx"))
clinical.disc["source"] <- rep("disc", dim(clinical.disc)[1])
clinical.valid["source"] <- rep("valid", dim(clinical.valid)[1])
clinical.all["source"] <- rep("all", dim(clinical.all)[1])
# Combine into one data frame. NOTE: This uses discovery data as default source
# for repeated clinical data (listed first), alternate therapy may be missing
BALL.outcome <- rbind(clinical.disc,clinical.valid,clinical.all) %>%
  distinct(`TARGET USI`, .keep_all=TRUE)
colnames(BALL.outcome)[1] <- "Patient_ID"
# Set MRD names from clinical data for later
BALL.MRD1.name <- "MRD Day 29"
BALL.MRD2.name <- "MRD End Consolidation"


## T-ALL (AALL0434) Data
# Deconvoluted leukemic cell fractions
TALL.b1.cell.frac <- read_delim("./scRNA/Scripts/CibersortX/All_TALL/Output/CIBERSORTx_Job1_Results_AllTALL_X01_batch1_detailPopulations_min0.5_rep50_sampling0.5.csv",
                                delim=",")
TALL.b2.cell.frac <- read_delim("./scRNA/Scripts/CibersortX/All_TALL/Output/CIBERSORTx_Job1_Results_AllTALL_X01_batch2_detailPopulations_min0.5_rep50_sampling0.5.csv",
                                delim=",")
TALL.cell.frac <- rbind(TALL.b1.cell.frac,TALL.b2.cell.frac)
colnames(TALL.cell.frac)[1] <- "Patient_ID"
# Clinical Data
TALL.outcome <- read_delim("./AllTALL_X01Cohort_ClinicalInformation.txt", delim="\t")
colnames(TALL.outcome)[1] <- "Patient_ID"
# Set MRD names from clinical data for later
TALL.MRD1.name <- "D29.MRD"
TALL.MRD2.name <- "EOC.MRD"


#### Build Data Frames with Fractions, Covariates and Outcomes ####

# Combine/Rename HSPC-like fractions and combine lineage-like fractions
# AML
hspc.like <- c("Progenitor-Multipotent-like")
lineage.like <- c("Progenitor-Myeloid-like", "Mature-Myeloid-like")
AML.cell.frac <- mutate(AML.cell.frac,
                          'HSPC-like'=rowSums(AML.cell.frac[,hspc.like]),
                          'Lineage-like'=rowSums(AML.cell.frac[,lineage.like]))
# B-ALL
hspc.like <- c("Progenitor-Multipotent-like")
lineage.like <- c("Progenitor-Lymphoid-like", "Progenitor-Myeloid-like", "B-like")
BALL.cell.frac <- mutate(BALL.cell.frac,
                           'HSPC-like'=rowSums(BALL.cell.frac[,hspc.like]),
                           'Lineage-like'=rowSums(BALL.cell.frac[,lineage.like]))
# T-ALL
hspc.like <- c("LMPP-like", "HSPC-like")
lineage.like <- c("Pro-T-like", "Pre-T-like", "CLP-like", "GMP-like", "DP-like", "ETP-like")
TALL.cell.frac <- mutate(TALL.cell.frac,
                           'HSPC-like'=rowSums(TALL.cell.frac[,hspc.like]),
                           'Lineage-like'=rowSums(TALL.cell.frac[,lineage.like]))

## Build Data Frame
# AML
AML.frac.out.df <- inner_join(AML.cell.frac, AML.outcome, by="Patient_ID") %>% 
  mutate_at(c(AML.MRD1.name, AML.MRD2.name), as.numeric)
# Convert select columns to factors and calculate normalized MRD
AML.frac.out.df <- rename(AML.frac.out.df, MRD1=sym(AML.MRD1.name),MRD2=sym(AML.MRD2.name)) %>%
  mutate(MRD1.norm=if_else(MRD1<0.01,log10(0.005),log10(MRD1)), 
         MRD2.norm=if_else(MRD2<0.01,log10(0.005),log10(MRD2)),
         Protocol=as.factor(Protocol))
# BALL
BALL.frac.out.df <- inner_join(BALL.cell.frac, BALL.outcome, by="Patient_ID") %>% 
  mutate_at(c(BALL.MRD1.name, BALL.MRD2.name), as.numeric)
# Convert select columns to factors and calculate normalized MRD
BALL.frac.out.df <- rename(BALL.frac.out.df, MRD1=sym(BALL.MRD1.name),MRD2=sym(BALL.MRD2.name)) %>%
  mutate(MRD1.norm=if_else(MRD1<0.01,log10(0.005),log10(MRD1)), 
         MRD2.norm=if_else(MRD2<0.01,log10(0.005),log10(MRD2)),
         `ALL Molecular Subtype`=as.factor(`ALL Molecular Subtype`),
         Protocol=as.factor(Protocol))
# TALL
TALL.frac.out.df <- inner_join(TALL.cell.frac, TALL.outcome, by="Patient_ID") %>% 
  mutate_at(c(TALL.MRD1.name, TALL.MRD2.name), as.numeric)
# Convert select columns to factors and calculate normalized MRD
TALL.frac.out.df <- rename(TALL.frac.out.df, MRD1=sym(TALL.MRD1.name),MRD2=sym(TALL.MRD2.name)) %>%
  mutate(MRD1.norm=if_else(MRD1<0.01,log10(0.005),log10(MRD1)), 
         MRD2.norm=if_else(MRD2<0.01,log10(0.005),log10(MRD2)),
         in.target=as.factor(in.target), risk.group=as.factor(risk.group), ETP=as.factor(ETP))


#### MULTIPLE REGRESSION ####

## Subset and filter data for linear regression (remove NAs)
cell.types=c("HSPC-like","Lineage-like")
# AML
AML.covariates=c("Age at Diagnosis in Days", "WBC at Diagnosis", "Protocol")
AML.frac.out.filt <- na.omit(AML.frac.out.df[,c("MRD1.norm", cell.types, AML.covariates)])
# BALL
BALL.covariates=c("Age at Diagnosis in Days", "WBC at Diagnosis", "Protocol")
BALL.frac.out.filt <- na.omit(BALL.frac.out.df[,c("MRD1.norm", cell.types, BALL.covariates)])
# TALL
TALL.covariates=c("dx.age.days", "dx.wbc", "in.target", "ETP")
TALL.frac.out.filt <- na.omit(TALL.frac.out.df[,c("MRD1.norm", cell.types, TALL.covariates)])

## Scale cell-fraction data
# AML
AML.frac.out.filt[,grep("like",colnames(AML.frac.out.filt))] <- 
  scale(AML.frac.out.filt[,grep("like",colnames(AML.frac.out.filt))])
# BALL
BALL.frac.out.filt[,grep("like",colnames(BALL.frac.out.filt))] <-
  scale(BALL.frac.out.filt[,grep("like",colnames(BALL.frac.out.filt))])
# TALL
TALL.frac.out.filt[,grep("like",colnames(TALL.frac.out.filt))] <-
  scale(TALL.frac.out.filt[,grep("like",colnames(TALL.frac.out.filt))])

## Run Regression and Plot (Fig. 3E)
# AML
lm.res <- lm(MRD1.norm ~  ., data=AML.frac.out.filt)
summary(lm.res)
## Create data frame
end.index <- length(cell.types)+length(AML.covariates)+1
results.df <- tibble(variable = gsub("`","",names(lm.res$coefficients[2:end.index])),
                     coeff = lm.res$coefficients[2:end.index],
                     ci.low = confint(lm.res,level=0.95)[2:end.index,1],
                     ci.high = confint(lm.res,level=0.95)[2:end.index,2],
                     p.value = summary(lm.res)$coefficients[2:end.index,4]) %>%
  dplyr::arrange(coeff)
y.order <- c("Lineage-like", "HSPC-like")
#pdf(file="./Figures/MultipleRegression_MRD_summed_cellFraction_adj_Age_WBC_Protocol_AML.pdf", 4, 2)
ggplot(results.df[results.df$variable %in% y.order,], aes(y=factor(variable,levels=y.order))) +
  geom_point(aes(x=coeff)) +
  geom_errorbar(aes(xmin=ci.low,xmax=ci.high), width=0.2) +
  geom_vline(xintercept=0, linetype="dotted") +
  ylab("") + xlab("Coefficient of Linear Regression") +
  ggtitle("AML (TARGET)")
#dev.off()

# BALL
lm.res <- lm(MRD1.norm ~  ., data=BALL.frac.out.filt)
summary(lm.res)
## Create data frame
end.index <- length(cell.types)+length(BALL.covariates)+1
results.df <- tibble(variable = gsub("`","",names(lm.res$coefficients[2:end.index])),
                     coeff = lm.res$coefficients[2:end.index],
                     ci.low = confint(lm.res,level=0.95)[2:end.index,1],
                     ci.high = confint(lm.res,level=0.95)[2:end.index,2],
                     p.value = summary(lm.res)$coefficients[2:end.index,4]) %>%
  dplyr::arrange(coeff)
y.order <- c("Lineage-like", "HSPC-like")
#pdf(file="./Figures/MultipleRegression_MRD_summed_cellFraction_adj_Age_WBC_Protocol_BALL.pdf", 4, 2)
ggplot(results.df[results.df$variable %in% y.order,], aes(y=factor(variable,levels=y.order))) +
  #ggplot(results.df, aes(y=variable)) +
  geom_point(aes(x=coeff)) +
  geom_errorbar(aes(xmin=ci.low,xmax=ci.high), width=0.2) +
  geom_vline(xintercept=0, linetype="dotted") +
  ylab("") + xlab("Coefficient of Linear Regression") +
  ggtitle("B-ALL (TARGET)")
#dev.off()


# TALL
lm.res <- lm(MRD1.norm ~  ., data=TALL.frac.out.filt)
#summary(lm.res)
## Create data frame
end.index <- length(cell.types)+length(TALL.covariates)+1
results.df <- tibble(variable = gsub("`","",names(lm.res$coefficients[2:end.index])),
                     coeff = lm.res$coefficients[2:end.index],
                     ci.low = confint(lm.res,level=0.95)[2:end.index,1],
                     ci.high = confint(lm.res,level=0.95)[2:end.index,2],
                     p.value = summary(lm.res)$coefficients[2:end.index,4]) %>%
  dplyr::arrange(coeff)
y.order <- c("Lineage-like", "HSPC-like")
#pdf(file="./Figures/MultipleRegression_MRD_summed_cellFraction_adj_Age_WBC_Protocol_ETP_TALL.pdf", 4, 2)
ggplot(results.df[results.df$variable %in% y.order,], aes(y=factor(variable,levels=y.order))) +
  geom_point(aes(x=coeff)) +
  geom_errorbar(aes(xmin=ci.low,xmax=ci.high), width=0.2) +
  geom_vline(xintercept=0, linetype="dotted") +
  ylab("") + xlab("Coefficient of Linear Regression") +
  ggtitle("T-ALL (AALL0434)")
#dev.off()


#### Repeat with Total XXI B-ALL Data Separated by Subtype ####

## Load Total XXI Data
# Cell fraction from deconvoluted bulk RNA
BALL.cell.frac <- read_delim("./CIBERSORTx_Job8_Results_BALL_withMyeloid_NMed_Pharmocotyping_min0.5_rep50_sampling0.5.csv", delim=",")
colnames(BALL.cell.frac)[1] <- "Patient_ID"
# Annotated outcome and LC50 values
outcome <- read_delim("./Pharmacotype_LC50_full.txt", delim="\t", skip=1)
colnames(outcome)[1] <- "Patient_ID"
MRD1.name = "Day 15 MRD (%)"
MRD2.name = "Day 42 or 46 MRD (%)"

## Build Data Frames with Fractions, Covariates and Outcomes
# Combine/Rename HSPC-like fractions and combine lineage-like fractions
hspc.like <- c("Progenitor-Multipotent-like")
lineage.like <- c("Progenitor-Lymphoid-like", "Progenitor-Myeloid-like", "B-like")
BALL.cell.frac <- mutate(BALL.cell.frac,
                         'HSPC-like'=rowSums(BALL.cell.frac[,hspc.like]),
                         'Lineage-like'=rowSums(BALL.cell.frac[,lineage.like]))
# Build Data Frame
BALL.frac.out.df <- inner_join(BALL.cell.frac, outcome, by="Patient_ID") %>%
  mutate_at(c(MRD1.name, MRD2.name), as.numeric)
BALL.frac.out.df <- BALL.frac.out.df[BALL.frac.out.df$Immunophenotype=="B",]
# Convert select columns to factors and calculate normalized MRD
BALL.frac.out.df <- rename(BALL.frac.out.df, MRD1=sym(MRD1.name),MRD2=sym(MRD2.name)) %>%
  mutate(MRD1.norm=if_else(MRD1<0.01,log10(0.005),log10(MRD1)), 
         MRD2.norm=if_else(MRD2<0.01,log10(0.005),log10(MRD2)),
         `Molecular subtype`=as.factor(`Molecular subtype`),
         Protocol=as.factor(Protocol), `NCI risk`=as.factor(`NCI risk`))


# Run Regression for each molecular subtype
cell.types=c("HSPC-like", "Lineage-like")
covariates = c("Age at diagnosis (years)", "WBC at diagnosis (x 109/L)", "Protocol")
# Make results df
results.df <- tibble(variable = c(), molecular.subtype = c(), N=c(),
                          coeff = c(), ci.low = c(), ci.high = c(), p.value = c())
# Separate df for results that have too low N, leads to NAs in coefficients and pvalues
lowN.results.df <- tibble(variable = c(), molecular.subtype = c(), N=c(),
                          coeff = c())
for (msub in unique(BALL.frac.out.df$`Molecular subtype`)) {
#for (msub in unique(frac.out.df$`ALL Molecular Subtype`)) {
  frac.out.filt <- na.omit(BALL.frac.out.df[BALL.frac.out.df$`Molecular subtype`==msub,
                                            c("MRD1.norm", cell.types, covariates)])
  # Get num of samples in this subtype
  m.n <- dim(frac.out.filt)[1]
  # In case where all patients had same treatment protocol, remove from covariates
  if (length(unique(frac.out.filt$Protocol))==1) frac.out.filt <- subset(frac.out.filt, select=-Protocol)
  if (m.n>1) {
    # Scale cell-fraction data
    frac.out.filt[,cell.types] <- scale(frac.out.filt[,cell.types])
    # Run Regression
    lm.res <- lm(MRD1.norm ~  ., data=frac.out.filt)
    # Create data frame
    end.index <- length(cell.types)+1
    # Check if n was sufficient, all coeffients and p values estimated
    if (dim(summary(lm.res)$coefficients)[1]>2 & !is.nan(summary(lm.res)$coefficients[2,4])) {
      mresults.df <- tibble(variable = gsub("`","",names(lm.res$coefficients[2:end.index])),
                            molecular.subtype = rep(paste0(msub," (N=",m.n,")"), end.index-1),
                            N = rep(m.n, end.index-1),
                            coeff = lm.res$coefficients[2:end.index],
                            ci.low = confint(lm.res,level=0.95)[2:end.index,1],
                            ci.high = confint(lm.res,level=0.95)[2:end.index,2],
                            p.value = summary(lm.res)$coefficients[2:end.index,4]) %>%
        dplyr::arrange(coeff)
      results.df <- rbind(results.df,mresults.df)
    }
    else {
      mresults.df <- tibble(variable = gsub("`","",names(lm.res$coefficients[2:end.index])),
                            molecular.subtype = rep(paste0(msub," (N=",m.n,")"), end.index-1),
                            N = rep(m.n, end.index-1),
                            coeff = lm.res$coefficients[2:end.index]) %>%
        dplyr::arrange(coeff)
      lowN.results.df <- rbind(lowN.results.df,mresults.df)
    }
  }
}

# Plot results,  excluding low N subtypes (Fig. S6D)
y.order <- c("Lineage-like", "HSPC-like")
#pdf(file="./Figures/MultipleRegression_MRD_summed_cell-fraction_adj_Age_WBC_Protocol_BALL_BCCA_byMolSubtype.pdf", 6, 4)
ggplot(results.df, aes(y=factor(variable,levels=y.order))) +
  geom_point(aes(x=coeff)) +
  geom_errorbar(aes(xmin=ci.low,xmax=ci.high), width=0.2) +
  geom_vline(xintercept=0, linetype="dotted") +
  ylab("Cell Type") +
  xlab("Coefficient of Linear Regression") +
  facet_wrap(~molecular.subtype, ncol=3, scales = "free_x",
             labeller=label_wrap_gen(width=20)) +
  ggtitle("Total XXI Day 15 MRD Correlation (B-ALL)")
#dev.off()
