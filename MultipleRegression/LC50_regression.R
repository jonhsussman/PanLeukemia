## Introduction ----
# This script analyzes correlation between leukemic cell type fraction and LC50
# values for for a panel of 18 conventional chemo and targeted drugs.
# Multiple regression is used to analyze correlation while adjusting for
# biological covariates. Figures 3F, 6G

## Load Packages ----
library(tidyverse)
library(readxl)
library(Matrix)
library(Hmisc)
library(ggplot2)
library(scales)
library(patchwork)

### Load in and Format Data ####

# Cell fraction from deconvoluted bulk RNA
BALL.cell.frac <- read_delim("./CIBERSORTx_Job8_Results_BALL_withMyeloid_NMed_Pharmocotyping_min0.5_rep50_sampling0.5.csv", delim=",")
colnames(BALL.cell.frac)[1] <- "Patient_ID"
TALL.cell.frac <- read_delim("./CibersortX_TALL_fraction.csv", delim=",")
colnames(TALL.cell.frac)[1] <- "Patient_ID"
# Annotated outcome and LC50 values
outcome <- read_delim("./Pharmacotype_LC50_full.txt", delim="\t", skip=1)
colnames(outcome)[1] <- "Patient_ID"

#### Build Data Frames with Fractions, Covariates and Outcomes ####

# Combine/Rename HSPC-like fractions and combine lineage-like fractions
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
# Build Data Frame
BALL.frac.out.df <- inner_join(BALL.cell.frac, outcome, by="Patient_ID")
BALL.frac.out.df <- BALL.frac.out.df[BALL.frac.out.df$Immunophenotype=="B",]
TALL.frac.out.df <- inner_join(TALL.cell.frac, outcome, by="Patient_ID")
TALL.frac.out.df <- TALL.frac.out.df[TALL.frac.out.df$Immunophenotype=="T",]
# Convert select columns to factors
BALL.frac.out.df <- mutate(BALL.frac.out.df,
                           `Molecular subtype`=as.factor(`Molecular subtype`),
                           Protocol=as.factor(Protocol), `NCI risk`=as.factor(`NCI risk`))
TALL.frac.out.df <- mutate(TALL.frac.out.df,
                           `Molecular subtype`=as.factor(`Molecular subtype`),
                           Protocol=as.factor(Protocol), `NCI risk`=as.factor(`NCI risk`))

#### MULTIPLE REGRESSION ####

## Set cell types and covariates
cell.types = c("HSPC-like", "Lineage-like")
covariates = c("Age at diagnosis (years)", "WBC at diagnosis (x 109/L)")

## B-ALL Regression
# Create full results df
B.LC50.res.df <- tibble(variable=c(),drug=c(),coeff=c(),ci.low=c(),ci.high=c(),p.value=c())
# Regress each drug LC50 against fractions and covariates
drug.names <- colnames(BALL.frac.out.df)[grep("normalized",colnames(BALL.frac.out.df))]
for (drug in drug.names) {
  # Subset for necessary variables and filter out na's
  frac.out.filt <- na.omit(BALL.frac.out.df[,c(drug,cell.types,covariates)])
  # Scale cell-fraction data
  frac.out.filt[,cell.types] <- scale(frac.out.filt[,cell.types])
  # Scale all numerical variables?
  frac.out.filt[,c("Age at diagnosis (years)", "WBC at diagnosis (x 109/L)")] <- 
    scale(frac.out.filt[,c("Age at diagnosis (years)", "WBC at diagnosis (x 109/L)")])
  # Define formula and run regression
  lm.formula <- as.formula(paste(drug,"~ ."))
  reg.res <- lm(lm.formula, data=frac.out.filt)
  # Summarize results (include covariates)
  reg.res.df <- tibble(gsub("`","",names(reg.res$coefficients[-1])), 
                       rep(strsplit(drug,split="_")[[1]][1],length(reg.res$coefficients[-1])),
                       reg.res$coefficients[-1], 
                       confint(reg.res,level=0.95)[-1,1],
                       confint(reg.res,level=0.95)[-1,2],
                       summary(reg.res)$coefficients[-1,4])
  colnames(reg.res.df) <- c("variable", "drug", "coeff",
                            "ci.low", "ci.high","p.value")
  B.LC50.res.df <- rbind(B.LC50.res.df,reg.res.df)
}

## T-ALL Regression
# Create full results df
T.LC50.res.df <- tibble(variable=c(),drug=c(),coeff=c(),ci.low=c(),ci.high=c(),p.value=c())
# Regress each drug LC50 against fractions and covariates
drug.names <- colnames(TALL.frac.out.df)[grep("normalized",colnames(TALL.frac.out.df))]
for (drug in drug.names) {
  # Subset for necessary variables and filter out na's
  frac.out.filt <- na.omit(TALL.frac.out.df[,c(drug,cell.types,covariates)])
  # Scale cell-fraction data
  frac.out.filt[,cell.types] <- scale(frac.out.filt[,cell.types])
  # Scale all numerical variables?
  frac.out.filt[,c("Age at diagnosis (years)", "WBC at diagnosis (x 109/L)")] <- 
    scale(frac.out.filt[,c("Age at diagnosis (years)", "WBC at diagnosis (x 109/L)")])
  # Define formula and run regression
  lm.formula <- as.formula(paste(drug,"~ ."))
  reg.res <- lm(lm.formula, data=frac.out.filt)
  # Summarize results (include covariates)
  reg.res.df <- tibble(gsub("`","",names(reg.res$coefficients[-1])), 
                       rep(strsplit(drug,split="_")[[1]][1],length(reg.res$coefficients[-1])),
                       reg.res$coefficients[-1], 
                       confint(reg.res,level=0.95)[-1,1],
                       confint(reg.res,level=0.95)[-1,2],
                       summary(reg.res)$coefficients[-1,4])
  colnames(reg.res.df) <- c("variable", "drug", "coeff",
                            "ci.low", "ci.high","p.value")
  T.LC50.res.df <- rbind(T.LC50.res.df,reg.res.df)
}

## Combine data and plot
B.LC50.res.df$subtype <- rep("B-ALL", dim(B.LC50.res.df)[1])
T.LC50.res.df$subtype <- rep("T-ALL", dim(T.LC50.res.df)[1])
LC50.frac.res.df <- rbind(B.LC50.res.df,T.LC50.res.df)

# Generate plot with only conventional chemotherapy drugs (Fig. 3F)
x.order <- c("Lineage-like", "HSPC-like")
y.order <- c("Asparaginase", "Cytarabine", "Mercaptopurine", "Thioguanine", "Nelarabine", "Daunorubicin",
             "Vincristine", "Dexamethasone","Prednisolone")
LC50.frac.plot.df <- LC50.frac.res.df[LC50.frac.res.df$variable %in% x.order &
                                        LC50.frac.res.df$drug %in% y.order,]
coe <- LC50.frac.plot.df$coeff
LC50.frac.plot.df$color <- (-1*(LC50.frac.plot.df$ci.high<0))+(1*(LC50.frac.plot.df$ci.low>0))
#pdf(file="./Figures/MultipleRegression_ChemoLC50_summed_cellFraction_adj_Age_WBC_byPhenotype.pdf", 8, 4)
ggplot(LC50.frac.plot.df, aes(y=factor(variable, level=x.order), x=drug)) +
  geom_point(aes(size=-log10(p.value), fill=coeff,
                 colour=factor(color), stroke=factor(color)), pch=21) +
  scale_fill_gradientn(colors = c("#00008b","#0000ff","#ffffff", "#ff0000"),
                       values=rescale(c(min(coe)-0.1,-1*max(coe),0,max(coe)),from=c(min(coe),max(coe)))) +
  scale_size(breaks=c(0.5,1,1.5,2),range=c(1,6)) +
  scale_color_manual(values=c("-1"="#000000", "0"="#808080", "1"="#8b0000"), guide="none") +
  scale_discrete_manual(aesthetics="stroke",values=c("-1"=2, "0"=0.5, "1"=2), breaks=c("-1"),
                        labels=c("p<0.05"), guide=guide_legend("")) +
  facet_grid(subtype~., scale="free", space="free") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), legend.spacing.y = unit(1,"mm")) +
  guides(fill=guide_colorbar("coefficient", order=1),size=guide_legend(order=2)) +
  labs(x = "", y = "", title="LC50 Multiple Regression")
#dev.off()

# Generate plot with targeted drugs (Fig. 6G)
x.order <- c("Lineage-like", "HSPC-like")
y.order <- c("Bortezomib", "CHZ868", "Dasatinib", "Ibrutinib", "Trametinib", "Ruxolitinib", "Venetoclax",
             "Panobinostat", "Vorinostat")
LC50.frac.plot.df <- LC50.frac.res.df[LC50.frac.res.df$variable %in% x.order &
                                        LC50.frac.res.df$drug %in% y.order,]
coe <- LC50.frac.plot.df$coeff
LC50.frac.plot.df$color <- (-1*(LC50.frac.plot.df$ci.high<0))+(1*(LC50.frac.plot.df$ci.low>0))
#pdf(file="./Figures/MultipleRegression_TargetLC50_summed_cellFraction_adj_Age_WBC_byPhenotype.pdf", 8, 4)
ggplot(LC50.frac.plot.df, aes(y=factor(variable, level=x.order), x=drug)) +
  geom_point(aes(size=-log10(p.value), fill=coeff,
                 colour=factor(color), stroke=factor(color)), pch=21) +
  scale_fill_gradientn(colors = c("#00008b","#0000ff","#ffffff", "#ff0000"),
                       values=rescale(c(min(coe)-0.1,-1*max(coe),0,max(coe)),from=c(min(coe),max(coe)))) +
  scale_size(breaks=c(0.5,1,1.5,2),range=c(1,6)) +
  scale_color_manual(values=c("-1"="#000000", "0"="#808080", "1"="#8b0000"), guide="none") +
  scale_discrete_manual(aesthetics="stroke",values=c("-1"=2, "0"=0.5, "1"=2), breaks=c("-1"),
                        labels=c("p<0.05"), guide=guide_legend("")) +
  facet_grid(subtype~., scale="free", space="free") +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), legend.spacing.y = unit(1,"mm")) +
  guides(fill=guide_colorbar("coefficient", order=1),size=guide_legend(order=2)) +
  labs(x = "", y = "", title="LC50 Multiple Regression")
#dev.off()
