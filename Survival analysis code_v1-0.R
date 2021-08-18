
# Version control

# v0-1:  First draft
# v0-2:  Minor miscellaneous edits
# v0-3:  Added in re-based curves for HR and piecewise approaches
# v0-4:  Re-formatted code to allow for alternative specification of cut-point
# v0-5:  Added in quantitative outputs (RMST and point estimate)
# v0-6:  Minor miscellaneous edits
# v0-7:  Added 95% CIs for Kaplan-Meier curves
# v0-8:  Edits to plots for legends with CI limits etc.
# v0-9:  Added code for 95% CI limits for Cox PH models
# v0-10: Adding in survHE Bayesian analysis
# v0-11: Adding in Bayesian muhaz plots
# v0-12: Adjustments to Bayesian models, plus formatting (e.g. N@risk for KM)
# v0-13: General tidying up of code, added some extra comments
# v1-0: Added instructions for use by others, and clarified where file directories are specified

# User inputs # ----

# Set working directory and read in .csv containing pseudo-IPD # ----
# Input directory
setwd("C:/Users/...") #<< Here, set the working directory for where you have saved the input file
# Input file
survival.data <- read.csv(file="pILD_v1-0.csv")

# Output directory
setwd("C:/Users/...") #<< Here, set the working directory for where you want to save outputs

# Set cut point
piecewise_cutpoint <- 12.8/12 #12.8 months, de Bono et al.

# Specify an error margin if re-based models yield errors due to survival times close to zero being included
error_margin <- 0.0001
#error_margin <- 0

# Choose whether to use the rebased models or not
use_rebased <- 1

# Reading in data and loading packages # ----
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c(
  "survival",  "flexsurv",  "survminer", "muhaz", "tidyverse", "scales", "PWEALL"
)

ipak(packages)

# devtools::install_github("giabaio/survHE", ref = "devel")
library("survHE")

# Trim data set to only include relevant studies # ----

survival.data.similar <- subset(survival.data,
                                Study == "AFFIRM" |
                                  Study == "COU-AA-301" |
                                  Study == "PREVAIL" |
                                  Study == "ALSYMPCA" |
                                  Study == "TROPIC" |
                                  Study == "TAX 327")

# Clean data set to remove unnecessary levels and record number of different arms available # ----
strata <- as.vector(unique(survival.data.similar$Strata))
number_of_arms <- length(strata)
# For the 22 potentially-relevant study arms, split out the data so each arm can be considered separately # ----

survival.data.by.arm <- list()
survival.data.by.arm[[01]] <- subset(survival.data.similar,Strata == strata[1])  # 1               AFFIRM (Scher 2012) - Enzalutamide
survival.data.by.arm[[02]] <- subset(survival.data.similar,Strata == strata[2])  # 2                    AFFIRM (Scher 2012) - Placebo
survival.data.by.arm[[03]] <- subset(survival.data.similar,Strata == strata[3])  # 3                 ALSYMPCA (Parker 2013) - Placebo
survival.data.by.arm[[04]] <- subset(survival.data.similar,Strata == strata[4])  # 4              ALSYMPCA (Parker 2013) - Radium-223
survival.data.by.arm[[05]] <- subset(survival.data.similar,Strata == strata[5])  # 5  COU-AA-301 (de Bono 2011) - AA*
survival.data.by.arm[[06]] <- subset(survival.data.similar,Strata == strata[6])  # 6              COU-AA-301 (de Bono 2011) - Placebo*
survival.data.by.arm[[07]] <- subset(survival.data.similar,Strata == strata[7])  # 7   COU-AA-301 (Fizazi 2012) - AA*
survival.data.by.arm[[08]] <- subset(survival.data.similar,Strata == strata[8])  # 8               COU-AA-301 (Fizazi 2012) - Placebo*
survival.data.by.arm[[09]] <- subset(survival.data.similar,Strata == strata[9])  # 9               PREVAIL (Beer 2014) - Enzalutamide**
survival.data.by.arm[[10]] <- subset(survival.data.similar,Strata == strata[10]) # 10                   PREVAIL (Beer 2014) - Placebo**
survival.data.by.arm[[11]] <- subset(survival.data.similar,Strata == strata[11]) # 11              PREVAIL (Beer 2017) - Enzalutamide**
survival.data.by.arm[[12]] <- subset(survival.data.similar,Strata == strata[12]) # 12                   PREVAIL (Beer 2017) - Placebo**
survival.data.by.arm[[13]] <- subset(survival.data.similar,Strata == strata[13]) # 13  TAX 327 (Berthold 2008) - Docetaxel every 3 wk
survival.data.by.arm[[14]] <- subset(survival.data.similar,Strata == strata[14]) # 14          TAX 327 (Berthold 2008) - Mitoxantrone
survival.data.by.arm[[15]] <- subset(survival.data.similar,Strata == strata[15]) # 15      TAX 327 (Berthold 2008) - Weekly docetaxel
survival.data.by.arm[[16]] <- subset(survival.data.similar,Strata == strata[16]) # 16   TAX 327 (Tannock 2004) - Docetaxel every 3 wk
survival.data.by.arm[[17]] <- subset(survival.data.similar,Strata == strata[17]) # 17           TAX 327 (Tannock 2004) - Mitoxantrone
survival.data.by.arm[[18]] <- subset(survival.data.similar,Strata == strata[18]) # 18       TAX 327 (Tannock 2004) - Weekly docetaxel
survival.data.by.arm[[19]] <- subset(survival.data.similar,Strata == strata[19]) # 19                TROPIC (Bahl 2012) - Cabazitaxel
survival.data.by.arm[[20]] <- subset(survival.data.similar,Strata == strata[20]) # 20               TROPIC (Bahl 2012) - Mitoxantrone
survival.data.by.arm[[21]] <- subset(survival.data.similar,Strata == strata[21]) # 21             TROPIC (de Bono 2010) - Cabazitaxel
survival.data.by.arm[[22]] <- subset(survival.data.similar,Strata == strata[22]) # 22            TROPIC (de Bono 2010) - Mitoxantrone

survival.data.by.arm.rebased <- list()
survival.data.by.arm.rebased[[01]] <- subset(survival.data.by.arm[[01]], Years >= piecewise_cutpoint + error_margin)  # 1               AFFIRM (Scher 2012) - Enzalutamide
survival.data.by.arm.rebased[[02]] <- subset(survival.data.by.arm[[02]], Years >= piecewise_cutpoint + error_margin)  # 2                    AFFIRM (Scher 2012) - Placebo
survival.data.by.arm.rebased[[03]] <- subset(survival.data.by.arm[[03]], Years >= piecewise_cutpoint + error_margin)  # 3                 ALSYMPCA (Parker 2013) - Placebo
survival.data.by.arm.rebased[[04]] <- subset(survival.data.by.arm[[04]], Years >= piecewise_cutpoint + error_margin)  # 4              ALSYMPCA (Parker 2013) - Radium-223
survival.data.by.arm.rebased[[05]] <- subset(survival.data.by.arm[[05]], Years >= piecewise_cutpoint + error_margin)  # 5  COU-AA-301 (de Bono 2011) - AA*
survival.data.by.arm.rebased[[06]] <- subset(survival.data.by.arm[[06]], Years >= piecewise_cutpoint + error_margin)  # 6              COU-AA-301 (de Bono 2011) - Placebo*
survival.data.by.arm.rebased[[07]] <- subset(survival.data.by.arm[[07]], Years >= piecewise_cutpoint + error_margin)  # 7   COU-AA-301 (Fizazi 2012) - AA*
survival.data.by.arm.rebased[[08]] <- subset(survival.data.by.arm[[08]], Years >= piecewise_cutpoint + error_margin)  # 8               COU-AA-301 (Fizazi 2012) - Placebo*
survival.data.by.arm.rebased[[09]] <- subset(survival.data.by.arm[[09]], Years >= piecewise_cutpoint + error_margin)  # 9               PREVAIL (Beer 2014) - Enzalutamide**
survival.data.by.arm.rebased[[10]] <- subset(survival.data.by.arm[[10]], Years >= piecewise_cutpoint + error_margin) # 10                   PREVAIL (Beer 2014) - Placebo**
survival.data.by.arm.rebased[[11]] <- subset(survival.data.by.arm[[11]], Years >= piecewise_cutpoint + error_margin) # 11              PREVAIL (Beer 2017) - Enzalutamide**
survival.data.by.arm.rebased[[12]] <- subset(survival.data.by.arm[[12]], Years >= piecewise_cutpoint + error_margin) # 12                   PREVAIL (Beer 2017) - Placebo**
survival.data.by.arm.rebased[[13]] <- subset(survival.data.by.arm[[13]], Years >= piecewise_cutpoint + error_margin) # 13  TAX 327 (Berthold 2008) - Docetaxel every 3 wk
survival.data.by.arm.rebased[[14]] <- subset(survival.data.by.arm[[14]], Years >= piecewise_cutpoint + error_margin) # 14          TAX 327 (Berthold 2008) - Mitoxantrone
survival.data.by.arm.rebased[[15]] <- subset(survival.data.by.arm[[15]], Years >= piecewise_cutpoint + error_margin) # 15      TAX 327 (Berthold 2008) - Weekly docetaxel
survival.data.by.arm.rebased[[16]] <- subset(survival.data.by.arm[[16]], Years >= piecewise_cutpoint + error_margin) # 16   TAX 327 (Tannock 2004) - Docetaxel every 3 wk
survival.data.by.arm.rebased[[17]] <- subset(survival.data.by.arm[[17]], Years >= piecewise_cutpoint + error_margin) # 17           TAX 327 (Tannock 2004) - Mitoxantrone
survival.data.by.arm.rebased[[18]] <- subset(survival.data.by.arm[[18]], Years >= piecewise_cutpoint + error_margin) # 18       TAX 327 (Tannock 2004) - Weekly docetaxel
survival.data.by.arm.rebased[[19]] <- subset(survival.data.by.arm[[19]], Years >= piecewise_cutpoint + error_margin) # 19                TROPIC (Bahl 2012) - Cabazitaxel
survival.data.by.arm.rebased[[20]] <- subset(survival.data.by.arm[[20]], Years >= piecewise_cutpoint + error_margin) # 20               TROPIC (Bahl 2012) - Mitoxantrone
survival.data.by.arm.rebased[[21]] <- subset(survival.data.by.arm[[21]], Years >= piecewise_cutpoint + error_margin) # 21             TROPIC (de Bono 2010) - Cabazitaxel
survival.data.by.arm.rebased[[22]] <- subset(survival.data.by.arm[[22]], Years >= piecewise_cutpoint + error_margin) # 22            TROPIC (de Bono 2010) - Mitoxantrone

survival.data.by.arm.rebased[[01]]$Years <- survival.data.by.arm.rebased[[01]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[02]]$Years <- survival.data.by.arm.rebased[[02]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[03]]$Years <- survival.data.by.arm.rebased[[03]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[04]]$Years <- survival.data.by.arm.rebased[[04]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[05]]$Years <- survival.data.by.arm.rebased[[05]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[06]]$Years <- survival.data.by.arm.rebased[[06]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[07]]$Years <- survival.data.by.arm.rebased[[07]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[08]]$Years <- survival.data.by.arm.rebased[[08]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[09]]$Years <- survival.data.by.arm.rebased[[09]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[10]]$Years <- survival.data.by.arm.rebased[[10]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[11]]$Years <- survival.data.by.arm.rebased[[11]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[12]]$Years <- survival.data.by.arm.rebased[[12]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[13]]$Years <- survival.data.by.arm.rebased[[13]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[14]]$Years <- survival.data.by.arm.rebased[[14]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[15]]$Years <- survival.data.by.arm.rebased[[15]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[16]]$Years <- survival.data.by.arm.rebased[[16]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[17]]$Years <- survival.data.by.arm.rebased[[17]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[18]]$Years <- survival.data.by.arm.rebased[[18]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[19]]$Years <- survival.data.by.arm.rebased[[19]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[20]]$Years <- survival.data.by.arm.rebased[[20]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[21]]$Years <- survival.data.by.arm.rebased[[21]]$Years - piecewise_cutpoint
survival.data.by.arm.rebased[[22]]$Years <- survival.data.by.arm.rebased[[22]]$Years - piecewise_cutpoint
  
#Notes (for above):
#* This is the main trial
#** This is the PREVAIL trial that was deemed inappropriate for consideration

# Set up specific colours for plots (used throughout) # ----

plot_red <- rgb(248, 118, 109, max=255)
plot_pink <- rgb(251, 97, 215, max=255)
plot_purple <- rgb(165, 138, 255, max=255)
plot_lightblue <- rgb(0, 182, 235, max=255)
plot_tropicalgreen <- rgb(0, 192, 148, max=255)
plot_forestgreen <- rgb(83, 180, 0, max=255)
plot_bronze <- rgb(196, 154, 0, max=255)

plot_red_2 <- rgb(248, 118, 109.1, max=255)
plot_pink_2 <- rgb(251, 97, 215.1, max=255)
plot_purple_2 <- rgb(165, 138, 254.9, max=255)
plot_lightblue_2 <- rgb(0, 182, 235.1, max=255)
plot_tropicalgreen_2 <- rgb(0, 192, 148.1, max=255)
plot_forestgreen_2 <- rgb(83, 180, 0.1, max=255)
plot_bronze_2 <- rgb(196, 154, 0.1, max=255)

# Near duplicates produced to avoid errors with ggplot reading the same colour in twice

# Relevant KM plot for active arms # ----

# Relevant studies
survival.data.similar.active <- subset(survival.data,
                                       Study == "AFFIRM" |
                                       Study == "COU-AA-301" |
                                       Study == "PREVAIL" |
                                       Study == "ALSYMPCA" |
                                       Study == "TROPIC" |
                                       Study == "TAX 327")
# Remove placebo and mitoxantrone arms
survival.data.similar.active <- subset(survival.data.similar.active, Arm != "Placebo")
survival.data.similar.active <- subset(survival.data.similar.active, Arm != "Mitoxantrone")

# Take the latest data cut (except for COU-AA-301)
survival.data.similar.active <- subset(survival.data.similar.active,
                                       Source == "Beer 2017" |
                                       Source == "Bahl 2012" |
                                       Source == "de Bono 2011" | #*This is the earlier data cut for COU-AA-301
                                       Source == "Scher 2012" |
                                       Source == "Parker 2013" |
                                       Source == "Berthold 2008")

KM.est.active <- survfit(Surv(Years,Event) ~ Strata, data = survival.data.similar.active, type = "kaplan-meier")
KM_active <- ggsurvplot(KM.est.active,
                     survival.data.similar.active,
                     palette = c(plot_red,plot_bronze,"black",plot_tropicalgreen,plot_lightblue,plot_lightblue_2,plot_pink),
                     linetype = c(1,1,1,1,1,2,1),
                     break.x.by = 1,
                     break.y.by = 0.1,
                     xlim=c(0,7),
                     surv.scale="percent",
                     ylab="Overall survival",
                     xlab="Time (years)",
                     font.x= c(10, "bold", "black"),
                     font.y= c(10, "bold", "black"),
                     font.tickslab= c(10, "plain", "black"),
                     font.legend = c(8,"plain","black"),
                     censor=T,
                     censor.size = 2,
                     legend="bottom",
                     legend.labs=c("AFFIRM (Scher 2012) - Enzalutamide",
                                   "ALSYMPCA (Parker 2013) - Radium-223",
                                   "COU-AA-301 (de Bono 2011) - AA",
                                   "PREVAIL (Beer 2017) - Enzalutamide",
                                   "TAX 327 (Berthold 2008) - Docetaxel every 3 wk",
                                   "TAX 327 (Berthold 2008) - Weekly docextaxel",
                                   "TROPIC (Bahl 2012) - Cabazitaxel")
) + guides(colour = guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4))

survival.data.pivotal.original <- subset(survival.data, Study == "COU-AA-301")
survival.data.pivotal.original <- subset(survival.data.pivotal.original, Source == "de Bono 2011")

KM.est.original.AA <- survfit(Surv(Years,Event) ~ Strata, data = survival.data.pivotal.original, type = "kaplan-meier")

KM_original_AA <- ggsurvplot(KM.est.original.AA,
                             survival.data.pivotal.original,
                             palette = c(plot_red,plot_lightblue),
                             break.x.by = 0.25,
                             break.y.by = .1,
                             lty=c(1),
                             xlim=c(0,2),
                             surv.scale="percent",
                             ylab="Overall survival",
                             xlab="Time (years)",
                             font.x= c(10, "bold", "black"),
                             font.y= c(10, "bold", "black"),
                             font.tickslab= c(10, "plain", "black"),
                             font.legend = c(8,"plain","black"),
                             censor=T,
                             risk.table = TRUE,
                             risk.table.y.text = FALSE,
                             tables.theme = theme_cleantable(),
                             #axes.offset = F,
                             conf.int=T,
                             fontsize=3.5,
                             censor.size = 2,
                             legend="bottom",
                             legend.labs=c("COU-AA-301 (de Bono 2011) - AA",
                                           "COU-AA-301 (de Bono 2011) - Placebo")
) + guides(colour = guide_legend(title="", nrow=1), fill = guide_legend(title="", nrow=1))

KM_original_AA$plot <- KM_original_AA$plot + scale_x_continuous(labels = scales::number_format(scale = 1, accuracy = 0.01, decimal.mark = '.', suffix = ""), breaks = c(seq(0,2,0.25)))
KM_original_AA$table <- KM_original_AA$table + theme(plot.title = element_text(size = 8, color = "black", face = "bold"))
KM_original_AA$table <- KM_original_AA$table + theme(axis.text.y = element_text(size = 8, face = "bold"))
KM_original_AA$plot <- KM_original_AA$plot + scale_y_continuous(labels = scales::number_format(scale = 100, accuracy = 1, decimal.mark = '.', suffix = "%"), breaks = c(seq(0,1,0.1)))

# Relevant KM plot for placebo arms # ----

# Relevant studies
survival.data.similar.placebo <- subset(survival.data,
                                        Study == "AFFIRM" |
                                        Study == "COU-AA-301" |
                                        Study == "PREVAIL" |
                                        Study == "ALSYMPCA" |
                                        Study == "TROPIC" |
                                        Study == "TAX 327")

# Use only placebo and mitoxantrone arms
survival.data.similar.placebo <- subset(survival.data.similar.placebo,
                                        Arm == "Placebo" |
                                        Arm == "Mitoxantrone")

# Take the latest data cut (except for COU-AA-301)
survival.data.similar.placebo <- subset(survival.data.similar.placebo,
                                        Source == "Beer 2017" |
                                        Source == "Bahl 2012" |
                                        Source == "de Bono 2011" |
                                        Source == "Scher 2012" |
                                        Source == "Parker 2013" |
                                        Source == "Berthold 2008")

KM.est.placebo <- survfit(Surv(Years,Event) ~ Strata, data = survival.data.similar.placebo, type = "kaplan-meier")
KM_placebo <- ggsurvplot(KM.est.placebo,
                      survival.data.similar.placebo,
                      palette = c(plot_red,plot_bronze,"black",plot_tropicalgreen,plot_lightblue,plot_pink),
                      linetype = c(1),
                      break.x.by = 1,
                      break.y.by = 0.1,
                      xlim=c(0,7),
                      surv.scale="percent",
                      ylab="Overall survival",
                      xlab="Time (years)",
                      font.x= c(10, "bold", "black"),
                      font.y= c(10, "bold", "black"),
                      font.tickslab= c(10, "plain", "black"),
                      font.legend = c(8,"plain","black"),
                      censor=T,
                      censor.size = 2,
                      legend="bottom",
                      legend.labs = c("AFFIRM (Scher 2012) - Placebo",
                                      "ALSYMPCA (Parker 2013) - Placebo",
                                      "COU-AA-301 (de Bono 2011) - Placebo",
                                      "PREVAIL (Beer 2017) - Placebo",
                                      "TAX 327 (Berthold 2008) - Mitoxantrone",
                                      "TROPIC (de Bono 2010) - Mitoxantrone")
) +  guides(colour = guide_legend(title="", nrow=4))

# Produce smoothed hazard estimates via muhaz # ----

# Multiplier set here in case wanting to explore similarity in shapes of curves
# Set to 1 in base-case analysis
muhaz_mult <- 1

hazTemp <- list()
hazTemp[[1]] = muhaz(time = survival.data.by.arm[[1]]$Years, delta = survival.data.by.arm[[1]]$Event, max.time = max(survival.data.by.arm[[1]]$Years)*muhaz_mult)
hazTemp[[2]] = muhaz(time = survival.data.by.arm[[2]]$Years, delta = survival.data.by.arm[[2]]$Event, max.time = max(survival.data.by.arm[[2]]$Years)*muhaz_mult)
hazTemp[[3]] = muhaz(time = survival.data.by.arm[[3]]$Years, delta = survival.data.by.arm[[3]]$Event, max.time = max(survival.data.by.arm[[3]]$Years)*muhaz_mult)
hazTemp[[4]] = muhaz(time = survival.data.by.arm[[4]]$Years, delta = survival.data.by.arm[[4]]$Event, max.time = max(survival.data.by.arm[[4]]$Years)*muhaz_mult)
hazTemp[[5]] = muhaz(time = survival.data.by.arm[[5]]$Years, delta = survival.data.by.arm[[5]]$Event, max.time = max(survival.data.by.arm[[5]]$Years)*muhaz_mult)
hazTemp[[6]] = muhaz(time = survival.data.by.arm[[6]]$Years, delta = survival.data.by.arm[[6]]$Event, max.time = max(survival.data.by.arm[[6]]$Years)*muhaz_mult)
hazTemp[[7]] = muhaz(time = survival.data.by.arm[[7]]$Years, delta = survival.data.by.arm[[7]]$Event, max.time = max(survival.data.by.arm[[7]]$Years)*muhaz_mult)
hazTemp[[8]] = muhaz(time = survival.data.by.arm[[8]]$Years, delta = survival.data.by.arm[[8]]$Event, max.time = max(survival.data.by.arm[[8]]$Years)*muhaz_mult)
hazTemp[[9]] = muhaz(time = survival.data.by.arm[[9]]$Years, delta = survival.data.by.arm[[9]]$Event, max.time = max(survival.data.by.arm[[9]]$Years)*muhaz_mult)
hazTemp[[10]] = muhaz(time = survival.data.by.arm[[10]]$Years, delta = survival.data.by.arm[[10]]$Event, max.time = max(survival.data.by.arm[[10]]$Years)*muhaz_mult)
hazTemp[[11]] = muhaz(time = survival.data.by.arm[[11]]$Years, delta = survival.data.by.arm[[11]]$Event, max.time = max(survival.data.by.arm[[11]]$Years)*muhaz_mult)
hazTemp[[12]] = muhaz(time = survival.data.by.arm[[12]]$Years, delta = survival.data.by.arm[[12]]$Event, max.time = max(survival.data.by.arm[[12]]$Years)*muhaz_mult)
hazTemp[[13]] = muhaz(time = survival.data.by.arm[[13]]$Years, delta = survival.data.by.arm[[13]]$Event, max.time = max(survival.data.by.arm[[13]]$Years)*muhaz_mult)
hazTemp[[14]] = muhaz(time = survival.data.by.arm[[14]]$Years, delta = survival.data.by.arm[[14]]$Event, max.time = max(survival.data.by.arm[[14]]$Years)*muhaz_mult)
hazTemp[[15]] = muhaz(time = survival.data.by.arm[[15]]$Years, delta = survival.data.by.arm[[15]]$Event, max.time = max(survival.data.by.arm[[15]]$Years)*muhaz_mult)
hazTemp[[16]] = muhaz(time = survival.data.by.arm[[16]]$Years, delta = survival.data.by.arm[[16]]$Event, max.time = max(survival.data.by.arm[[16]]$Years)*muhaz_mult)
hazTemp[[17]] = muhaz(time = survival.data.by.arm[[17]]$Years, delta = survival.data.by.arm[[17]]$Event, max.time = max(survival.data.by.arm[[17]]$Years)*muhaz_mult)
hazTemp[[18]] = muhaz(time = survival.data.by.arm[[18]]$Years, delta = survival.data.by.arm[[18]]$Event, max.time = max(survival.data.by.arm[[18]]$Years)*muhaz_mult)
hazTemp[[19]] = muhaz(time = survival.data.by.arm[[19]]$Years, delta = survival.data.by.arm[[19]]$Event, max.time = max(survival.data.by.arm[[19]]$Years)*muhaz_mult)
hazTemp[[20]] = muhaz(time = survival.data.by.arm[[20]]$Years, delta = survival.data.by.arm[[20]]$Event, max.time = max(survival.data.by.arm[[20]]$Years)*muhaz_mult)
hazTemp[[21]] = muhaz(time = survival.data.by.arm[[21]]$Years, delta = survival.data.by.arm[[21]]$Event, max.time = max(survival.data.by.arm[[21]]$Years)*muhaz_mult)
hazTemp[[22]] = muhaz(time = survival.data.by.arm[[22]]$Years, delta = survival.data.by.arm[[22]]$Event, max.time = max(survival.data.by.arm[[22]]$Years)*muhaz_mult)

hazTemp.rebased <- list()
hazTemp.rebased[[1]] = muhaz(time = survival.data.by.arm.rebased[[1]]$Years, delta = survival.data.by.arm.rebased[[1]]$Event, max.time = max(survival.data.by.arm.rebased[[1]]$Years)*muhaz_mult)
hazTemp.rebased[[2]] = muhaz(time = survival.data.by.arm.rebased[[2]]$Years, delta = survival.data.by.arm.rebased[[2]]$Event, max.time = max(survival.data.by.arm.rebased[[2]]$Years)*muhaz_mult)
hazTemp.rebased[[3]] = muhaz(time = survival.data.by.arm.rebased[[3]]$Years, delta = survival.data.by.arm.rebased[[3]]$Event, max.time = max(survival.data.by.arm.rebased[[3]]$Years)*muhaz_mult)
hazTemp.rebased[[4]] = muhaz(time = survival.data.by.arm.rebased[[4]]$Years, delta = survival.data.by.arm.rebased[[4]]$Event, max.time = max(survival.data.by.arm.rebased[[4]]$Years)*muhaz_mult)
hazTemp.rebased[[5]] = muhaz(time = survival.data.by.arm.rebased[[5]]$Years, delta = survival.data.by.arm.rebased[[5]]$Event, max.time = max(survival.data.by.arm.rebased[[5]]$Years)*muhaz_mult)
hazTemp.rebased[[6]] = muhaz(time = survival.data.by.arm.rebased[[6]]$Years, delta = survival.data.by.arm.rebased[[6]]$Event, max.time = max(survival.data.by.arm.rebased[[6]]$Years)*muhaz_mult)
hazTemp.rebased[[7]] = muhaz(time = survival.data.by.arm.rebased[[7]]$Years, delta = survival.data.by.arm.rebased[[7]]$Event, max.time = max(survival.data.by.arm.rebased[[7]]$Years)*muhaz_mult)
hazTemp.rebased[[8]] = muhaz(time = survival.data.by.arm.rebased[[8]]$Years, delta = survival.data.by.arm.rebased[[8]]$Event, max.time = max(survival.data.by.arm.rebased[[8]]$Years)*muhaz_mult)
hazTemp.rebased[[9]] = muhaz(time = survival.data.by.arm.rebased[[9]]$Years, delta = survival.data.by.arm.rebased[[9]]$Event, max.time = max(survival.data.by.arm.rebased[[9]]$Years)*muhaz_mult)
hazTemp.rebased[[10]] = muhaz(time = survival.data.by.arm.rebased[[10]]$Years, delta = survival.data.by.arm.rebased[[10]]$Event, max.time = max(survival.data.by.arm.rebased[[10]]$Years)*muhaz_mult)
hazTemp.rebased[[11]] = muhaz(time = survival.data.by.arm.rebased[[11]]$Years, delta = survival.data.by.arm.rebased[[11]]$Event, max.time = max(survival.data.by.arm.rebased[[11]]$Years)*muhaz_mult)
hazTemp.rebased[[12]] = muhaz(time = survival.data.by.arm.rebased[[12]]$Years, delta = survival.data.by.arm.rebased[[12]]$Event, max.time = max(survival.data.by.arm.rebased[[12]]$Years)*muhaz_mult)
hazTemp.rebased[[13]] = muhaz(time = survival.data.by.arm.rebased[[13]]$Years, delta = survival.data.by.arm.rebased[[13]]$Event, max.time = max(survival.data.by.arm.rebased[[13]]$Years)*muhaz_mult)
hazTemp.rebased[[14]] = muhaz(time = survival.data.by.arm.rebased[[14]]$Years, delta = survival.data.by.arm.rebased[[14]]$Event, max.time = max(survival.data.by.arm.rebased[[14]]$Years)*muhaz_mult)
hazTemp.rebased[[15]] = muhaz(time = survival.data.by.arm.rebased[[15]]$Years, delta = survival.data.by.arm.rebased[[15]]$Event, max.time = max(survival.data.by.arm.rebased[[15]]$Years)*muhaz_mult)
hazTemp.rebased[[16]] = muhaz(time = survival.data.by.arm.rebased[[16]]$Years, delta = survival.data.by.arm.rebased[[16]]$Event, max.time = max(survival.data.by.arm.rebased[[16]]$Years)*muhaz_mult)
hazTemp.rebased[[17]] = muhaz(time = survival.data.by.arm.rebased[[17]]$Years, delta = survival.data.by.arm.rebased[[17]]$Event, max.time = max(survival.data.by.arm.rebased[[17]]$Years)*muhaz_mult)
hazTemp.rebased[[18]] = muhaz(time = survival.data.by.arm.rebased[[18]]$Years, delta = survival.data.by.arm.rebased[[18]]$Event, max.time = max(survival.data.by.arm.rebased[[18]]$Years)*muhaz_mult)
hazTemp.rebased[[19]] = muhaz(time = survival.data.by.arm.rebased[[19]]$Years, delta = survival.data.by.arm.rebased[[19]]$Event, max.time = max(survival.data.by.arm.rebased[[19]]$Years)*muhaz_mult)
hazTemp.rebased[[20]] = muhaz(time = survival.data.by.arm.rebased[[20]]$Years, delta = survival.data.by.arm.rebased[[20]]$Event, max.time = max(survival.data.by.arm.rebased[[20]]$Years)*muhaz_mult)
hazTemp.rebased[[21]] = muhaz(time = survival.data.by.arm.rebased[[21]]$Years, delta = survival.data.by.arm.rebased[[21]]$Event, max.time = max(survival.data.by.arm.rebased[[21]]$Years)*muhaz_mult)
hazTemp.rebased[[22]] = muhaz(time = survival.data.by.arm.rebased[[22]]$Years, delta = survival.data.by.arm.rebased[[22]]$Event, max.time = max(survival.data.by.arm.rebased[[22]]$Years)*muhaz_mult)

# Relevant muhaz plot for active arms # ----

# Produce a list of the relevant curves for active arms
# The relevant arms are marked with '++++' as follows:

# 1               AFFIRM (Scher 2012) - Enzalutamide ++++
# 2                    AFFIRM (Scher 2012) - Placebo
# 3                 ALSYMPCA (Parker 2013) - Placebo
# 4              ALSYMPCA (Parker 2013) - Radium-223 ++++ 
# 5  COU-AA-301 (de Bono 2011) - AA ++++
# 6              COU-AA-301 (de Bono 2011) - Placebo
# 7   COU-AA-301 (Fizazi 2012) - AA
# 8               COU-AA-301 (Fizazi 2012) - Placebo
# 9               PREVAIL (Beer 2014) - Enzalutamide
# 10                   PREVAIL (Beer 2014) - Placebo
# 11              PREVAIL (Beer 2017) - Enzalutamide ++++
# 12                   PREVAIL (Beer 2017) - Placebo
# 13  TAX 327 (Berthold 2008) - Docetaxel every 3 wk ++++
# 14          TAX 327 (Berthold 2008) - Mitoxantrone
# 15      TAX 327 (Berthold 2008) - Weekly docetaxel ++++
# 16   TAX 327 (Tannock 2004) - Docetaxel every 3 wk
# 17           TAX 327 (Tannock 2004) - Mitoxantrone
# 18       TAX 327 (Tannock 2004) - Weekly docetaxel
# 19                TROPIC (Bahl 2012) - Cabazitaxel ++++
# 20               TROPIC (Bahl 2012) - Mitoxantrone
# 21             TROPIC (de Bono 2010) - Cabazitaxel
# 22            TROPIC (de Bono 2010) - Mitoxantrone

# Put items into a list
active.muhaz.plotting <- list()
active.muhaz.plotting[[01]] <- data.frame(Hazard = hazTemp[[01]]$haz.est, Time = hazTemp[[01]]$est.grid, Strata = "AFFIRM (Scher 2012) - Enzalutamide")
active.muhaz.plotting[[04]] <- data.frame(Hazard = hazTemp[[04]]$haz.est, Time = hazTemp[[04]]$est.grid, Strata = "ALSYMPCA (Parker 2013) - Radium-223")
active.muhaz.plotting[[05]] <- data.frame(Hazard = hazTemp[[05]]$haz.est, Time = hazTemp[[05]]$est.grid, Strata = "COU-AA-301 (de Bono 2011) - AA")
active.muhaz.plotting[[11]] <- data.frame(Hazard = hazTemp[[11]]$haz.est, Time = hazTemp[[11]]$est.grid, Strata = "PREVAIL (Beer 2017) - Enzalutamide")
active.muhaz.plotting[[13]] <- data.frame(Hazard = hazTemp[[13]]$haz.est, Time = hazTemp[[13]]$est.grid, Strata = "TAX 327 (Berthold 2008) - Docetaxel every 3 weeks")
active.muhaz.plotting[[15]] <- data.frame(Hazard = hazTemp[[15]]$haz.est, Time = hazTemp[[15]]$est.grid, Strata = "TAX 327 (Berthold 2008) - Weekly docetaxel")
active.muhaz.plotting[[19]] <- data.frame(Hazard = hazTemp[[19]]$haz.est, Time = hazTemp[[19]]$est.grid, Strata = "TROPIC (Bahl 2012) - Cabazitaxel")

# Combine into a single dataset
haztempplot.active <- rbind(active.muhaz.plotting[[01]],
                            active.muhaz.plotting[[04]],
                            active.muhaz.plotting[[05]],
                            active.muhaz.plotting[[11]],
                            active.muhaz.plotting[[13]],
                            active.muhaz.plotting[[15]],
                            active.muhaz.plotting[[19]])

# Ensure levels are ordered correctly
summary(haztempplot.active$Strata)
haztempplot.active$Strata <- factor(haztempplot.active$Strata)
summary(haztempplot.active$Strata)
haztempplot.active$Strata <- factor(haztempplot.active$Strata,levels(haztempplot.active$Strata)[c(1,2,3,4,5,6,7)])

# Produce list of relevant colours and line types for plotting
active_colours_list <- c("AFFIRM (Scher 2012) - Enzalutamide" = plot_red,
                         "ALSYMPCA (Parker 2013) - Radium-223" = plot_bronze,
                         "COU-AA-301 (de Bono 2011) - AA" = "black",
                         "PREVAIL (Beer 2017) - Enzalutamide" = plot_tropicalgreen,
                         "TAX 327 (Berthold 2008) - Docetaxel every 3 weeks" = plot_lightblue,
                         "TAX 327 (Berthold 2008) - Weekly docetaxel" = plot_lightblue_2,
                         "TROPIC (Bahl 2012) - Cabazitaxel" = plot_pink)

active_lines_list <- c("AFFIRM (Scher 2012) - Enzalutamide" = "solid",
                       "ALSYMPCA (Parker 2013) - Radium-223" = "solid",
                       "COU-AA-301 (de Bono 2011) - AA" = "solid",
                       "PREVAIL (Beer 2017) - Enzalutamide" = "solid",
                       "TAX 327 (Berthold 2008) - Docetaxel every 3 weeks" = "solid",
                       "TAX 327 (Berthold 2008) - Weekly docetaxel" = "dashed",
                       "TROPIC (Bahl 2012) - Cabazitaxel" = "solid")

# Produce plot
muhaz_active <- ggplot(haztempplot.active, aes(x=Time, y=Hazard, colour= Strata, linetype = Strata)) +
  geom_line(size=1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,4), oob = rescale_none, breaks=c(0,1,2,3,4)) +
  scale_x_continuous(limits = c(0,7), oob = rescale_none, breaks=c(0,1,2,3,4,5,6,7)) +
  scale_colour_manual(values=active_colours_list)+
  scale_linetype_manual(values=active_lines_list)+
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8))

# Put items into a list
active.muhaz.plotting.rebased <- list()
active.muhaz.plotting.rebased[[01]] <- data.frame(Hazard = hazTemp.rebased[[01]]$haz.est, Time = hazTemp.rebased[[01]]$est.grid, Strata = "AFFIRM (Scher 2012) - Enzalutamide")
active.muhaz.plotting.rebased[[04]] <- data.frame(Hazard = hazTemp.rebased[[04]]$haz.est, Time = hazTemp.rebased[[04]]$est.grid, Strata = "ALSYMPCA (Parker 2013) - Radium-223")
active.muhaz.plotting.rebased[[05]] <- data.frame(Hazard = hazTemp.rebased[[05]]$haz.est, Time = hazTemp.rebased[[05]]$est.grid, Strata = "COU-AA-301 (de Bono 2011) - AA")
active.muhaz.plotting.rebased[[11]] <- data.frame(Hazard = hazTemp.rebased[[11]]$haz.est, Time = hazTemp.rebased[[11]]$est.grid, Strata = "PREVAIL (Beer 2017) - Enzalutamide")
active.muhaz.plotting.rebased[[13]] <- data.frame(Hazard = hazTemp.rebased[[13]]$haz.est, Time = hazTemp.rebased[[13]]$est.grid, Strata = "TAX 327 (Berthold 2008) - Docetaxel every 3 weeks")
active.muhaz.plotting.rebased[[15]] <- data.frame(Hazard = hazTemp.rebased[[15]]$haz.est, Time = hazTemp.rebased[[15]]$est.grid, Strata = "TAX 327 (Berthold 2008) - Weekly docetaxel")
active.muhaz.plotting.rebased[[19]] <- data.frame(Hazard = hazTemp.rebased[[19]]$haz.est, Time = hazTemp.rebased[[19]]$est.grid, Strata = "TROPIC (Bahl 2012) - Cabazitaxel")

# Combine into a single dataset
haztempplot.active.rebased <- rbind(active.muhaz.plotting.rebased[[01]],
                                    active.muhaz.plotting.rebased[[04]],
                                    active.muhaz.plotting.rebased[[05]],
                                    active.muhaz.plotting.rebased[[11]],
                                    active.muhaz.plotting.rebased[[13]],
                                    active.muhaz.plotting.rebased[[15]],
                                    active.muhaz.plotting.rebased[[19]])

# Ensure levels are ordered correctly
summary(haztempplot.active.rebased$Strata)
haztempplot.active.rebased$Strata <- factor(haztempplot.active.rebased$Strata)
summary(haztempplot.active.rebased$Strata)
haztempplot.active.rebased$Strata <- factor(haztempplot.active.rebased$Strata,levels(haztempplot.active.rebased$Strata)[c(1,2,3,4,5,6,7)])

# Produce plot
muhaz_active_rebased <- ggplot(haztempplot.active.rebased, aes(x=Time, y=Hazard, colour= Strata, linetype = Strata)) +
  geom_line(size=1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,4), oob = rescale_none, breaks=c(0,1,2,3,4)) +
  scale_x_continuous(limits = c(0,7), oob = rescale_none, breaks=c(0,1,2,3,4,5,6,7)) +
  scale_colour_manual(values=active_colours_list)+
  scale_linetype_manual(values=active_lines_list)+
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8))

# Combine into a single dataset
haztempplot.data.active <- rbind(survival.data.by.arm[[01]],
                                 survival.data.by.arm[[04]],
                                 survival.data.by.arm[[11]],
                                 survival.data.by.arm[[13]],
                                 survival.data.by.arm[[15]],
                                 survival.data.by.arm[[19]])

haztempplot.data.active.rebased <- rbind(survival.data.by.arm.rebased[[01]],
                                         survival.data.by.arm.rebased[[04]],
                                         survival.data.by.arm.rebased[[11]],
                                         survival.data.by.arm.rebased[[13]],
                                         survival.data.by.arm.rebased[[15]],
                                         survival.data.by.arm.rebased[[19]])

hazTemp.pooled.active = muhaz(time = haztempplot.data.active$Years, delta = haztempplot.data.active$Event, max.time = max(haztempplot.data.active$Years)*muhaz_mult)
hazTemp.pooled.active.plotting <- data.frame(Hazard = hazTemp.pooled.active$haz.est, Time = hazTemp.pooled.active$est.grid, Strata = "External data (pooled)")

hazTemp.pooled.active.rebased = muhaz(time = haztempplot.data.active.rebased$Years, delta = haztempplot.data.active.rebased$Event, max.time = max(haztempplot.data.active.rebased$Years)*muhaz_mult)
hazTemp.pooled.active.rebased.plotting <- data.frame(Hazard = hazTemp.pooled.active.rebased$haz.est, Time = hazTemp.pooled.active.rebased$est.grid, Strata = "External data (pooled, rebased)")

muhaz_active_pooled <- ggplot(hazTemp.pooled.active.plotting, aes(x=Time, y=Hazard)) +
  geom_line(size=1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,7), oob = rescale_none, breaks=c(0,1,2,3,4,5,6,7)) +
  scale_colour_manual(values=active_colours_list)+
  scale_linetype_manual(values=active_lines_list)+
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8))

muhaz_active_pooled_rebased <- ggplot(hazTemp.pooled.active.rebased.plotting, aes(x=Time, y=Hazard)) +
  geom_line(size=1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,7), oob = rescale_none, breaks=c(0,1,2,3,4,5,6,7)) +
  scale_colour_manual(values=active_colours_list)+
  scale_linetype_manual(values=active_lines_list)+
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8))

# Relevant muhaz plot for placebo arms # ----

# Produce a list of the relevant curves for active arms
# The relevant arms are marked with 'xxxx' as follows:

# 1               AFFIRM (Scher 2012) - Enzalutamide
# 2                    AFFIRM (Scher 2012) - Placebo xxxx
# 3                 ALSYMPCA (Parker 2013) - Placebo xxxx
# 4              ALSYMPCA (Parker 2013) - Radium-223
# 5  COU-AA-301 (de Bono 2011) - AA
# 6              COU-AA-301 (de Bono 2011) - Placebo xxxx
# 7   COU-AA-301 (Fizazi 2012) - AA
# 8               COU-AA-301 (Fizazi 2012) - Placebo
# 9               PREVAIL (Beer 2014) - Enzalutamide
# 10                   PREVAIL (Beer 2014) - Placebo
# 11              PREVAIL (Beer 2017) - Enzalutamide
# 12                   PREVAIL (Beer 2017) - Placebo xxxx
# 13  TAX 327 (Berthold 2008) - Docetaxel every 3 wk
# 14          TAX 327 (Berthold 2008) - Mitoxantrone xxxx
# 15      TAX 327 (Berthold 2008) - Weekly docetaxel
# 16   TAX 327 (Tannock 2004) - Docetaxel every 3 wk
# 17           TAX 327 (Tannock 2004) - Mitoxantrone
# 18       TAX 327 (Tannock 2004) - Weekly docetaxel
# 19                TROPIC (Bahl 2012) - Cabazitaxel
# 20               TROPIC (Bahl 2012) - Mitoxantrone xxxx
# 21             TROPIC (de Bono 2010) - Cabazitaxel
# 22            TROPIC (de Bono 2010) - Mitoxantrone

# Put items into a list
placebo.muhaz.plotting <- list()
placebo.muhaz.plotting[[02]] <- data.frame(Hazard = hazTemp[[02]]$haz.est, Time = hazTemp[[02]]$est.grid, Strata = "AFFIRM (Scher 2012) - Placebo")
placebo.muhaz.plotting[[03]] <- data.frame(Hazard = hazTemp[[03]]$haz.est, Time = hazTemp[[03]]$est.grid, Strata = "ALSYMPCA (Parker 2013) - Placebo")
placebo.muhaz.plotting[[06]] <- data.frame(Hazard = hazTemp[[06]]$haz.est, Time = hazTemp[[06]]$est.grid, Strata = "COU-AA-301 (de Bono 2011) - Placebo")
placebo.muhaz.plotting[[12]] <- data.frame(Hazard = hazTemp[[12]]$haz.est, Time = hazTemp[[12]]$est.grid, Strata = "PREVAIL (Beer 2017) - Placebo")
placebo.muhaz.plotting[[14]] <- data.frame(Hazard = hazTemp[[14]]$haz.est, Time = hazTemp[[14]]$est.grid, Strata = "TAX 327 (Berthold 2008) - Mitoxantrone")
placebo.muhaz.plotting[[20]] <- data.frame(Hazard = hazTemp[[20]]$haz.est, Time = hazTemp[[20]]$est.grid, Strata = "TROPIC (Bahl 2012) - Mitoxantrone")

# Combine into a single dataset
haztempplot.placebo <- rbind(placebo.muhaz.plotting[[02]],
                             placebo.muhaz.plotting[[03]],
                             placebo.muhaz.plotting[[06]],
                             placebo.muhaz.plotting[[12]],
                             placebo.muhaz.plotting[[14]],
                             placebo.muhaz.plotting[[20]])

# Ensure levels are ordered correctly
summary(haztempplot.placebo$Strata)
haztempplot.placebo$Strata <- factor(haztempplot.placebo$Strata)
summary(haztempplot.placebo$Strata)
haztempplot.placebo$Strata <- factor(haztempplot.placebo$Strata,levels(haztempplot.placebo$Strata)[c(1,2,3,4,5,6)])

# Produce list of relevant colours and line types for plotting
placebo_colours_list = c("AFFIRM (Scher 2012) - Placebo" = plot_red,
                         "ALSYMPCA (Parker 2013) - Placebo" = plot_bronze,
                         "COU-AA-301 (de Bono 2011) - Placebo" = "black",
                         "PREVAIL (Beer 2017) - Placebo" = plot_tropicalgreen,
                         "TAX 327 (Berthold 2008) - Mitoxantrone" = plot_lightblue,
                         "TROPIC (Bahl 2012) - Mitoxantrone" = plot_pink)

placebo_lines_list = c("AFFIRM (Scher 2012) - Placebo" = "solid",
                         "ALSYMPCA (Parker 2013) - Placebo" = "solid",
                         "COU-AA-301 (de Bono 2011) - Placebo" = "solid",
                         "PREVAIL (Beer 2017) - Placebo" = "solid",
                         "TAX 327 (Berthold 2008) - Mitoxantrone" = "solid",
                         "TROPIC (Bahl 2012) - Mitoxantrone" = "solid")

# Produce plot
muhaz_placebo <- ggplot(haztempplot.placebo, aes(x=Time, y=Hazard, colour= Strata, linetype = Strata)) +
  geom_line(size=1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,4), oob = rescale_none, breaks=c(0,1,2,3,4)) +
  scale_x_continuous(limits = c(0,7), oob = rescale_none, breaks=c(0,1,2,3,4,5,6,7)) +
  scale_colour_manual(values=placebo_colours_list)+
  scale_linetype_manual(values=placebo_lines_list)+
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8))

# Put items into a list
placebo.muhaz.plotting.rebased <- list()
placebo.muhaz.plotting.rebased[[02]] <- data.frame(Hazard = hazTemp.rebased[[02]]$haz.est, Time = hazTemp.rebased[[02]]$est.grid, Strata = "AFFIRM (Scher 2012) - Placebo")
placebo.muhaz.plotting.rebased[[03]] <- data.frame(Hazard = hazTemp.rebased[[03]]$haz.est, Time = hazTemp.rebased[[03]]$est.grid, Strata = "ALSYMPCA (Parker 2013) - Placebo")
placebo.muhaz.plotting.rebased[[06]] <- data.frame(Hazard = hazTemp.rebased[[06]]$haz.est, Time = hazTemp.rebased[[06]]$est.grid, Strata = "COU-AA-301 (de Bono 2011) - Placebo")
placebo.muhaz.plotting.rebased[[12]] <- data.frame(Hazard = hazTemp.rebased[[12]]$haz.est, Time = hazTemp.rebased[[12]]$est.grid, Strata = "PREVAIL (Beer 2017) - Placebo")
placebo.muhaz.plotting.rebased[[14]] <- data.frame(Hazard = hazTemp.rebased[[14]]$haz.est, Time = hazTemp.rebased[[14]]$est.grid, Strata = "TAX 327 (Berthold 2008) - Mitoxantrone")
placebo.muhaz.plotting.rebased[[20]] <- data.frame(Hazard = hazTemp.rebased[[20]]$haz.est, Time = hazTemp.rebased[[20]]$est.grid, Strata = "TROPIC (Bahl 2012) - Mitoxantrone")

# Combine into a single dataset
haztempplot.placebo.rebased <- rbind(placebo.muhaz.plotting.rebased[[02]],
                             placebo.muhaz.plotting.rebased[[03]],
                             placebo.muhaz.plotting.rebased[[06]],
                             placebo.muhaz.plotting.rebased[[12]],
                             placebo.muhaz.plotting.rebased[[14]],
                             placebo.muhaz.plotting.rebased[[20]])

# Ensure levels are ordered correctly
summary(haztempplot.placebo.rebased$Strata)
haztempplot.placebo.rebased$Strata <- factor(haztempplot.placebo.rebased$Strata)
summary(haztempplot.placebo.rebased$Strata)
haztempplot.placebo.rebased$Strata <- factor(haztempplot.placebo.rebased$Strata,levels(haztempplot.placebo.rebased$Strata)[c(1,2,3,4,5,6)])

# Produce plot
muhaz_placebo_rebased <- ggplot(haztempplot.placebo.rebased, aes(x=Time, y=Hazard, colour= Strata, linetype = Strata)) +
  geom_line(size=1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,4), oob = rescale_none, breaks=c(0,1,2,3,4)) +
  scale_x_continuous(limits = c(0,7), oob = rescale_none, breaks=c(0,1,2,3,4,5,6,7)) +
  scale_colour_manual(values=placebo_colours_list)+
  scale_linetype_manual(values=placebo_lines_list)+
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8))

# Combine into a single dataset
haztempplot.data.placebo <- rbind(survival.data.by.arm[[02]],
                                 survival.data.by.arm[[03]],
                                 survival.data.by.arm[[14]],
                                 survival.data.by.arm[[12]],
                                 survival.data.by.arm[[20]])

haztempplot.data.placebo.rebased <- rbind(survival.data.by.arm.rebased[[02]],
                                         survival.data.by.arm.rebased[[03]],
                                         survival.data.by.arm.rebased[[12]],
                                         survival.data.by.arm.rebased[[14]],
                                         survival.data.by.arm.rebased[[20]])

hazTemp.pooled.placebo = muhaz(time = haztempplot.data.placebo$Years, delta = haztempplot.data.placebo$Event, max.time = max(haztempplot.data.placebo$Years)*muhaz_mult)
hazTemp.pooled.placebo.plotting <- data.frame(Hazard = hazTemp.pooled.placebo$haz.est, Time = hazTemp.pooled.placebo$est.grid, Strata = "External data (pooled)")

hazTemp.pooled.placebo.rebased = muhaz(time = haztempplot.data.placebo.rebased$Years, delta = haztempplot.data.placebo.rebased$Event, max.time = max(haztempplot.data.placebo.rebased$Years)*muhaz_mult)
hazTemp.pooled.placebo.rebased.plotting <- data.frame(Hazard = hazTemp.pooled.placebo.rebased$haz.est, Time = hazTemp.pooled.placebo.rebased$est.grid, Strata = "External data (pooled, rebased)")

muhaz_placebo_pooled <- ggplot(hazTemp.pooled.placebo.plotting, aes(x=Time, y=Hazard)) +
  geom_line(size=1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,7), oob = rescale_none, breaks=c(0,1,2,3,4,5,6,7)) +
  scale_colour_manual(values=active_colours_list)+
  scale_linetype_manual(values=active_lines_list)+
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8))

muhaz_placebo_pooled_rebased <- ggplot(hazTemp.pooled.placebo.rebased.plotting, aes(x=Time, y=Hazard)) +
  geom_line(size=1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,7), oob = rescale_none, breaks=c(0,1,2,3,4,5,6,7)) +
  scale_colour_manual(values=active_colours_list)+
  scale_linetype_manual(values=active_lines_list)+
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8))

# Muhaz and pehaz plots for COU-AA-301 interim analysis # ----

aa_muhaz_active_original = muhaz(time = survival.data.by.arm[[5]]$Years, delta = survival.data.by.arm[[5]]$Event, max.time = max(survival.data.by.arm[[5]]$Years)*muhaz_mult)
aa_muhaz_placebo_original = muhaz(time = survival.data.by.arm[[6]]$Years, delta = survival.data.by.arm[[6]]$Event, max.time = max(survival.data.by.arm[[6]]$Years)*muhaz_mult)
aa_pehaz_active_original = pehaz(time = survival.data.by.arm[[5]]$Years, delta = survival.data.by.arm[[5]]$Event, max.time = max(survival.data.by.arm[[5]]$Years)*muhaz_mult)
aa_pehaz_placebo_original = pehaz(time = survival.data.by.arm[[6]]$Years, delta = survival.data.by.arm[[6]]$Event, max.time = max(survival.data.by.arm[[6]]$Years)*muhaz_mult)

pehaz_active_rows <- length(aa_pehaz_active_original$Cuts)*100
aa_pehaz_active_original_extract <- matrix(0, ncol=3, nrow = pehaz_active_rows)

for (i in 1:(pehaz_active_rows)){
  aa_pehaz_active_original_extract[i,1] <- (i-1)*(aa_pehaz_active_original$Width)/100
  aa_pehaz_active_original_extract[i,2] <- floor((i-1)/100)
}

for (i in 1:(pehaz_active_rows)){
  extr <- as.numeric(aa_pehaz_active_original_extract[i,2])
  aa_pehaz_active_original_extract[i,3] <- aa_pehaz_active_original$Hazard[extr+1]
}

pehaz_placebo_rows <- length(aa_pehaz_placebo_original$Cuts)*100
aa_pehaz_placebo_original_extract <- matrix(0, ncol=3, nrow = pehaz_placebo_rows)

for (i in 1:(pehaz_placebo_rows)){
  aa_pehaz_placebo_original_extract[i,1] <- (i-1)*(aa_pehaz_placebo_original$Width)/100
  aa_pehaz_placebo_original_extract[i,2] <- floor((i-1)/100)
}

for (i in 1:(pehaz_placebo_rows)){
  extr <- as.numeric(aa_pehaz_placebo_original_extract[i,2])
  aa_pehaz_placebo_original_extract[i,3] <- aa_pehaz_placebo_original$Hazard[extr+1]
}

aa_muhaz_active_original_times <- data.frame(Hazard = aa_muhaz_active_original$haz.est, Time = aa_muhaz_active_original$est.grid, Strata = "AA, smoothed")
aa_muhaz_placebo_original_times <- data.frame(Hazard = aa_muhaz_placebo_original$haz.est, Time = aa_muhaz_placebo_original$est.grid, Strata = "Placebo, smoothed")
aa_pehaz_active_original_times <- data.frame(Hazard = aa_pehaz_active_original_extract[,3], Time = aa_pehaz_active_original_extract[,1], Strata = "AA, piecewise")
aa_pehaz_placebo_original_times <- data.frame(Hazard = aa_pehaz_placebo_original_extract[,3], Time = aa_pehaz_placebo_original_extract[,1], Strata = "Placebo, piecewise")

muhaz_vs_pehaz <- rbind(aa_muhaz_active_original_times,
                        aa_muhaz_placebo_original_times,
                        aa_pehaz_active_original_times,
                        aa_pehaz_placebo_original_times)

# Produce list of relevant colours and line types for plotting
muhaz_pehaz_colours_list <- c("AA, smoothed" = plot_red,
                              "Placebo, smoothed" = plot_lightblue,
                              "AA, piecewise" = plot_red_2,
                              "Placebo, piecewise" = plot_lightblue_2)

muhaz_pehaz_lines_list <- c("AA, smoothed" = "dashed",
                            "Placebo, smoothed" = "dashed",
                            "AA, piecewise" = "solid",
                            "Placebo, piecewise" = "solid")
# Produce plot
muhaz_pehaz_plot <- ggplot(muhaz_vs_pehaz, aes(x=Time, y=Hazard, colour= Strata, linetype = Strata)) +
  geom_line(size=1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,5), oob = rescale_none, breaks=seq(0,5,0.5)) +
  scale_x_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_colour_manual(values=muhaz_pehaz_colours_list)+
  scale_linetype_manual(values=muhaz_pehaz_lines_list)+
  guides(colour=guide_legend(title="", nrow=1),linetype = guide_legend(title="", nrow=1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8))

# Combine datasets to allow for comparison of survival times # ----

# Combine active studies into a single dataset
# At this stage going forward, the PREVAIL study is dropped (which would be [[11]])
survival.data.external.active <- rbind(survival.data.by.arm[[01]],
                                       survival.data.by.arm[[04]],
                                       # survival.data.by.arm[[13]],
                                       # survival.data.by.arm[[15]],
                                       survival.data.by.arm[[19]])
# Separately, label pivotal trial sets
survival.data.pivotal.original.active <- survival.data.by.arm[[05]]
survival.data.pivotal.updated.active <- survival.data.by.arm[[07]]

# Combine placebo studies into a single dataset
# At this stage going forward, the PREVAIL study is dropped (which would be [[12]])
survival.data.external.placebo <- rbind(survival.data.by.arm[[02]],
                                        survival.data.by.arm[[03]],
                                        # survival.data.by.arm[[14]],
                                        survival.data.by.arm[[20]])

# Separately, label pivotal trial sets
survival.data.pivotal.original.placebo <- survival.data.by.arm[[06]]
survival.data.pivotal.updated.placebo <- survival.data.by.arm[[08]]

# Use an additional strata to separate possible sources
# Abi == 0 for external data, Abi == 1 for updated COU-AA-301 datacut, and Abi == 2 for original COU-AA-301 datacut

survival.data.external.active$Abi <- 0
survival.data.pivotal.updated.active$Abi <- 1
survival.data.pivotal.original.active$Abi <- 2

survival.data.external.placebo$Abi <- 0
survival.data.pivotal.updated.placebo$Abi <- 1
survival.data.pivotal.original.placebo$Abi <- 2

# Combine into one dataset

survival.data.active.all <- rbind(survival.data.external.active,
                                   survival.data.pivotal.updated.active,
                                   survival.data.pivotal.original.active)

survival.data.placebo.all <- rbind(survival.data.external.placebo,
                                   survival.data.pivotal.updated.placebo,
                                   survival.data.pivotal.original.placebo)

# Produce external dataset taking only time after a specified cut off # ----

# Perform on combined data set

survival.data.active.all.rebased <- survival.data.active.all
survival.data.placebo.all.rebased <- survival.data.placebo.all

survival.data.active.all.rebased <- subset(survival.data.active.all.rebased,survival.data.active.all.rebased$Years > piecewise_cutpoint + error_margin)
survival.data.active.all.rebased$Years <- survival.data.active.all.rebased$Years - piecewise_cutpoint

survival.data.placebo.all.rebased <- subset(survival.data.placebo.all.rebased,survival.data.placebo.all.rebased$Years > piecewise_cutpoint + error_margin)
survival.data.placebo.all.rebased$Years <- survival.data.placebo.all.rebased$Years - piecewise_cutpoint

# Repeat just for the external data

survival.data.external.active.rebased <- survival.data.external.active
survival.data.external.active.rebased <- subset(survival.data.external.active.rebased,survival.data.external.active.rebased$Years >= piecewise_cutpoint + error_margin)
survival.data.external.active.rebased$Years <- survival.data.external.active.rebased$Years - piecewise_cutpoint

survival.data.external.placebo.rebased <- survival.data.external.placebo
survival.data.external.placebo.rebased <- subset(survival.data.external.placebo.rebased,survival.data.external.placebo.rebased$Years >= piecewise_cutpoint + error_margin)
survival.data.external.placebo.rebased$Years <- survival.data.external.placebo.rebased$Years - piecewise_cutpoint

# Parametric survival models fitted to original COU-AA-301 dataset # ----

# Active arm
fit.expo.pivotal.original.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.active, dist="exp")
fit.weib.pivotal.original.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.active, dist="weibull")
fit.ggam.pivotal.original.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.active, dist="gengamma")
fit.gomp.pivotal.original.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.active, dist="gompertz")
fit.lnor.pivotal.original.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.active, dist="lnorm")
fit.llog.pivotal.original.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.active, dist="llogis")

# BIC = -2 * LL + log(N) * k
fit.expo.pivotal.original.active$BIC <- -2 * fit.expo.pivotal.original.active$loglik + log(fit.expo.pivotal.original.active$N)*fit.expo.pivotal.original.active$npars
fit.weib.pivotal.original.active$BIC <- -2 * fit.weib.pivotal.original.active$loglik + log(fit.weib.pivotal.original.active$N)*fit.weib.pivotal.original.active$npars
fit.ggam.pivotal.original.active$BIC <- -2 * fit.ggam.pivotal.original.active$loglik + log(fit.ggam.pivotal.original.active$N)*fit.ggam.pivotal.original.active$npars
fit.gomp.pivotal.original.active$BIC <- -2 * fit.gomp.pivotal.original.active$loglik + log(fit.gomp.pivotal.original.active$N)*fit.gomp.pivotal.original.active$npars
fit.lnor.pivotal.original.active$BIC <- -2 * fit.lnor.pivotal.original.active$loglik + log(fit.lnor.pivotal.original.active$N)*fit.lnor.pivotal.original.active$npars
fit.llog.pivotal.original.active$BIC <- -2 * fit.llog.pivotal.original.active$loglik + log(fit.llog.pivotal.original.active$N)*fit.llog.pivotal.original.active$npars

# Placebo arm
fit.expo.pivotal.original.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.placebo, dist="exp")
fit.weib.pivotal.original.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.placebo, dist="weibull")
fit.ggam.pivotal.original.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.placebo, dist="gengamma")
fit.gomp.pivotal.original.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.placebo, dist="gompertz")
fit.lnor.pivotal.original.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.placebo, dist="lnorm")
fit.llog.pivotal.original.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.original.placebo, dist="llogis")

# BIC = -2 * LL + log(N) * k
fit.expo.pivotal.original.placebo$BIC <- -2 * fit.expo.pivotal.original.placebo$loglik + log(fit.expo.pivotal.original.placebo$N)*fit.expo.pivotal.original.placebo$npars
fit.weib.pivotal.original.placebo$BIC <- -2 * fit.weib.pivotal.original.placebo$loglik + log(fit.weib.pivotal.original.placebo$N)*fit.weib.pivotal.original.placebo$npars
fit.ggam.pivotal.original.placebo$BIC <- -2 * fit.ggam.pivotal.original.placebo$loglik + log(fit.ggam.pivotal.original.placebo$N)*fit.ggam.pivotal.original.placebo$npars
fit.gomp.pivotal.original.placebo$BIC <- -2 * fit.gomp.pivotal.original.placebo$loglik + log(fit.gomp.pivotal.original.placebo$N)*fit.gomp.pivotal.original.placebo$npars
fit.lnor.pivotal.original.placebo$BIC <- -2 * fit.lnor.pivotal.original.placebo$loglik + log(fit.lnor.pivotal.original.placebo$N)*fit.lnor.pivotal.original.placebo$npars
fit.llog.pivotal.original.placebo$BIC <- -2 * fit.llog.pivotal.original.placebo$loglik + log(fit.llog.pivotal.original.placebo$N)*fit.llog.pivotal.original.placebo$npars

# Combine statistical fit scores

AICandBIC.pivotal.original <- cbind(rbind(fit.expo.pivotal.original.active$AIC,
                                          fit.weib.pivotal.original.active$AIC,
                                          fit.ggam.pivotal.original.active$AIC,
                                          fit.gomp.pivotal.original.active$AIC,
                                          fit.lnor.pivotal.original.active$AIC,
                                          fit.llog.pivotal.original.active$AIC),
                                    rbind(fit.expo.pivotal.original.active$BIC,
                                          fit.weib.pivotal.original.active$BIC,
                                          fit.ggam.pivotal.original.active$BIC,
                                          fit.gomp.pivotal.original.active$BIC,
                                          fit.lnor.pivotal.original.active$BIC,
                                          fit.llog.pivotal.original.active$BIC),
                                    rbind(fit.expo.pivotal.original.placebo$AIC,
                                          fit.weib.pivotal.original.placebo$AIC,
                                          fit.ggam.pivotal.original.placebo$AIC,
                                          fit.gomp.pivotal.original.placebo$AIC,
                                          fit.lnor.pivotal.original.placebo$AIC,
                                          fit.llog.pivotal.original.placebo$AIC),
                                    rbind(fit.expo.pivotal.original.placebo$BIC,
                                          fit.weib.pivotal.original.placebo$BIC,
                                          fit.ggam.pivotal.original.placebo$BIC,
                                          fit.gomp.pivotal.original.placebo$BIC,
                                          fit.lnor.pivotal.original.placebo$BIC,
                                          fit.llog.pivotal.original.placebo$BIC))

rownames(AICandBIC.pivotal.original) <- rbind("expo","weib","ggam","gomp","lnor","llog")
colnames(AICandBIC.pivotal.original) <- rbind("active, AIC","active, BIC","placebo, AIC","placebo, BIC")

# Parametric survival models fitted to updated COU-AA-301 dataset # ----

# Active arm
fit.expo.pivotal.updated.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.active, dist="exp")
fit.weib.pivotal.updated.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.active, dist="weibull")
fit.ggam.pivotal.updated.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.active, dist="gengamma")
fit.gomp.pivotal.updated.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.active, dist="gompertz")
fit.lnor.pivotal.updated.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.active, dist="lnorm")
fit.llog.pivotal.updated.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.active, dist="llogis")

# BIC = -2 * LL + log(N) * k
fit.expo.pivotal.updated.active$BIC <- -2 * fit.expo.pivotal.updated.active$loglik + log(fit.expo.pivotal.updated.active$N)*fit.expo.pivotal.updated.active$npars
fit.weib.pivotal.updated.active$BIC <- -2 * fit.weib.pivotal.updated.active$loglik + log(fit.weib.pivotal.updated.active$N)*fit.weib.pivotal.updated.active$npars
fit.ggam.pivotal.updated.active$BIC <- -2 * fit.ggam.pivotal.updated.active$loglik + log(fit.ggam.pivotal.updated.active$N)*fit.ggam.pivotal.updated.active$npars
fit.gomp.pivotal.updated.active$BIC <- -2 * fit.gomp.pivotal.updated.active$loglik + log(fit.gomp.pivotal.updated.active$N)*fit.gomp.pivotal.updated.active$npars
fit.lnor.pivotal.updated.active$BIC <- -2 * fit.lnor.pivotal.updated.active$loglik + log(fit.lnor.pivotal.updated.active$N)*fit.lnor.pivotal.updated.active$npars
fit.llog.pivotal.updated.active$BIC <- -2 * fit.llog.pivotal.updated.active$loglik + log(fit.llog.pivotal.updated.active$N)*fit.llog.pivotal.updated.active$npars

# Placebo arm
fit.expo.pivotal.updated.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.placebo, dist="exp")
fit.weib.pivotal.updated.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.placebo, dist="weibull")
fit.ggam.pivotal.updated.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.placebo, dist="gengamma")
fit.gomp.pivotal.updated.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.placebo, dist="gompertz")
fit.lnor.pivotal.updated.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.placebo, dist="lnorm")
fit.llog.pivotal.updated.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.pivotal.updated.placebo, dist="llogis")

# BIC = -2 * LL + log(N) * k
fit.expo.pivotal.updated.placebo$BIC <- -2 * fit.expo.pivotal.updated.placebo$loglik + log(fit.expo.pivotal.updated.placebo$N)*fit.expo.pivotal.updated.placebo$npars
fit.weib.pivotal.updated.placebo$BIC <- -2 * fit.weib.pivotal.updated.placebo$loglik + log(fit.weib.pivotal.updated.placebo$N)*fit.weib.pivotal.updated.placebo$npars
fit.ggam.pivotal.updated.placebo$BIC <- -2 * fit.ggam.pivotal.updated.placebo$loglik + log(fit.ggam.pivotal.updated.placebo$N)*fit.ggam.pivotal.updated.placebo$npars
fit.gomp.pivotal.updated.placebo$BIC <- -2 * fit.gomp.pivotal.updated.placebo$loglik + log(fit.gomp.pivotal.updated.placebo$N)*fit.gomp.pivotal.updated.placebo$npars
fit.lnor.pivotal.updated.placebo$BIC <- -2 * fit.lnor.pivotal.updated.placebo$loglik + log(fit.lnor.pivotal.updated.placebo$N)*fit.lnor.pivotal.updated.placebo$npars
fit.llog.pivotal.updated.placebo$BIC <- -2 * fit.llog.pivotal.updated.placebo$loglik + log(fit.llog.pivotal.updated.placebo$N)*fit.llog.pivotal.updated.placebo$npars

# Combine statistical fit scores

AICandBIC.pivotal.updated <- cbind(rbind(fit.expo.pivotal.updated.active$AIC,
                                          fit.weib.pivotal.updated.active$AIC,
                                          fit.ggam.pivotal.updated.active$AIC,
                                          fit.gomp.pivotal.updated.active$AIC,
                                          fit.lnor.pivotal.updated.active$AIC,
                                          fit.llog.pivotal.updated.active$AIC),
                                    rbind(fit.expo.pivotal.updated.active$BIC,
                                          fit.weib.pivotal.updated.active$BIC,
                                          fit.ggam.pivotal.updated.active$BIC,
                                          fit.gomp.pivotal.updated.active$BIC,
                                          fit.lnor.pivotal.updated.active$BIC,
                                          fit.llog.pivotal.updated.active$BIC),
                                    rbind(fit.expo.pivotal.updated.placebo$AIC,
                                          fit.weib.pivotal.updated.placebo$AIC,
                                          fit.ggam.pivotal.updated.placebo$AIC,
                                          fit.gomp.pivotal.updated.placebo$AIC,
                                          fit.lnor.pivotal.updated.placebo$AIC,
                                          fit.llog.pivotal.updated.placebo$AIC),
                                    rbind(fit.expo.pivotal.updated.placebo$BIC,
                                          fit.weib.pivotal.updated.placebo$BIC,
                                          fit.ggam.pivotal.updated.placebo$BIC,
                                          fit.gomp.pivotal.updated.placebo$BIC,
                                          fit.lnor.pivotal.updated.placebo$BIC,
                                          fit.llog.pivotal.updated.placebo$BIC))

rownames(AICandBIC.pivotal.updated) <- rbind("expo","weib","ggam","gomp","lnor","llog")
colnames(AICandBIC.pivotal.updated) <- rbind("active, AIC","active, BIC","placebo, AIC","placebo, BIC")

# Parametric survival models fitted to external (non-rebased) dataset # ----

# Active arm
fit.expo.external.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active, dist="exp")
fit.weib.external.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active, dist="weibull")
fit.ggam.external.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active, dist="gengamma")
fit.gomp.external.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active, dist="gompertz")
fit.lnor.external.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active, dist="lnorm")
fit.llog.external.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active, dist="llogis")

# BIC = -2 * LL + log(N) * k
fit.expo.external.active$BIC <- -2 * fit.expo.external.active$loglik + log(fit.expo.external.active$N)*fit.expo.external.active$npars
fit.weib.external.active$BIC <- -2 * fit.weib.external.active$loglik + log(fit.weib.external.active$N)*fit.weib.external.active$npars
fit.ggam.external.active$BIC <- -2 * fit.ggam.external.active$loglik + log(fit.ggam.external.active$N)*fit.ggam.external.active$npars
fit.gomp.external.active$BIC <- -2 * fit.gomp.external.active$loglik + log(fit.gomp.external.active$N)*fit.gomp.external.active$npars
fit.lnor.external.active$BIC <- -2 * fit.lnor.external.active$loglik + log(fit.lnor.external.active$N)*fit.lnor.external.active$npars
fit.llog.external.active$BIC <- -2 * fit.llog.external.active$loglik + log(fit.llog.external.active$N)*fit.llog.external.active$npars

# Placebo arm
fit.expo.external.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo, dist="exp")
fit.weib.external.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo, dist="weibull")
fit.ggam.external.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo, dist="gengamma")
fit.gomp.external.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo, dist="gompertz")
fit.lnor.external.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo, dist="lnorm")
fit.llog.external.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo, dist="llogis")

# BIC = -2 * LL + log(N) * k
fit.expo.external.placebo$BIC <- -2 * fit.expo.external.placebo$loglik + log(fit.expo.external.placebo$N)*fit.expo.external.placebo$npars
fit.weib.external.placebo$BIC <- -2 * fit.weib.external.placebo$loglik + log(fit.weib.external.placebo$N)*fit.weib.external.placebo$npars
fit.ggam.external.placebo$BIC <- -2 * fit.ggam.external.placebo$loglik + log(fit.ggam.external.placebo$N)*fit.ggam.external.placebo$npars
fit.gomp.external.placebo$BIC <- -2 * fit.gomp.external.placebo$loglik + log(fit.gomp.external.placebo$N)*fit.gomp.external.placebo$npars
fit.lnor.external.placebo$BIC <- -2 * fit.lnor.external.placebo$loglik + log(fit.lnor.external.placebo$N)*fit.lnor.external.placebo$npars
fit.llog.external.placebo$BIC <- -2 * fit.llog.external.placebo$loglik + log(fit.llog.external.placebo$N)*fit.llog.external.placebo$npars

# Combine statistical fit scores

AICandBIC.external <- cbind(rbind(fit.expo.external.active$AIC,
                                         fit.weib.external.active$AIC,
                                         fit.ggam.external.active$AIC,
                                         fit.gomp.external.active$AIC,
                                         fit.lnor.external.active$AIC,
                                         fit.llog.external.active$AIC),
                                   rbind(fit.expo.external.active$BIC,
                                         fit.weib.external.active$BIC,
                                         fit.ggam.external.active$BIC,
                                         fit.gomp.external.active$BIC,
                                         fit.lnor.external.active$BIC,
                                         fit.llog.external.active$BIC),
                                   rbind(fit.expo.external.placebo$AIC,
                                         fit.weib.external.placebo$AIC,
                                         fit.ggam.external.placebo$AIC,
                                         fit.gomp.external.placebo$AIC,
                                         fit.lnor.external.placebo$AIC,
                                         fit.llog.external.placebo$AIC),
                                   rbind(fit.expo.external.placebo$BIC,
                                         fit.weib.external.placebo$BIC,
                                         fit.ggam.external.placebo$BIC,
                                         fit.gomp.external.placebo$BIC,
                                         fit.lnor.external.placebo$BIC,
                                         fit.llog.external.placebo$BIC))

rownames(AICandBIC.external) <- rbind("expo","weib","ggam","gomp","lnor","llog")
colnames(AICandBIC.external) <- rbind("active, AIC","active, BIC","placebo, AIC","placebo, BIC")

# Parametric survival models fitted to external (rebased) dataset # ----

# Active arm
fit.expo.external.rebased.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active.rebased, dist="exp")
fit.weib.external.rebased.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active.rebased, dist="weibull")
fit.ggam.external.rebased.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active.rebased, dist="gengamma")
fit.gomp.external.rebased.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active.rebased, dist="gompertz")
fit.lnor.external.rebased.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active.rebased, dist="lnorm")
fit.llog.external.rebased.active <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.active.rebased, dist="llogis")
# BIC = -2 * LL + log(N) * k
fit.expo.external.rebased.active$BIC <- -2 * fit.expo.external.rebased.active$loglik + log(fit.expo.external.rebased.active$N)*fit.expo.external.rebased.active$npars
fit.weib.external.rebased.active$BIC <- -2 * fit.weib.external.rebased.active$loglik + log(fit.weib.external.rebased.active$N)*fit.weib.external.rebased.active$npars
fit.ggam.external.rebased.active$BIC <- -2 * fit.ggam.external.rebased.active$loglik + log(fit.ggam.external.rebased.active$N)*fit.ggam.external.rebased.active$npars
fit.gomp.external.rebased.active$BIC <- -2 * fit.gomp.external.rebased.active$loglik + log(fit.gomp.external.rebased.active$N)*fit.gomp.external.rebased.active$npars
fit.lnor.external.rebased.active$BIC <- -2 * fit.lnor.external.rebased.active$loglik + log(fit.lnor.external.rebased.active$N)*fit.lnor.external.rebased.active$npars
fit.llog.external.rebased.active$BIC <- -2 * fit.llog.external.rebased.active$loglik + log(fit.llog.external.rebased.active$N)*fit.llog.external.rebased.active$npars

# Placebo arm
fit.expo.external.rebased.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo.rebased, dist="exp")
fit.weib.external.rebased.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo.rebased, dist="weibull")
fit.ggam.external.rebased.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo.rebased, dist="gengamma")
fit.gomp.external.rebased.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo.rebased, dist="gompertz")
fit.lnor.external.rebased.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo.rebased, dist="lnorm")
fit.llog.external.rebased.placebo <- flexsurvreg(formula = Surv(Years, Event) ~ 1, data = survival.data.external.placebo.rebased, dist="llogis")

# BIC = -2 * LL + log(N) * k
fit.expo.external.rebased.placebo$BIC <- -2 * fit.expo.external.rebased.placebo$loglik + log(fit.expo.external.rebased.placebo$N)*fit.expo.external.rebased.placebo$npars
fit.weib.external.rebased.placebo$BIC <- -2 * fit.weib.external.rebased.placebo$loglik + log(fit.weib.external.rebased.placebo$N)*fit.weib.external.rebased.placebo$npars
fit.ggam.external.rebased.placebo$BIC <- -2 * fit.ggam.external.rebased.placebo$loglik + log(fit.ggam.external.rebased.placebo$N)*fit.ggam.external.rebased.placebo$npars
fit.gomp.external.rebased.placebo$BIC <- -2 * fit.gomp.external.rebased.placebo$loglik + log(fit.gomp.external.rebased.placebo$N)*fit.gomp.external.rebased.placebo$npars
fit.lnor.external.rebased.placebo$BIC <- -2 * fit.lnor.external.rebased.placebo$loglik + log(fit.lnor.external.rebased.placebo$N)*fit.lnor.external.rebased.placebo$npars
fit.llog.external.rebased.placebo$BIC <- -2 * fit.llog.external.rebased.placebo$loglik + log(fit.llog.external.rebased.placebo$N)*fit.llog.external.rebased.placebo$npars

# Combine statistical fit scores

AICandBIC.external.rebased <- cbind(rbind(fit.expo.external.rebased.active$AIC,
                                  fit.weib.external.rebased.active$AIC,
                                  fit.ggam.external.rebased.active$AIC,
                                  fit.gomp.external.rebased.active$AIC,
                                  fit.lnor.external.rebased.active$AIC,
                                  fit.llog.external.rebased.active$AIC),
                            rbind(fit.expo.external.rebased.active$BIC,
                                  fit.weib.external.rebased.active$BIC,
                                  fit.ggam.external.rebased.active$BIC,
                                  fit.gomp.external.rebased.active$BIC,
                                  fit.lnor.external.rebased.active$BIC,
                                  fit.llog.external.rebased.active$BIC),
                            rbind(fit.expo.external.rebased.placebo$AIC,
                                  fit.weib.external.rebased.placebo$AIC,
                                  fit.ggam.external.rebased.placebo$AIC,
                                  fit.gomp.external.rebased.placebo$AIC,
                                  fit.lnor.external.rebased.placebo$AIC,
                                  fit.llog.external.rebased.placebo$AIC),
                            rbind(fit.expo.external.rebased.placebo$BIC,
                                  fit.weib.external.rebased.placebo$BIC,
                                  fit.ggam.external.rebased.placebo$BIC,
                                  fit.gomp.external.rebased.placebo$BIC,
                                  fit.lnor.external.rebased.placebo$BIC,
                                  fit.llog.external.rebased.placebo$BIC))

rownames(AICandBIC.external.rebased) <- rbind("expo","weib","ggam","gomp","lnor","llog")
colnames(AICandBIC.external.rebased) <- rbind("active, AIC","active, BIC","placebo, AIC","placebo, BIC")


# Plot KM curves by source (external, original, or updated) # ----

survival.data.external.active$Abi <- "External data (pooled)"
survival.data.external.active.rebased$Abi <- "External data (pooled)"
survival.data.external.placebo$Abi <- "External data (pooled)"
survival.data.external.placebo.rebased$Abi <- "External data (pooled)"

KM.est.active.pivotal.original<- survfit(Surv(Years,Event) ~ Abi, data = survival.data.pivotal.original.active, type = "kaplan-meier")
KM.est.active.pivotal.updated<- survfit(Surv(Years,Event) ~ Abi, data = survival.data.pivotal.updated.active, type = "kaplan-meier")
KM.est.active.external<- survfit(Surv(Years,Event) ~ Abi, data = survival.data.external.active, type = "kaplan-meier")
KM.est.active.external.rebased<- survfit(Surv(Years,Event) ~ Abi, data = survival.data.external.active.rebased, type = "kaplan-meier")

KM.est.placebo.pivotal.original<- survfit(Surv(Years,Event) ~ Abi, data = survival.data.pivotal.original.placebo, type = "kaplan-meier")
KM.est.placebo.pivotal.updated<- survfit(Surv(Years,Event) ~ Abi, data = survival.data.pivotal.updated.placebo, type = "kaplan-meier")
KM.est.placebo.external<- survfit(Surv(Years,Event) ~ Abi, data = survival.data.external.placebo, type = "kaplan-meier")
KM.est.placebo.external.rebased<- survfit(Surv(Years,Event) ~ Abi, data = survival.data.external.placebo.rebased, type = "kaplan-meier")

KM.list.active <-list(Original = KM.est.active.pivotal.original,
                      Updated = KM.est.active.pivotal.updated,
                      External = KM.est.active.external)

KM_active_compare <- ggsurvplot_combine(KM.list.active,
                                        survival.data.active.all,
                                        palette=c("black","grey40","grey80"),
                                        break.x.by = 1,
                                        xlim=c(0,10),
                                        linetype=c(1,2,3),
                                        surv.scale="percent",
                                        ylab="Overall survival",
                                        xlab="Time (years)",
                                        fontsize = 3.5,
                                        font.x= c(10, "bold", "black"),
                                        font.y= c(10, "bold", "black"),
                                        censor.size=2,
                                        break.y.by=0.1,
                                        font.tickslab= c(10, "plain", "black"),
                                        legend="bottom",
                                        legend.labs=c("COU-AA-301: Original","COU-AA-301: Updated","External data (pooled)"),
                                        legend.title="",
                                        conf.int = TRUE
) + guides(fill = guide_legend(title="", nrow=4))

KM.list.placebo <-list(Original = KM.est.placebo.pivotal.original,
                       Updated = KM.est.placebo.pivotal.updated,
                       External = KM.est.placebo.external)

KM_placebo_compare <- ggsurvplot_combine(KM.list.placebo,
                                         survival.data.placebo.all,
                                         palette=c("black","grey40","grey80"),
                                         break.x.by = 1,
                                         xlim=c(0,10),
                                         linetype=c(1,2,3),
                                         surv.scale="percent",
                                         ylab="Overall survival",
                                         xlab="Time (years)",
                                         fontsize = 3.5,
                                         font.x= c(10, "bold", "black"),
                                         font.y= c(10, "bold", "black"),
                                         censor.size=2,
                                         break.y.by=0.1,
                                         font.tickslab= c(10, "plain", "black"),
                                         legend="bottom",
                                         legend.labs=c("COU-AA-301: Original","COU-AA-301: Updated","External data (pooled)"),
                                         legend.title="",
                                         conf.int = TRUE
) + guides(fill = guide_legend(title="", nrow=4))

KM.list.active.newold <-list(Original = KM.est.active.pivotal.original,
                      Updated = KM.est.active.pivotal.updated)

KM_active_newold_ci <- ggsurvplot_combine(KM.list.active.newold,
                                        survival.data.active.all,
                                        palette=c("black","grey40"),
                                        break.x.by = 0.5,
                                        xlim=c(0,3),
                                        linetype=c(1,2),
                                        surv.scale="percent",
                                        ylab="Overall survival",
                                        xlab="Time (years)",
                                        fontsize = 3.5,
                                        font.x= c(10, "bold", "black"),
                                        font.y= c(10, "bold", "black"),
                                        censor.size=2,
                                        break.y.by=0.1,
                                        font.tickslab= c(10, "plain", "black"),
                                        legend="bottom",
                                        legend.labs=c("COU-AA-301: Original","COU-AA-301: Updated"),
                                        legend.title="",
                                        conf.int = TRUE
)

KM.list.placebo.newold <-list(Original = KM.est.placebo.pivotal.original,
                       Updated = KM.est.placebo.pivotal.updated)

KM_placebo_newold_ci <- ggsurvplot_combine(KM.list.placebo.newold,
                                         survival.data.placebo.all,
                                         palette=c("black","grey40"),
                                         break.x.by = 0.5,
                                         xlim=c(0,3),
                                         linetype=c(1,2),
                                         surv.scale="percent",
                                         ylab="Overall survival",
                                         xlab="Time (years)",
                                         fontsize = 3.5,
                                         font.x= c(10, "bold", "black"),
                                         font.y= c(10, "bold", "black"),
                                         censor.size=2,
                                         break.y.by=0.1,
                                         font.tickslab= c(10, "plain", "black"),
                                         legend="bottom",
                                         legend.labs=c("COU-AA-301: Original","COU-AA-301: Updated"),
                                         legend.title="",
                                         conf.int = TRUE
)

KM.list.active.oldexternal <-list(Original = KM.est.active.pivotal.original,
                                   External = KM.est.active.external)

KM_active_oldexternal_ci <- ggsurvplot_combine(KM.list.active.oldexternal,
                                                survival.data.active.all,
                                                palette=c("black","grey60"),
                                                break.x.by = 0.5,
                                                xlim=c(0,4),
                                                linetype=c(1,2),
                                                fill = FALSE,
                                                surv.scale="percent",
                                                ylab="Overall survival",
                                                xlab="Time (years)",
                                                fontsize = 3.5,
                                                font.x= c(10, "bold", "black"),
                                                font.y= c(10, "bold", "black"),
                                                censor.size=2,
                                                break.y.by=0.1,
                                                font.tickslab= c(10, "plain", "black"),
                                                legend="bottom",
                                                legend.labs=c("COU-AA-301: Original","External data (pooled)"),
                                                legend.title="",
                                                conf.int = TRUE
)

KM_active_oldexternal_ci_10 <- ggsurvplot_combine(KM.list.active.oldexternal,
                                               survival.data.active.all,
                                               palette=c("black","grey60"),
                                               break.x.by = 1,
                                               xlim=c(0,10),
                                               linetype=c(1,2),
                                               fill = FALSE,
                                               surv.scale="percent",
                                               ylab="Overall survival",
                                               xlab="Time (years)",
                                               fontsize = 3.5,
                                               font.x= c(10, "bold", "black"),
                                               font.y= c(10, "bold", "black"),
                                               censor.size=2,
                                               break.y.by=0.1,
                                               font.tickslab= c(10, "plain", "black"),
                                               legend="bottom",
                                               legend.labs=c("COU-AA-301: Original","External data (pooled)"),
                                               legend.title="",
                                               conf.int = TRUE
)

KM.list.placebo.oldexternal <-list(Original = KM.est.placebo.pivotal.original,
                              External = KM.est.placebo.external)

KM_placebo_oldexternal_ci <- ggsurvplot_combine(KM.list.placebo.oldexternal,
                                           survival.data.placebo.all,
                                           palette=c("black","grey60"),
                                           break.x.by = 0.5,
                                           xlim=c(0,4),
                                           linetype=c(1,2),
                                           surv.scale="percent",
                                           ylab="Overall survival",
                                           xlab="Time (years)",
                                           fontsize = 3.5,
                                           font.x= c(10, "bold", "black"),
                                           font.y= c(10, "bold", "black"),
                                           censor.size=2,
                                           break.y.by=0.1,
                                           font.tickslab= c(10, "plain", "black"),
                                           legend="bottom",
                                           legend.labs=c("COU-AA-301: Original","External data (pooled)"),
                                           legend.title="",
                                           conf.int = TRUE
)

KM_placebo_oldexternal_ci_10 <- ggsurvplot_combine(KM.list.placebo.oldexternal,
                                                survival.data.placebo.all,
                                                palette=c("black","grey60"),
                                                break.x.by = 1,
                                                xlim=c(0,10),
                                                linetype=c(1,2),
                                                surv.scale="percent",
                                                ylab="Overall survival",
                                                xlab="Time (years)",
                                                fontsize = 3.5,
                                                font.x= c(10, "bold", "black"),
                                                font.y= c(10, "bold", "black"),
                                                censor.size=2,
                                                break.y.by=0.1,
                                                font.tickslab= c(10, "plain", "black"),
                                                legend="bottom",
                                                legend.labs=c("COU-AA-301: Original","External data (pooled)"),
                                                legend.title="",
                                                conf.int = TRUE
)

curves_colours_list_external <- c("External data (pooled)" = "grey80",
                                  "Extrapolation: Weibull" = plot_pink,
                                  "Extrapolation: Gompertz" = plot_forestgreen,
                                  "Extrapolation: Generalised gamma" = plot_bronze,
                                  "Extrapolation: Exponential" = plot_red,
                                  "Extrapolation: Log-logistic" = plot_lightblue,
                                  "Extrapolation: Lognormal" = plot_purple)

curves_lines_list_external <- c("External data (pooled)" = "solid",
                                "Extrapolation: Weibull" = "dashed",
                                "Extrapolation: Gompertz" = "dashed",
                                "Extrapolation: Generalised gamma" = "dashed",
                                "Extrapolation: Exponential" = "dashed",
                                "Extrapolation: Log-logistic" = "dashed",
                                "Extrapolation: Lognormal" = "dashed")

curves_fill_list_external <- c("External data (pooled)" = "grey80",
                               "Extrapolation: Weibull" = "white",
                               "Extrapolation: Gompertz" = "white",
                               "Extrapolation: Generalised gamma" = "white",
                               "Extrapolation: Exponential" = "white",
                               "Extrapolation: Log-logistic" = "white",
                               "Extrapolation: Lognormal" = "white")

KM_active_external <- ggsurvplot(KM.est.active.external,
                                 survival.data.external.active,
                                 palette=c("grey80"),
                                 break.x.by = 1,
                                 xlim=c(0,10),
                                 surv.scale="percent",
                                 ylab="Overall survival",
                                 xlab="Time (years)",
                                 fontsize = 3.5,
                                 font.x= c(10, "bold", "black"),
                                 font.y= c(10, "bold", "black"),
                                 censor.size=2,
                                 break.y.by=0.1,
                                 linetype=curves_lines_list_external,
                                 font.tickslab= c(10, "plain", "black"),
                                 legend="bottom",
                                 legend.labs=c("External data (pooled)"),
                                 legend.title=""
)

KM_placebo_external <- ggsurvplot(KM.est.placebo.external,
                                  survival.data.external.placebo,
                                  palette=c("grey80"),
                                  break.x.by = 1,
                                  xlim=c(0,10),
                                  surv.scale="percent",
                                  ylab="Overall survival",
                                  xlab="Time (years)",
                                  fontsize = 3.5,
                                  font.x= c(10, "bold", "black"),
                                  font.y= c(10, "bold", "black"),
                                  censor.size=2,
                                  break.y.by=0.1,
                                  linetype=curves_lines_list_external,
                                  font.tickslab= c(10, "plain", "black"),
                                  legend="bottom",
                                  legend.labs=c("External data (pooled)"),
                                  legend.title=""
)

KM_active_external_rebased <- ggsurvplot(KM.est.active.external.rebased,
                                 survival.data.external.active.rebased,
                                 palette=c("grey80"),
                                 break.x.by = 1,
                                 xlim=c(0,10),
                                 surv.scale="percent",
                                 ylab="Overall survival",
                                 xlab="Time (years)",
                                 fontsize = 3.5,
                                 font.x= c(10, "bold", "black"),
                                 font.y= c(10, "bold", "black"),
                                 censor.size=2,
                                 break.y.by=0.1,
                                 linetype=curves_lines_list_external,
                                 font.tickslab= c(10, "plain", "black"),
                                 legend="bottom",
                                 legend.labs=c("External data (pooled)"),
                                 legend.title=""
)

KM_placebo_external_rebased <- ggsurvplot(KM.est.placebo.external.rebased,
                                  survival.data.external.placebo.rebased,
                                  palette=c("grey80"),
                                  break.x.by = 1,
                                  xlim=c(0,10),
                                  surv.scale="percent",
                                  ylab="Overall survival",
                                  xlab="Time (years)",
                                  fontsize = 3.5,
                                  font.x= c(10, "bold", "black"),
                                  font.y= c(10, "bold", "black"),
                                  censor.size=2,
                                  break.y.by=0.1,
                                  linetype=curves_lines_list_external,
                                  font.tickslab= c(10, "plain", "black"),
                                  legend="bottom",
                                  legend.labs=c("External data (pooled)"),
                                  legend.title=""
)

# Extract parametric model information for plotting # ----

lower_limit <- 0
upper_limit <- 20
increment <- 0.01

times.expo.pivotal.original.active <- summary(fit.expo.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.weib.pivotal.original.active <- summary(fit.weib.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.ggam.pivotal.original.active <- summary(fit.ggam.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.gomp.pivotal.original.active <- summary(fit.gomp.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.lnor.pivotal.original.active <- summary(fit.lnor.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.llog.pivotal.original.active <- summary(fit.llog.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)

times.expo.pivotal.original.active$Abi <- "Extrapolation: Exponential"       
times.weib.pivotal.original.active$Abi <- "Extrapolation: Weibull"          
times.ggam.pivotal.original.active$Abi <- "Extrapolation: Generalised gamma"
times.gomp.pivotal.original.active$Abi <- "Extrapolation: Gompertz"              
times.lnor.pivotal.original.active$Abi <- "Extrapolation: Lognormal"             
times.llog.pivotal.original.active$Abi <- "Extrapolation: Log-logistic"         

times.expo.pivotal.original.placebo <- summary(fit.expo.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.weib.pivotal.original.placebo <- summary(fit.weib.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.ggam.pivotal.original.placebo <- summary(fit.ggam.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.gomp.pivotal.original.placebo <- summary(fit.gomp.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.lnor.pivotal.original.placebo <- summary(fit.lnor.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.llog.pivotal.original.placebo <- summary(fit.llog.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)

times.expo.pivotal.original.placebo$Abi <- "Extrapolation: Exponential"      
times.weib.pivotal.original.placebo$Abi <- "Extrapolation: Weibull"          
times.ggam.pivotal.original.placebo$Abi <- "Extrapolation: Generalised gamma"
times.gomp.pivotal.original.placebo$Abi <- "Extrapolation: Gompertz"         
times.lnor.pivotal.original.placebo$Abi <- "Extrapolation: Lognormal"        
times.llog.pivotal.original.placebo$Abi <- "Extrapolation: Log-logistic"     

times.expo.pivotal.updated.active <- summary(fit.expo.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.weib.pivotal.updated.active <- summary(fit.weib.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.ggam.pivotal.updated.active <- summary(fit.ggam.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.gomp.pivotal.updated.active <- summary(fit.gomp.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.lnor.pivotal.updated.active <- summary(fit.lnor.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.llog.pivotal.updated.active <- summary(fit.llog.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)

times.expo.pivotal.updated.active$Abi <- "Extrapolation: Exponential"      
times.weib.pivotal.updated.active$Abi <- "Extrapolation: Weibull"          
times.ggam.pivotal.updated.active$Abi <- "Extrapolation: Generalised gamma"
times.gomp.pivotal.updated.active$Abi <- "Extrapolation: Gompertz"         
times.lnor.pivotal.updated.active$Abi <- "Extrapolation: Lognormal"        
times.llog.pivotal.updated.active$Abi <- "Extrapolation: Log-logistic"     

times.expo.pivotal.updated.placebo <- summary(fit.expo.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.weib.pivotal.updated.placebo <- summary(fit.weib.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.ggam.pivotal.updated.placebo <- summary(fit.ggam.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.gomp.pivotal.updated.placebo <- summary(fit.gomp.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.lnor.pivotal.updated.placebo <- summary(fit.lnor.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.llog.pivotal.updated.placebo <- summary(fit.llog.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)

times.expo.pivotal.updated.placebo$Abi <- "Extrapolation: Exponential"      
times.weib.pivotal.updated.placebo$Abi <- "Extrapolation: Weibull"          
times.ggam.pivotal.updated.placebo$Abi <- "Extrapolation: Generalised gamma"
times.gomp.pivotal.updated.placebo$Abi <- "Extrapolation: Gompertz"         
times.lnor.pivotal.updated.placebo$Abi <- "Extrapolation: Lognormal"        
times.llog.pivotal.updated.placebo$Abi <- "Extrapolation: Log-logistic"     

times.expo.external.active <- summary(fit.expo.external.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.weib.external.active <- summary(fit.weib.external.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.ggam.external.active <- summary(fit.ggam.external.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.gomp.external.active <- summary(fit.gomp.external.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.lnor.external.active <- summary(fit.lnor.external.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.llog.external.active <- summary(fit.llog.external.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)

times.expo.external.active$Abi <- "Extrapolation: Exponential"      
times.weib.external.active$Abi <- "Extrapolation: Weibull"          
times.ggam.external.active$Abi <- "Extrapolation: Generalised gamma"
times.gomp.external.active$Abi <- "Extrapolation: Gompertz"         
times.lnor.external.active$Abi <- "Extrapolation: Lognormal"        
times.llog.external.active$Abi <- "Extrapolation: Log-logistic"     

times.expo.external.placebo <- summary(fit.expo.external.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.weib.external.placebo <- summary(fit.weib.external.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.ggam.external.placebo <- summary(fit.ggam.external.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.gomp.external.placebo <- summary(fit.gomp.external.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.lnor.external.placebo <- summary(fit.lnor.external.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.llog.external.placebo <- summary(fit.llog.external.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)

times.expo.external.placebo$Abi <- "Extrapolation: Exponential"      
times.weib.external.placebo$Abi <- "Extrapolation: Weibull"          
times.ggam.external.placebo$Abi <- "Extrapolation: Generalised gamma"
times.gomp.external.placebo$Abi <- "Extrapolation: Gompertz"         
times.lnor.external.placebo$Abi <- "Extrapolation: Lognormal"        
times.llog.external.placebo$Abi <- "Extrapolation: Log-logistic"     

times.expo.external.rebased.active <- summary(fit.expo.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.weib.external.rebased.active <- summary(fit.weib.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.ggam.external.rebased.active <- summary(fit.ggam.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.gomp.external.rebased.active <- summary(fit.gomp.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.lnor.external.rebased.active <- summary(fit.lnor.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.llog.external.rebased.active <- summary(fit.llog.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)

times.expo.external.rebased.active$Abi <- "Extrapolation: Exponential"      
times.weib.external.rebased.active$Abi <- "Extrapolation: Weibull"          
times.ggam.external.rebased.active$Abi <- "Extrapolation: Generalised gamma"
times.gomp.external.rebased.active$Abi <- "Extrapolation: Gompertz"         
times.lnor.external.rebased.active$Abi <- "Extrapolation: Lognormal"        
times.llog.external.rebased.active$Abi <- "Extrapolation: Log-logistic"     

times.expo.external.rebased.placebo <- summary(fit.expo.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.weib.external.rebased.placebo <- summary(fit.weib.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.ggam.external.rebased.placebo <- summary(fit.ggam.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.gomp.external.rebased.placebo <- summary(fit.gomp.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.lnor.external.rebased.placebo <- summary(fit.lnor.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)
times.llog.external.rebased.placebo <- summary(fit.llog.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="survival",ci=FALSE,tidy=TRUE)

times.expo.external.rebased.placebo$Abi <- "Extrapolation: Exponential"      
times.weib.external.rebased.placebo$Abi <- "Extrapolation: Weibull"          
times.ggam.external.rebased.placebo$Abi <- "Extrapolation: Generalised gamma"
times.gomp.external.rebased.placebo$Abi <- "Extrapolation: Gompertz"         
times.lnor.external.rebased.placebo$Abi <- "Extrapolation: Lognormal"        
times.llog.external.rebased.placebo$Abi <- "Extrapolation: Log-logistic"     

hazards.expo.pivotal.original.active <- summary(fit.expo.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.weib.pivotal.original.active <- summary(fit.weib.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.ggam.pivotal.original.active <- summary(fit.ggam.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.gomp.pivotal.original.active <- summary(fit.gomp.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.lnor.pivotal.original.active <- summary(fit.lnor.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.llog.pivotal.original.active <- summary(fit.llog.pivotal.original.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)

hazards.expo.pivotal.original.active$Abi <- "Extrapolation: Exponential"       
hazards.weib.pivotal.original.active$Abi <- "Extrapolation: Weibull"          
hazards.ggam.pivotal.original.active$Abi <- "Extrapolation: Generalised gamma"
hazards.gomp.pivotal.original.active$Abi <- "Extrapolation: Gompertz"              
hazards.lnor.pivotal.original.active$Abi <- "Extrapolation: Lognormal"             
hazards.llog.pivotal.original.active$Abi <- "Extrapolation: Log-logistic"         

hazards.expo.pivotal.original.placebo <- summary(fit.expo.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.weib.pivotal.original.placebo <- summary(fit.weib.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.ggam.pivotal.original.placebo <- summary(fit.ggam.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.gomp.pivotal.original.placebo <- summary(fit.gomp.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.lnor.pivotal.original.placebo <- summary(fit.lnor.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.llog.pivotal.original.placebo <- summary(fit.llog.pivotal.original.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)

hazards.expo.pivotal.original.placebo$Abi <- "Extrapolation: Exponential"      
hazards.weib.pivotal.original.placebo$Abi <- "Extrapolation: Weibull"          
hazards.ggam.pivotal.original.placebo$Abi <- "Extrapolation: Generalised gamma"
hazards.gomp.pivotal.original.placebo$Abi <- "Extrapolation: Gompertz"         
hazards.lnor.pivotal.original.placebo$Abi <- "Extrapolation: Lognormal"        
hazards.llog.pivotal.original.placebo$Abi <- "Extrapolation: Log-logistic"     

hazards.expo.pivotal.updated.active <- summary(fit.expo.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.weib.pivotal.updated.active <- summary(fit.weib.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.ggam.pivotal.updated.active <- summary(fit.ggam.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.gomp.pivotal.updated.active <- summary(fit.gomp.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.lnor.pivotal.updated.active <- summary(fit.lnor.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.llog.pivotal.updated.active <- summary(fit.llog.pivotal.updated.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)

hazards.expo.pivotal.updated.active$Abi <- "Extrapolation: Exponential"      
hazards.weib.pivotal.updated.active$Abi <- "Extrapolation: Weibull"          
hazards.ggam.pivotal.updated.active$Abi <- "Extrapolation: Generalised gamma"
hazards.gomp.pivotal.updated.active$Abi <- "Extrapolation: Gompertz"         
hazards.lnor.pivotal.updated.active$Abi <- "Extrapolation: Lognormal"        
hazards.llog.pivotal.updated.active$Abi <- "Extrapolation: Log-logistic"     

hazards.expo.pivotal.updated.placebo <- summary(fit.expo.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.weib.pivotal.updated.placebo <- summary(fit.weib.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.ggam.pivotal.updated.placebo <- summary(fit.ggam.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.gomp.pivotal.updated.placebo <- summary(fit.gomp.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.lnor.pivotal.updated.placebo <- summary(fit.lnor.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.llog.pivotal.updated.placebo <- summary(fit.llog.pivotal.updated.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)

hazards.expo.pivotal.updated.placebo$Abi <- "Extrapolation: Exponential"      
hazards.weib.pivotal.updated.placebo$Abi <- "Extrapolation: Weibull"          
hazards.ggam.pivotal.updated.placebo$Abi <- "Extrapolation: Generalised gamma"
hazards.gomp.pivotal.updated.placebo$Abi <- "Extrapolation: Gompertz"         
hazards.lnor.pivotal.updated.placebo$Abi <- "Extrapolation: Lognormal"        
hazards.llog.pivotal.updated.placebo$Abi <- "Extrapolation: Log-logistic"     

hazards.expo.external.active <- summary(fit.expo.external.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.weib.external.active <- summary(fit.weib.external.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.ggam.external.active <- summary(fit.ggam.external.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.gomp.external.active <- summary(fit.gomp.external.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.lnor.external.active <- summary(fit.lnor.external.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.llog.external.active <- summary(fit.llog.external.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)

hazards.expo.external.active$Abi <- "Extrapolation: Exponential"      
hazards.weib.external.active$Abi <- "Extrapolation: Weibull"          
hazards.ggam.external.active$Abi <- "Extrapolation: Generalised gamma"
hazards.gomp.external.active$Abi <- "Extrapolation: Gompertz"         
hazards.lnor.external.active$Abi <- "Extrapolation: Lognormal"        
hazards.llog.external.active$Abi <- "Extrapolation: Log-logistic"     

hazards.expo.external.placebo <- summary(fit.expo.external.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.weib.external.placebo <- summary(fit.weib.external.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.ggam.external.placebo <- summary(fit.ggam.external.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.gomp.external.placebo <- summary(fit.gomp.external.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.lnor.external.placebo <- summary(fit.lnor.external.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.llog.external.placebo <- summary(fit.llog.external.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)

hazards.expo.external.placebo$Abi <- "Extrapolation: Exponential"      
hazards.weib.external.placebo$Abi <- "Extrapolation: Weibull"          
hazards.ggam.external.placebo$Abi <- "Extrapolation: Generalised gamma"
hazards.gomp.external.placebo$Abi <- "Extrapolation: Gompertz"         
hazards.lnor.external.placebo$Abi <- "Extrapolation: Lognormal"        
hazards.llog.external.placebo$Abi <- "Extrapolation: Log-logistic"     

hazards.expo.external.rebased.active <- summary(fit.expo.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.weib.external.rebased.active <- summary(fit.weib.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.ggam.external.rebased.active <- summary(fit.ggam.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.gomp.external.rebased.active <- summary(fit.gomp.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.lnor.external.rebased.active <- summary(fit.lnor.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.llog.external.rebased.active <- summary(fit.llog.external.rebased.active, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)

hazards.expo.external.rebased.active$Abi <- "Extrapolation: Exponential"      
hazards.weib.external.rebased.active$Abi <- "Extrapolation: Weibull"          
hazards.ggam.external.rebased.active$Abi <- "Extrapolation: Generalised gamma"
hazards.gomp.external.rebased.active$Abi <- "Extrapolation: Gompertz"         
hazards.lnor.external.rebased.active$Abi <- "Extrapolation: Lognormal"        
hazards.llog.external.rebased.active$Abi <- "Extrapolation: Log-logistic"     

hazards.expo.external.rebased.placebo <- summary(fit.expo.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.weib.external.rebased.placebo <- summary(fit.weib.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.ggam.external.rebased.placebo <- summary(fit.ggam.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.gomp.external.rebased.placebo <- summary(fit.gomp.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.lnor.external.rebased.placebo <- summary(fit.lnor.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)
hazards.llog.external.rebased.placebo <- summary(fit.llog.external.rebased.placebo, t=seq(lower_limit, upper_limit, increment),type="hazard",ci=FALSE,tidy=TRUE)

hazards.expo.external.rebased.placebo$Abi <- "Extrapolation: Exponential"      
hazards.weib.external.rebased.placebo$Abi <- "Extrapolation: Weibull"          
hazards.ggam.external.rebased.placebo$Abi <- "Extrapolation: Generalised gamma"
hazards.gomp.external.rebased.placebo$Abi <- "Extrapolation: Gompertz"         
hazards.lnor.external.rebased.placebo$Abi <- "Extrapolation: Lognormal"        
hazards.llog.external.rebased.placebo$Abi <- "Extrapolation: Log-logistic"   

# Plot basic models fitted to original COU-AA-301 data on top of KM curves # ----

curves_colours_list <- c("COU-AA-301: Original" = "black",
                                "COU-AA-301: Updated" = "grey40",
                                "External data (pooled)" = "grey80",
                                "Extrapolation: Weibull" = plot_pink,
                                "Extrapolation: Gompertz" = plot_forestgreen,
                                "Extrapolation: Generalised gamma" = plot_bronze,
                                "Extrapolation: Exponential" = plot_red,
                                "Extrapolation: Log-logistic" = plot_lightblue,
                                "Extrapolation: Lognormal" = plot_purple)

curves_lines_list <- c("COU-AA-301: Original" = "solid",
                              "COU-AA-301: Updated" = "solid",
                              "External data (pooled)" = "solid",
                              "Extrapolation: Weibull" = "dashed",
                              "Extrapolation: Gompertz" = "dashed",
                              "Extrapolation: Generalised gamma" = "dashed",
                              "Extrapolation: Exponential" = "dashed",
                              "Extrapolation: Log-logistic" = "dashed",
                              "Extrapolation: Lognormal" = "dashed")

curves_fill_list <- c("COU-AA-301: Original" = "black",
                       "COU-AA-301: Updated" = "grey40",
                       "External data (pooled)" = "grey80",
                      "Extrapolation: Weibull" = "white",
                      "Extrapolation: Gompertz" = "white",
                      "Extrapolation: Generalised gamma" = "white",
                      "Extrapolation: Exponential" = "white",
                      "Extrapolation: Log-logistic" = "white",
                      "Extrapolation: Lognormal" = "white")

KM_active_compare_basicmodels <- KM_active_compare$plot +
  geom_line(data = times.expo.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) +
  scale_linetype_manual(values=curves_lines_list) +
  scale_colour_manual(values=curves_colours_list) +
  scale_fill_manual(values=curves_fill_list)

curves_colours_list_v2 <- c("COU-AA-301: Original" = "black",
                         "External data (pooled)" = "grey60",
                         "Extrapolation: Weibull" = plot_pink,
                         "Extrapolation: Gompertz" = plot_forestgreen,
                         "Extrapolation: Generalised gamma" = plot_bronze,
                         "Extrapolation: Exponential" = plot_red,
                         "Extrapolation: Log-logistic" = plot_lightblue,
                         "Extrapolation: Lognormal" = plot_purple)

curves_lines_list_v2 <- c("COU-AA-301: Original" = "solid",
                       "External data (pooled)" = "solid",
                       "Extrapolation: Weibull" = "dashed",
                       "Extrapolation: Gompertz" = "dashed",
                       "Extrapolation: Generalised gamma" = "dashed",
                       "Extrapolation: Exponential" = "dashed",
                       "Extrapolation: Log-logistic" = "dashed",
                       "Extrapolation: Lognormal" = "dashed")

curves_fill_list_v2 <- c("COU-AA-301: Original" = "black",
                      "External data (pooled)" = "grey60",
                      "Extrapolation: Weibull" = "white",
                      "Extrapolation: Gompertz" = "white",
                      "Extrapolation: Generalised gamma" = "white",
                      "Extrapolation: Exponential" = "white",
                      "Extrapolation: Log-logistic" = "white",
                      "Extrapolation: Lognormal" = "white")


KM_active_compareoldext_basicmodels <- KM_active_oldexternal_ci_10$plot +
  geom_line(data = times.expo.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) +
  scale_linetype_manual(values=curves_lines_list_v2) +
  scale_colour_manual(values=curves_colours_list_v2) +
  scale_fill_manual(values=curves_fill_list_v2)

KM_placebo_compareoldext_basicmodels <- KM_placebo_oldexternal_ci_10$plot +
  geom_line(data = times.expo.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) +
  scale_linetype_manual(values=curves_lines_list_v2) +
  scale_colour_manual(values=curves_colours_list_v2) +
  scale_fill_manual(values=curves_fill_list_v2)

KM.list.active.single <-list(Original = KM.est.active.pivotal.original)

survival.data.pivotal.original.active$Abi <- "COU-AA-301: Original"

curves_colours_list2 <- c("COU-AA-301: Original" = "black",
                         "Extrapolation: Weibull" = plot_pink,
                         "Extrapolation: Gompertz" = plot_forestgreen,
                         "Extrapolation: Generalised gamma" = plot_bronze,
                         "Extrapolation: Exponential" = plot_red,
                         "Extrapolation: Log-logistic" = plot_lightblue,
                         "Extrapolation: Lognormal" = plot_purple)

curves_lines_list2 <- c("COU-AA-301: Original" = "solid",
                        "Extrapolation: Weibull" = "dashed",
                        "Extrapolation: Gompertz" = "dashed",
                        "Extrapolation: Generalised gamma" = "dashed",
                        "Extrapolation: Exponential" = "dashed",
                        "Extrapolation: Log-logistic" = "dashed",
                        "Extrapolation: Lognormal" = "dashed")

curves_fill_list2 <- c("COU-AA-301: Original" = "black",
                        "Extrapolation: Weibull" = "white",
                        "Extrapolation: Gompertz" = "white",
                        "Extrapolation: Generalised gamma" = "white",
                        "Extrapolation: Exponential" = "white",
                        "Extrapolation: Log-logistic" = "white",
                        "Extrapolation: Lognormal" = "white")

KM_active_original <- ggsurvplot(KM.est.active.pivotal.original,
                                 survival.data.pivotal.original.active,
                                 palette=c("black"),
                                 linetype=curves_lines_list2,
                                 break.x.by = 0.5,
                                 xlim=c(0,4),
                                 surv.scale="percent",
                                 ylab="Overall survival",
                                 xlab="Time (years)",
                                 fontsize = 3.5,
                                 font.x= c(10, "bold", "black"),
                                 font.y= c(10, "bold", "black"),
                                 censor.size=2,
                                 break.y.by=0.1,
                                 font.tickslab= c(10, "plain", "black"),
                                 legend="bottom",
                                 legend.title="",
                                 legend.labs=c("COU-AA-301: Original")
)

KM_active_original_compare_basicmodels <- KM_active_original$plot +
  geom_line(data = times.expo.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype=guide_legend(title="", nrow=3), fill = guide_legend(title="", nrow=3)) +
  scale_colour_manual(values=curves_colours_list2) +
  scale_linetype_manual(values=curves_lines_list2)+
  scale_fill_manual(values=curves_fill_list2)

KM_placebo_original <- ggsurvplot(KM.est.placebo.pivotal.original,
                                 survival.data.pivotal.original.placebo,
                                 palette=c("black"),
                                 linetype=curves_lines_list2,
                                 break.x.by = 0.5,
                                 xlim=c(0,4),
                                 surv.scale="percent",
                                 ylab="Overall survival",
                                 xlab="Time (years)",
                                 fontsize = 3.5,
                                 font.x= c(10, "bold", "black"),
                                 font.y= c(10, "bold", "black"),
                                 censor.size=2,
                                 break.y.by=0.1,
                                 font.tickslab= c(10, "plain", "black"),
                                 legend="bottom",
                                 legend.title="",
                                 legend.labs=c("COU-AA-301: Original")
)

KM_placebo_original_compare_basicmodels <- KM_placebo_original$plot +
  geom_line(data = times.expo.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype=guide_legend(title="", nrow=3), fill=guide_legend(title="", nrow=3)) +
  scale_colour_manual(values=curves_colours_list2) +
  scale_linetype_manual(values=curves_lines_list2) +
  scale_fill_manual(values=curves_fill_list2)

KM_active_compare_basicmodels <- KM_active_compare$plot +
  geom_line(data = times.expo.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) +
  scale_linetype_manual(values=curves_lines_list) +
  scale_colour_manual(values=curves_colours_list) +
  scale_fill_manual(values=curves_fill_list)

KM_placebo_compare_basicmodels <- KM_placebo_compare$plot +
  geom_line(data = times.expo.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.ggam.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.gomp.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.lnor.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.llog.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) + 
  scale_colour_manual(values=curves_colours_list) +
  scale_linetype_manual(values=curves_lines_list) +
  scale_fill_manual(values=curves_fill_list)

KM_active_compare_basicmodels_reduced <- KM_active_compare$plot +
  #geom_line(data = times.expo.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  #geom_line(data = times.gomp.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  #geom_line(data = times.lnor.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) +
  scale_linetype_manual(values=curves_lines_list) +
  scale_colour_manual(values=curves_colours_list) +
  scale_fill_manual(values=curves_fill_list)

KM_placebo_compare_basicmodels_reduced <- KM_placebo_compare$plot +
  #geom_line(data = times.expo.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  #geom_line(data = times.weib.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.ggam.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  #geom_line(data = times.gomp.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.lnor.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.llog.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) + 
  scale_colour_manual(values=curves_colours_list) +
  scale_linetype_manual(values=curves_lines_list) +
  scale_fill_manual(values=curves_fill_list)

# Plot basic models fitted to naively pooled external (non-rebased) data on top of KM curve # ----

KM_active_external_basicmodels <- KM_active_external$plot +
  geom_line(data = times.expo.external.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.external.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.external.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.external.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.external.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.external.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  scale_colour_manual(values=curves_colours_list_external) +
  scale_linetype_manual(values=curves_lines_list_external) +
  scale_fill_manual(values=curves_fill_list_external) +
  guides(colour=guide_legend(title="", nrow=4), linetype = guide_legend(title="", nrow=4), fill = guide_legend(title="", nrow=4))

KM_placebo_external_basicmodels <- KM_placebo_external$plot +
  geom_line(data = times.expo.external.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.external.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.external.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.external.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.external.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.external.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  scale_colour_manual(values=curves_colours_list_external) +
  scale_linetype_manual(values=curves_lines_list_external) +
  scale_fill_manual(values=curves_fill_list_external) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4), fill = guide_legend(title="", nrow=4))

# Plot basic models fitted to naively pooled external (rebased) data on top of KM curve # ----

# Taking same colour and line lists from non-rebased fitting

KM_active_external_rebased_basicmodels <- KM_active_external_rebased$plot +
  geom_line(data = times.expo.external.rebased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.external.rebased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.external.rebased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.external.rebased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.external.rebased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.external.rebased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  scale_colour_manual(values=curves_colours_list_external) +
  scale_linetype_manual(values=curves_lines_list_external) +
  scale_fill_manual(values=curves_fill_list_external) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4), fill = guide_legend(title="", nrow=4))

KM_placebo_external_rebased_basicmodels <- KM_placebo_external_rebased$plot +
  geom_line(data = times.expo.external.rebased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.external.rebased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.external.rebased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.external.rebased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.external.rebased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.external.rebased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  scale_colour_manual(values=curves_colours_list_external) +
  scale_linetype_manual(values=curves_lines_list_external) +
  scale_fill_manual(values=curves_fill_list_external) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4), fill = guide_legend(title="", nrow=4))

# Plot smoothed hazards versus basic models fitted to original COU-AA-301 data # ----

hazards.muhaz.pivotal.original.active <- data.frame(est = aa_muhaz_active_original$haz.est, time = aa_muhaz_active_original$est.grid, Abi = "COU-AA-301 (de Bono 2011) - AA, muhaz")
hazards.muhaz.pivotal.original.placebo <- data.frame(est = aa_muhaz_placebo_original$haz.est, time = aa_muhaz_placebo_original$est.grid, Abi = "COU-AA-301 (de Bono 2011) - Placebo, muhaz")

curves_colours_list_hazards <- c("COU-AA-301 (de Bono 2011) - AA, muhaz" = "black",
                                 "COU-AA-301 (de Bono 2011) - Placebo, muhaz" = "black",
                                 "Extrapolation: Weibull" = plot_pink,
                                 "Extrapolation: Gompertz" = plot_forestgreen,
                                 "Extrapolation: Generalised gamma" = plot_bronze,
                                 "Extrapolation: Exponential" = plot_red,
                                 "Extrapolation: Log-logistic" = plot_lightblue,
                                 "Extrapolation: Lognormal" = plot_purple)

curves_lines_list_hazards <- c("COU-AA-301 (de Bono 2011) - AA, muhaz" = "solid",
                               "COU-AA-301 (de Bono 2011) - Placebo, muhaz" = "solid",
                               "Extrapolation: Weibull" = "dashed",
                               "Extrapolation: Gompertz" = "dashed",
                               "Extrapolation: Generalised gamma" = "dashed",
                               "Extrapolation: Exponential" = "dashed",
                               "Extrapolation: Log-logistic" = "dashed",
                               "Extrapolation: Lognormal" = "dashed")

KM_active_compare_hazards <- ggplot(hazards.muhaz.pivotal.original.active, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.gomp.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.lnor.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards) +
  scale_linetype_manual(values=curves_lines_list_hazards)

KM_placebo_compare_hazards <- ggplot(hazards.muhaz.pivotal.original.placebo, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.gomp.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.lnor.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards) +
  scale_linetype_manual(values=curves_lines_list_hazards)

KM_active_compare_basicmodels_hazards_reduced <- ggplot(hazards.muhaz.pivotal.original.active, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  #geom_line(data = hazards.expo.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  #geom_line(data = hazards.gomp.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  #geom_line(data = hazards.lnor.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards) +
  scale_linetype_manual(values=curves_lines_list_hazards)

KM_placebo_compare_basicmodels_hazards_reduced <- ggplot(hazards.muhaz.pivotal.original.placebo, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  #geom_line(data = hazards.expo.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  #geom_line(data = hazards.weib.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  #geom_line(data = hazards.gomp.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.lnor.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards) +
  scale_linetype_manual(values=curves_lines_list_hazards)

# Plot smoothed hazards for updated or external data versus basic models fitted to original COU-AA-301 data # ----

aa_muhaz_active_updated = muhaz(time = survival.data.by.arm[[7]]$Years, delta = survival.data.by.arm[[7]]$Event, max.time = max(survival.data.by.arm[[7]]$Years)*muhaz_mult)
aa_muhaz_placebo_updated = muhaz(time = survival.data.by.arm[[8]]$Years, delta = survival.data.by.arm[[8]]$Event, max.time = max(survival.data.by.arm[[8]]$Years)*muhaz_mult)

hazards.muhaz.pivotal.updated.active <- data.frame(est = aa_muhaz_active_updated$haz.est, time = aa_muhaz_active_updated$est.grid, Abi = "COU-AA-301 (Fizazi 2012) - AA, muhaz")
hazards.muhaz.pivotal.updated.placebo <- data.frame(est = aa_muhaz_placebo_updated$haz.est, time = aa_muhaz_placebo_updated$est.grid, Abi = "COU-AA-301 (Fizazi 2012) - Placebo, muhaz")

curves_colours_list_hazards_updated <- c("COU-AA-301 (Fizazi 2012) - AA, muhaz" = "black",
                                 "COU-AA-301 (Fizazi 2012) - Placebo, muhaz" = "black",
                                 "Extrapolation: Weibull" = plot_pink,
                                 "Extrapolation: Gompertz" = plot_forestgreen,
                                 "Extrapolation: Generalised gamma" = plot_bronze,
                                 "Extrapolation: Exponential" = plot_red,
                                 "Extrapolation: Log-logistic" = plot_lightblue,
                                 "Extrapolation: Lognormal" = plot_purple)

curves_lines_list_hazards_updated <- c("COU-AA-301 (Fizazi 2012) - AA, muhaz" = "solid",
                               "COU-AA-301 (Fizazi 2012) - Placebo, muhaz" = "solid",
                               "Extrapolation: Weibull" = "dashed",
                               "Extrapolation: Gompertz" = "dashed",
                               "Extrapolation: Generalised gamma" = "dashed",
                               "Extrapolation: Exponential" = "dashed",
                               "Extrapolation: Log-logistic" = "dashed",
                               "Extrapolation: Lognormal" = "dashed")

KM_active_compare_hazards_updated <- ggplot(hazards.muhaz.pivotal.updated.active, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.gomp.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.lnor.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards_updated) +
  scale_linetype_manual(values=curves_lines_list_hazards_updated)

KM_placebo_compare_hazards_updated <- ggplot(hazards.muhaz.pivotal.updated.placebo, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.gomp.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.lnor.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards_updated) +
  scale_linetype_manual(values=curves_lines_list_hazards_updated)

hazTemp.pooled.active.plotting.v2 <- data.frame(est = hazTemp.pooled.active$haz.est, time = hazTemp.pooled.active$est.grid, Abi = "External data (pooled), muhaz")
hazTemp.pooled.placebo.plotting.v2 <- data.frame(est = hazTemp.pooled.active$haz.est, time = hazTemp.pooled.active$est.grid, Abi = "External data (pooled), muhaz")

curves_colours_list_hazards_external <- c("External data (pooled), muhaz" = "black",
                                         "Extrapolation: Weibull" = plot_pink,
                                         "Extrapolation: Gompertz" = plot_forestgreen,
                                         "Extrapolation: Generalised gamma" = plot_bronze,
                                         "Extrapolation: Exponential" = plot_red,
                                         "Extrapolation: Log-logistic" = plot_lightblue,
                                         "Extrapolation: Lognormal" = plot_purple)

curves_lines_list_hazards_external <- c("External data (pooled), muhaz" = "solid",
                                       "Extrapolation: Weibull" = "dashed",
                                       "Extrapolation: Gompertz" = "dashed",
                                       "Extrapolation: Generalised gamma" = "dashed",
                                       "Extrapolation: Exponential" = "dashed",
                                       "Extrapolation: Log-logistic" = "dashed",
                                       "Extrapolation: Lognormal" = "dashed")

KM_active_compare_hazards_external <- ggplot(hazTemp.pooled.active.plotting.v2, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.gomp.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.lnor.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.pivotal.original.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards_external) +
  scale_linetype_manual(values=curves_lines_list_hazards_external)

KM_placebo_compare_hazards_external <- ggplot(hazTemp.pooled.placebo.plotting.v2, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.gomp.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.lnor.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.pivotal.original.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards_external) +
  scale_linetype_manual(values=curves_lines_list_hazards_external)

# Fitting a piecewise approach # ----

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.expo.piecewise.active <- times.expo.pivotal.original.active

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
  times.expo.piecewise.active.external <- times.expo.external.rebased.active,
  times.expo.piecewise.active.external <- times.expo.external.active))

# Add a column to make a time check for when the cut point is reached
times.expo.piecewise.active$check <- ifelse(times.expo.piecewise.active$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.expo.piecewise.active)){
  ifelse(times.expo.piecewise.active$check[i] == 0,
         times.expo.piecewise.active$est[i] <- times.expo.piecewise.active$est[i],
         times.expo.piecewise.active$est[i] <- times.expo.piecewise.active$est[i - 1] * times.expo.piecewise.active.external$est[i]/times.expo.piecewise.active.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.ggam.piecewise.active <- times.ggam.pivotal.original.active

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.ggam.piecewise.active.external <- times.ggam.external.rebased.active,
                 times.ggam.piecewise.active.external <- times.ggam.external.active))

# Add a column to make a time check for when the cut point is reached
times.ggam.piecewise.active$check <- ifelse(times.ggam.piecewise.active$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.ggam.piecewise.active)){
  ifelse(times.ggam.piecewise.active$check[i] == 0,
         times.ggam.piecewise.active$est[i] <- times.ggam.piecewise.active$est[i],
         times.ggam.piecewise.active$est[i] <- times.ggam.piecewise.active$est[i - 1] * times.ggam.piecewise.active.external$est[i]/times.ggam.piecewise.active.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.gomp.piecewise.active <- times.gomp.pivotal.original.active

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.gomp.piecewise.active.external <- times.gomp.external.rebased.active,
                 times.gomp.piecewise.active.external <- times.gomp.external.active))

# Add a column to make a time check for when the cut point is reached
times.gomp.piecewise.active$check <- ifelse(times.gomp.piecewise.active$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.gomp.piecewise.active)){
  ifelse(times.gomp.piecewise.active$check[i] == 0,
         times.gomp.piecewise.active$est[i] <- times.gomp.piecewise.active$est[i],
         times.gomp.piecewise.active$est[i] <- times.gomp.piecewise.active$est[i - 1] * times.gomp.piecewise.active.external$est[i]/times.gomp.piecewise.active.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.llog.piecewise.active <- times.llog.pivotal.original.active

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.llog.piecewise.active.external <- times.llog.external.rebased.active,
                 times.llog.piecewise.active.external <- times.llog.external.active))

# Add a column to make a time check for when the cut point is reached
times.llog.piecewise.active$check <- ifelse(times.llog.piecewise.active$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.llog.piecewise.active)){
  ifelse(times.llog.piecewise.active$check[i] == 0,
         times.llog.piecewise.active$est[i] <- times.llog.piecewise.active$est[i],
         times.llog.piecewise.active$est[i] <- times.llog.piecewise.active$est[i - 1] * times.llog.piecewise.active.external$est[i]/times.llog.piecewise.active.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.lnor.piecewise.active <- times.lnor.pivotal.original.active

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.lnor.piecewise.active.external <- times.lnor.external.rebased.active,
                 times.lnor.piecewise.active.external <- times.lnor.external.active))

# Add a column to make a time check for when the cut point is reached
times.lnor.piecewise.active$check <- ifelse(times.lnor.piecewise.active$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.lnor.piecewise.active)){
  ifelse(times.lnor.piecewise.active$check[i] == 0,
         times.lnor.piecewise.active$est[i] <- times.lnor.piecewise.active$est[i],
         times.lnor.piecewise.active$est[i] <- times.lnor.piecewise.active$est[i - 1] * times.lnor.piecewise.active.external$est[i]/times.lnor.piecewise.active.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.weib.piecewise.active <- times.weib.pivotal.original.active

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.weib.piecewise.active.external <- times.weib.external.rebased.active,
                 times.weib.piecewise.active.external <- times.weib.external.active))

# Add a column to make a time check for when the cut point is reached
times.weib.piecewise.active$check <- ifelse(times.weib.piecewise.active$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.weib.piecewise.active)){
  ifelse(times.weib.piecewise.active$check[i] == 0,
         times.weib.piecewise.active$est[i] <- times.weib.piecewise.active$est[i],
         times.weib.piecewise.active$est[i] <- times.weib.piecewise.active$est[i - 1] * times.weib.piecewise.active.external$est[i]/times.weib.piecewise.active.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.expo.piecewise.placebo <- times.expo.pivotal.original.placebo

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.expo.piecewise.placebo.external <- times.expo.external.rebased.placebo,
                 times.expo.piecewise.placebo.external <- times.expo.external.placebo))

# Add a column to make a time check for when the cut point is reached
times.expo.piecewise.placebo$check <- ifelse(times.expo.piecewise.placebo$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.expo.piecewise.placebo)){
  ifelse(times.expo.piecewise.placebo$check[i] == 0,
         times.expo.piecewise.placebo$est[i] <- times.expo.piecewise.placebo$est[i],
         times.expo.piecewise.placebo$est[i] <- times.expo.piecewise.placebo$est[i - 1] * times.expo.piecewise.placebo.external$est[i]/times.expo.piecewise.placebo.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.ggam.piecewise.placebo <- times.ggam.pivotal.original.placebo

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.ggam.piecewise.placebo.external <- times.ggam.external.rebased.placebo,
                 times.ggam.piecewise.placebo.external <- times.ggam.external.placebo))

# Add a column to make a time check for when the cut point is reached
times.ggam.piecewise.placebo$check <- ifelse(times.ggam.piecewise.placebo$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.ggam.piecewise.placebo)){
  ifelse(times.ggam.piecewise.placebo$check[i] == 0,
         times.ggam.piecewise.placebo$est[i] <- times.ggam.piecewise.placebo$est[i],
         times.ggam.piecewise.placebo$est[i] <- times.ggam.piecewise.placebo$est[i - 1] * times.ggam.piecewise.placebo.external$est[i]/times.ggam.piecewise.placebo.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.gomp.piecewise.placebo <- times.gomp.pivotal.original.placebo

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.gomp.piecewise.placebo.external <- times.gomp.external.rebased.placebo,
                 times.gomp.piecewise.placebo.external <- times.gomp.external.placebo))

# Add a column to make a time check for when the cut point is reached
times.gomp.piecewise.placebo$check <- ifelse(times.gomp.piecewise.placebo$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.gomp.piecewise.placebo)){
  ifelse(times.gomp.piecewise.placebo$check[i] == 0,
         times.gomp.piecewise.placebo$est[i] <- times.gomp.piecewise.placebo$est[i],
         times.gomp.piecewise.placebo$est[i] <- times.gomp.piecewise.placebo$est[i - 1] * times.gomp.piecewise.placebo.external$est[i]/times.gomp.piecewise.placebo.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.llog.piecewise.placebo <- times.llog.pivotal.original.placebo

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.llog.piecewise.placebo.external <- times.llog.external.rebased.placebo,
                 times.llog.piecewise.placebo.external <- times.llog.external.placebo))

# Add a column to make a time check for when the cut point is reached
times.llog.piecewise.placebo$check <- ifelse(times.llog.piecewise.placebo$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.llog.piecewise.placebo)){
  ifelse(times.llog.piecewise.placebo$check[i] == 0,
         times.llog.piecewise.placebo$est[i] <- times.llog.piecewise.placebo$est[i],
         times.llog.piecewise.placebo$est[i] <- times.llog.piecewise.placebo$est[i - 1] * times.llog.piecewise.placebo.external$est[i]/times.llog.piecewise.placebo.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.lnor.piecewise.placebo <- times.lnor.pivotal.original.placebo

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.lnor.piecewise.placebo.external <- times.lnor.external.rebased.placebo,
                 times.lnor.piecewise.placebo.external <- times.lnor.external.placebo))

# Add a column to make a time check for when the cut point is reached
times.lnor.piecewise.placebo$check <- ifelse(times.lnor.piecewise.placebo$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.lnor.piecewise.placebo)){
  ifelse(times.lnor.piecewise.placebo$check[i] == 0,
         times.lnor.piecewise.placebo$est[i] <- times.lnor.piecewise.placebo$est[i],
         times.lnor.piecewise.placebo$est[i] <- times.lnor.piecewise.placebo$est[i - 1] * times.lnor.piecewise.placebo.external$est[i]/times.lnor.piecewise.placebo.external$est[i - 1])}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.weib.piecewise.placebo <- times.weib.pivotal.original.placebo

# Make selection based on whether to use rebased curves or not
# Invisible stops console outputting values
invisible(ifelse(use_rebased == 1,
                 times.weib.piecewise.placebo.external <- times.weib.external.rebased.placebo,
                 times.weib.piecewise.placebo.external <- times.weib.external.placebo))

# Add a column to make a time check for when the cut point is reached
times.weib.piecewise.placebo$check <- ifelse(times.weib.piecewise.placebo$time <= piecewise_cutpoint,0,1)

# Estimate piecewise model times using a check for time, and then informing hazard by taking relative % of deaths from other source
for (i in 1:nrow(times.weib.piecewise.placebo)){
  ifelse(times.weib.piecewise.placebo$check[i] == 0,
         times.weib.piecewise.placebo$est[i] <- times.weib.piecewise.placebo$est[i],
         times.weib.piecewise.placebo$est[i] <- times.weib.piecewise.placebo$est[i - 1] * times.weib.piecewise.placebo.external$est[i]/times.weib.piecewise.placebo.external$est[i - 1])}

times.weib.piecewise.active$Abi <- "Piecewise: Weibull"
times.gomp.piecewise.active$Abi <- "Piecewise: Gompertz"
times.ggam.piecewise.active$Abi <- "Piecewise: Generalised gamma"
times.expo.piecewise.active$Abi <- "Piecewise: Exponential"
times.llog.piecewise.active$Abi <- "Piecewise: Log-logistic"
times.lnor.piecewise.active$Abi <- "Piecewise: Lognormal"

times.weib.piecewise.placebo$Abi <- "Piecewise: Weibull"
times.gomp.piecewise.placebo$Abi <- "Piecewise: Gompertz"
times.ggam.piecewise.placebo$Abi <- "Piecewise: Generalised gamma"
times.expo.piecewise.placebo$Abi <- "Piecewise: Exponential"
times.llog.piecewise.placebo$Abi <- "Piecewise: Log-logistic"
times.lnor.piecewise.placebo$Abi <- "Piecewise: Lognormal"

hazards.expo.piecewise.active <- hazards.expo.pivotal.original.active
hazards.weib.piecewise.active <- hazards.weib.pivotal.original.active
hazards.ggam.piecewise.active <- hazards.ggam.pivotal.original.active
hazards.gomp.piecewise.active <- hazards.gomp.pivotal.original.active
hazards.lnor.piecewise.active <- hazards.lnor.pivotal.original.active
hazards.llog.piecewise.active <- hazards.llog.pivotal.original.active

hazards.expo.piecewise.active$check <- ifelse(hazards.expo.piecewise.active$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.expo.piecewise.active.external <- hazards.expo.external.rebased.active,
                 hazards.expo.piecewise.active.external <- hazards.expo.external.active))
for (i in 1:nrow(hazards.expo.piecewise.active)){
  ifelse(hazards.expo.piecewise.active$check[i] == 0,
         hazards.expo.piecewise.active$est[i] <- hazards.expo.piecewise.active$est[i],
         hazards.expo.piecewise.active$est[i] <- hazards.expo.piecewise.active.external$est[i])}

hazards.weib.piecewise.active$check <- ifelse(hazards.weib.piecewise.active$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.weib.piecewise.active.external <- hazards.weib.external.rebased.active,
                 hazards.weib.piecewise.active.external <- hazards.weib.external.active))
for (i in 1:nrow(hazards.weib.piecewise.active)){
  ifelse(hazards.weib.piecewise.active$check[i] == 0,
         hazards.weib.piecewise.active$est[i] <- hazards.weib.piecewise.active$est[i],
         hazards.weib.piecewise.active$est[i] <- hazards.weib.piecewise.active.external$est[i])}

hazards.gomp.piecewise.active$check <- ifelse(hazards.gomp.piecewise.active$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.gomp.piecewise.active.external <- hazards.gomp.external.rebased.active,
                 hazards.gomp.piecewise.active.external <- hazards.gomp.external.active))
for (i in 1:nrow(hazards.gomp.piecewise.active)){
  ifelse(hazards.gomp.piecewise.active$check[i] == 0,
         hazards.gomp.piecewise.active$est[i] <- hazards.gomp.piecewise.active$est[i],
         hazards.gomp.piecewise.active$est[i] <- hazards.gomp.piecewise.active.external$est[i])}

hazards.ggam.piecewise.active$check <- ifelse(hazards.ggam.piecewise.active$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.ggam.piecewise.active.external <- hazards.ggam.external.rebased.active,
                 hazards.ggam.piecewise.active.external <- hazards.ggam.external.active))
for (i in 1:nrow(hazards.ggam.piecewise.active)){
  ifelse(hazards.ggam.piecewise.active$check[i] == 0,
         hazards.ggam.piecewise.active$est[i] <- hazards.ggam.piecewise.active$est[i],
         hazards.ggam.piecewise.active$est[i] <- hazards.ggam.piecewise.active.external$est[i])}

hazards.lnor.piecewise.active$check <- ifelse(hazards.lnor.piecewise.active$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.lnor.piecewise.active.external <- hazards.lnor.external.rebased.active,
                 hazards.lnor.piecewise.active.external <- hazards.lnor.external.active))
for (i in 1:nrow(hazards.lnor.piecewise.active)){
  ifelse(hazards.lnor.piecewise.active$check[i] == 0,
         hazards.lnor.piecewise.active$est[i] <- hazards.lnor.piecewise.active$est[i],
         hazards.lnor.piecewise.active$est[i] <- hazards.lnor.piecewise.active.external$est[i])}

hazards.llog.piecewise.active$check <- ifelse(hazards.llog.piecewise.active$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.llog.piecewise.active.external <- hazards.llog.external.rebased.active,
                 hazards.llog.piecewise.active.external <- hazards.llog.external.active))
for (i in 1:nrow(hazards.llog.piecewise.active)){
  ifelse(hazards.llog.piecewise.active$check[i] == 0,
         hazards.llog.piecewise.active$est[i] <- hazards.llog.piecewise.active$est[i],
         hazards.llog.piecewise.active$est[i] <- hazards.llog.piecewise.active.external$est[i])}

hazards.expo.piecewise.active$Abi <- "Piecewise: Exponential"       
hazards.weib.piecewise.active$Abi <- "Piecewise: Weibull"          
hazards.ggam.piecewise.active$Abi <- "Piecewise: Generalised gamma"
hazards.gomp.piecewise.active$Abi <- "Piecewise: Gompertz"              
hazards.lnor.piecewise.active$Abi <- "Piecewise: Lognormal"             
hazards.llog.piecewise.active$Abi <- "Piecewise: Log-logistic"

hazards.expo.piecewise.placebo <- hazards.expo.pivotal.original.placebo
hazards.weib.piecewise.placebo <- hazards.weib.pivotal.original.placebo
hazards.ggam.piecewise.placebo <- hazards.ggam.pivotal.original.placebo
hazards.gomp.piecewise.placebo <- hazards.gomp.pivotal.original.placebo
hazards.lnor.piecewise.placebo <- hazards.lnor.pivotal.original.placebo
hazards.llog.piecewise.placebo <- hazards.llog.pivotal.original.placebo

hazards.expo.piecewise.placebo$check <- ifelse(hazards.expo.piecewise.placebo$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.expo.piecewise.placebo.external <- hazards.expo.external.rebased.placebo,
                 hazards.expo.piecewise.placebo.external <- hazards.expo.external.placebo))
for (i in 1:nrow(hazards.expo.piecewise.placebo)){
  ifelse(hazards.expo.piecewise.placebo$check[i] == 0,
         hazards.expo.piecewise.placebo$est[i] <- hazards.expo.piecewise.placebo$est[i],
         hazards.expo.piecewise.placebo$est[i] <- hazards.expo.piecewise.placebo.external$est[i])}

hazards.weib.piecewise.placebo$check <- ifelse(hazards.weib.piecewise.placebo$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.weib.piecewise.placebo.external <- hazards.weib.external.rebased.placebo,
                 hazards.weib.piecewise.placebo.external <- hazards.weib.external.placebo))
for (i in 1:nrow(hazards.weib.piecewise.placebo)){
  ifelse(hazards.weib.piecewise.placebo$check[i] == 0,
         hazards.weib.piecewise.placebo$est[i] <- hazards.weib.piecewise.placebo$est[i],
         hazards.weib.piecewise.placebo$est[i] <- hazards.weib.piecewise.placebo.external$est[i])}

hazards.gomp.piecewise.placebo$check <- ifelse(hazards.gomp.piecewise.placebo$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.gomp.piecewise.placebo.external <- hazards.gomp.external.rebased.placebo,
                 hazards.gomp.piecewise.placebo.external <- hazards.gomp.external.placebo))
for (i in 1:nrow(hazards.gomp.piecewise.placebo)){
  ifelse(hazards.gomp.piecewise.placebo$check[i] == 0,
         hazards.gomp.piecewise.placebo$est[i] <- hazards.gomp.piecewise.placebo$est[i],
         hazards.gomp.piecewise.placebo$est[i] <- hazards.gomp.piecewise.placebo.external$est[i])}

hazards.ggam.piecewise.placebo$check <- ifelse(hazards.ggam.piecewise.placebo$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.ggam.piecewise.placebo.external <- hazards.ggam.external.rebased.placebo,
                 hazards.ggam.piecewise.placebo.external <- hazards.ggam.external.placebo))
for (i in 1:nrow(hazards.ggam.piecewise.placebo)){
  ifelse(hazards.ggam.piecewise.placebo$check[i] == 0,
         hazards.ggam.piecewise.placebo$est[i] <- hazards.ggam.piecewise.placebo$est[i],
         hazards.ggam.piecewise.placebo$est[i] <- hazards.ggam.piecewise.placebo.external$est[i])}

hazards.lnor.piecewise.placebo$check <- ifelse(hazards.lnor.piecewise.placebo$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.lnor.piecewise.placebo.external <- hazards.lnor.external.rebased.placebo,
                 hazards.lnor.piecewise.placebo.external <- hazards.lnor.external.placebo))
for (i in 1:nrow(hazards.lnor.piecewise.placebo)){
  ifelse(hazards.lnor.piecewise.placebo$check[i] == 0,
         hazards.lnor.piecewise.placebo$est[i] <- hazards.lnor.piecewise.placebo$est[i],
         hazards.lnor.piecewise.placebo$est[i] <- hazards.lnor.piecewise.placebo.external$est[i])}

hazards.llog.piecewise.placebo$check <- ifelse(hazards.llog.piecewise.placebo$time <= piecewise_cutpoint,0,1)
invisible(ifelse(use_rebased == 1,
                 hazards.llog.piecewise.placebo.external <- hazards.llog.external.rebased.placebo,
                 hazards.llog.piecewise.placebo.external <- hazards.llog.external.placebo))
for (i in 1:nrow(hazards.llog.piecewise.placebo)){
  ifelse(hazards.llog.piecewise.placebo$check[i] == 0,
         hazards.llog.piecewise.placebo$est[i] <- hazards.llog.piecewise.placebo$est[i],
         hazards.llog.piecewise.placebo$est[i] <- hazards.llog.piecewise.placebo.external$est[i])}

hazards.expo.piecewise.placebo$Abi <- "Piecewise: Exponential"       
hazards.weib.piecewise.placebo$Abi <- "Piecewise: Weibull"          
hazards.ggam.piecewise.placebo$Abi <- "Piecewise: Generalised gamma"
hazards.gomp.piecewise.placebo$Abi <- "Piecewise: Gompertz"              
hazards.lnor.piecewise.placebo$Abi <- "Piecewise: Lognormal"             
hazards.llog.piecewise.placebo$Abi <- "Piecewise: Log-logistic"

# Plot piecewise models # ----

piecewise_colours_list <- c("COU-AA-301: Original" = "black",
                         "COU-AA-301: Updated" = "grey40",
                         "External data (pooled)" = "grey80",
                         "Piecewise: Weibull" = plot_pink,
                         "Piecewise: Gompertz" = plot_forestgreen,
                         "Piecewise: Generalised gamma" = plot_bronze,
                         "Piecewise: Exponential" = plot_red,
                         "Piecewise: Log-logistic" = plot_lightblue,
                         "Piecewise: Lognormal" = plot_purple)

piecewise_lines_list <- c("COU-AA-301: Original" = "solid",
                       "COU-AA-301: Updated" = "solid",
                       "External data (pooled)" = "solid",
                       "Piecewise: Weibull" = "dashed",
                       "Piecewise: Gompertz" = "dashed",
                       "Piecewise: Generalised gamma" = "dashed",
                       "Piecewise: Exponential" = "dashed",
                       "Piecewise: Log-logistic" = "dashed",
                       "Piecewise: Lognormal" = "dashed")

piecewise_fill_list <- c("COU-AA-301: Original" = "black",
                            "COU-AA-301: Updated" = "grey40",
                            "External data (pooled)" = "grey80",
                            "Piecewise: Weibull" = "white",
                            "Piecewise: Gompertz" = "white",
                            "Piecewise: Generalised gamma" = "white",
                            "Piecewise: Exponential" = "white",
                            "Piecewise: Log-logistic" = "white",
                            "Piecewise: Lognormal" = "white")

KM_active_compare_piecewisemodels <- KM_active_compare$plot +
  geom_line(data = times.expo.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.ggam.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.gomp.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.lnor.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.llog.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill = guide_legend(title="", nrow=3)) + 
  scale_colour_manual(values=piecewise_colours_list) +
  scale_linetype_manual(values=piecewise_lines_list)+
  scale_fill_manual(values=piecewise_fill_list)

KM_placebo_compare_piecewisemodels <- KM_placebo_compare$plot +
  geom_line(data = times.expo.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.ggam.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.gomp.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.lnor.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.llog.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill = guide_legend(title="", nrow=3)) + 
  scale_colour_manual(values=piecewise_colours_list) +
  scale_linetype_manual(values=piecewise_lines_list)+
  scale_fill_manual(values=piecewise_fill_list)

KM_active_compareoldext_piecewisemodels <- KM_active_oldexternal_ci_10$plot +
  geom_line(data = times.expo.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) +
  scale_colour_manual(values=piecewise_colours_list) +
  scale_linetype_manual(values=piecewise_lines_list)+
  scale_fill_manual(values=piecewise_fill_list)

KM_placebo_compareoldext_piecewisemodels <- KM_placebo_oldexternal_ci_10$plot +
  geom_line(data = times.expo.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) +
  scale_colour_manual(values=piecewise_colours_list) +
  scale_linetype_manual(values=piecewise_lines_list)+
  scale_fill_manual(values=piecewise_fill_list)

# Plot smoothed hazards versus piecewise models # ----

curves_colours_list_hazards_piecewise <- c("COU-AA-301 (de Bono 2011) - AA, muhaz" = "black",
                                 "COU-AA-301 (de Bono 2011) - Placebo, muhaz" = "black",
                                 "Piecewise: Weibull" = plot_pink,
                                 "Piecewise: Gompertz" = plot_forestgreen,
                                 "Piecewise: Generalised gamma" = plot_bronze,
                                 "Piecewise: Exponential" = plot_red,
                                 "Piecewise: Log-logistic" = plot_lightblue,
                                 "Piecewise: Lognormal" = plot_purple)

curves_lines_list_hazards_piecewise <- c("COU-AA-301 (de Bono 2011) - AA, muhaz" = "solid",
                               "COU-AA-301 (de Bono 2011) - Placebo, muhaz" = "solid",
                               "Piecewise: Weibull" = "dashed",
                               "Piecewise: Gompertz" = "dashed",
                               "Piecewise: Generalised gamma" = "dashed",
                               "Piecewise: Exponential" = "dashed",
                               "Piecewise: Log-logistic" = "dashed",
                               "Piecewise: Lognormal" = "dashed")

KM_active_compare_hazards_piecewise <- ggplot(hazards.muhaz.pivotal.original.active, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.ggam.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.gomp.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.lnor.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.llog.piecewise.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards_piecewise) +
  scale_linetype_manual(values=curves_lines_list_hazards_piecewise)

KM_placebo_compare_hazards_piecewise <- ggplot(hazards.muhaz.pivotal.original.placebo, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.gomp.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.lnor.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.piecewise.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards_piecewise) +
  scale_linetype_manual(values=curves_lines_list_hazards_piecewise)

# Fitting HR-based models (rebased) # ----

# Obtain HR for active arm
survival.data.active.original.external.rebased <- rbind(survival.data.pivotal.original.active, survival.data.external.active)
survival.data.active.original.external.rebased <- subset(survival.data.active.original.external.rebased, survival.data.active.original.external.rebased$Years > piecewise_cutpoint + error_margin)
survival.data.active.original.external.rebased$Years <- survival.data.active.original.external.rebased$Years - piecewise_cutpoint
PH.test.active <- coxph(Surv(Years,Event) ~ Abi, survival.data.active.original.external.rebased)
HR <- as.numeric(exp(coef(PH.test.active)))
HR <- 1/HR
HR.active <- HR

ph_plot_test_active <- survfit(Surv(Years,Event) ~ Abi, survival.data.active.original.external.rebased)
plot(ph_plot_test_active, col=c("red","blue"))
plot(ph_plot_test_active, fun = "cloglog", col=c("red","blue"))

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.expo.HRbased.active <- times.expo.pivotal.original.active
# Identify external source
times.expo.HRbased.active.external <- times.expo.external.rebased.active
# Add a column to make a time check for when the cut point is reached
times.expo.HRbased.active$check <- ifelse(times.expo.HRbased.active$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.expo.HRbased.active)){
  ifelse(times.expo.HRbased.active$check[i] == 0,
         times.expo.HRbased.active$est[i] <- times.expo.HRbased.active$est[i],
         times.expo.HRbased.active$est[i] <- times.expo.HRbased.active$est[i - 1] * (times.expo.HRbased.active.external$est[i]^HR)/(times.expo.HRbased.active.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.weib.HRbased.active <- times.weib.pivotal.original.active
# Identify external source
times.weib.HRbased.active.external <- times.weib.external.rebased.active
# Add a column to make a time check for when the cut point is reached
times.weib.HRbased.active$check <- ifelse(times.weib.HRbased.active$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.weib.HRbased.active)){
  ifelse(times.weib.HRbased.active$check[i] == 0,
         times.weib.HRbased.active$est[i] <- times.weib.HRbased.active$est[i],
         times.weib.HRbased.active$est[i] <- times.weib.HRbased.active$est[i - 1] * (times.weib.HRbased.active.external$est[i]^HR)/(times.weib.HRbased.active.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.gomp.HRbased.active <- times.gomp.pivotal.original.active
# Identify external source
times.gomp.HRbased.active.external <- times.gomp.external.rebased.active
# Add a column to make a time check for when the cut point is reached
times.gomp.HRbased.active$check <- ifelse(times.gomp.HRbased.active$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.gomp.HRbased.active)){
  ifelse(times.gomp.HRbased.active$check[i] == 0,
         times.gomp.HRbased.active$est[i] <- times.gomp.HRbased.active$est[i],
         times.gomp.HRbased.active$est[i] <- times.gomp.HRbased.active$est[i - 1] * (times.gomp.HRbased.active.external$est[i]^HR)/(times.gomp.HRbased.active.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.ggam.HRbased.active <- times.ggam.pivotal.original.active
# Identify external source
times.ggam.HRbased.active.external <- times.ggam.external.rebased.active
# Add a column to make a time check for when the cut point is reached
times.ggam.HRbased.active$check <- ifelse(times.ggam.HRbased.active$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.ggam.HRbased.active)){
  ifelse(times.ggam.HRbased.active$check[i] == 0,
         times.ggam.HRbased.active$est[i] <- times.ggam.HRbased.active$est[i],
         times.ggam.HRbased.active$est[i] <- times.ggam.HRbased.active$est[i - 1] * (times.ggam.HRbased.active.external$est[i]^HR)/(times.ggam.HRbased.active.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.lnor.HRbased.active <- times.lnor.pivotal.original.active
# Identify external source
times.lnor.HRbased.active.external <- times.lnor.external.rebased.active
# Add a column to make a time check for when the cut point is reached
times.lnor.HRbased.active$check <- ifelse(times.lnor.HRbased.active$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.lnor.HRbased.active)){
  ifelse(times.lnor.HRbased.active$check[i] == 0,
         times.lnor.HRbased.active$est[i] <- times.lnor.HRbased.active$est[i],
         times.lnor.HRbased.active$est[i] <- times.lnor.HRbased.active$est[i - 1] * (times.lnor.HRbased.active.external$est[i]^HR)/(times.lnor.HRbased.active.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.llog.HRbased.active <- times.llog.pivotal.original.active
# Identify external source
times.llog.HRbased.active.external <- times.llog.external.rebased.active
# Add a column to make a time check for when the cut point is reached
times.llog.HRbased.active$check <- ifelse(times.llog.HRbased.active$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.llog.HRbased.active)){
  ifelse(times.llog.HRbased.active$check[i] == 0,
         times.llog.HRbased.active$est[i] <- times.llog.HRbased.active$est[i],
         times.llog.HRbased.active$est[i] <- times.llog.HRbased.active$est[i - 1] * (times.llog.HRbased.active.external$est[i]^HR)/(times.llog.HRbased.active.external$est[i - 1]^HR))}

# Obtain HR for placebo arm
survival.data.placebo.original.external.rebased <- rbind(survival.data.pivotal.original.placebo, survival.data.external.placebo)
survival.data.placebo.original.external.rebased <- subset(survival.data.placebo.original.external.rebased, survival.data.placebo.original.external.rebased$Years > piecewise_cutpoint + error_margin)
survival.data.placebo.original.external.rebased$Years <- survival.data.placebo.original.external.rebased$Years - piecewise_cutpoint
PH.test.placebo <- coxph(Surv(Years,Event) ~ Abi, survival.data.placebo.original.external.rebased)
HR <- as.numeric(exp(coef(PH.test.placebo)))
HR <- 1/HR
HR.placebo <- HR

ph_plot_test_placebo <- survfit(Surv(Years,Event) ~ Abi, survival.data.placebo.original.external.rebased)
plot(ph_plot_test_placebo, col=c("red","blue"))
plot(ph_plot_test_placebo, fun = "cloglog", col=c("red","blue"))

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.expo.HRbased.placebo <- times.expo.pivotal.original.placebo
# Identify external source
times.expo.HRbased.placebo.external <- times.expo.external.rebased.placebo
# Add a column to make a time check for when the cut point is reached
times.expo.HRbased.placebo$check <- ifelse(times.expo.HRbased.placebo$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.expo.HRbased.placebo)){
  ifelse(times.expo.HRbased.placebo$check[i] == 0,
         times.expo.HRbased.placebo$est[i] <- times.expo.HRbased.placebo$est[i],
         times.expo.HRbased.placebo$est[i] <- times.expo.HRbased.placebo$est[i - 1] * (times.expo.HRbased.placebo.external$est[i]^HR)/(times.expo.HRbased.placebo.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.weib.HRbased.placebo <- times.weib.pivotal.original.placebo
# Identify external source
times.weib.HRbased.placebo.external <- times.weib.external.rebased.placebo
# Add a column to make a time check for when the cut point is reached
times.weib.HRbased.placebo$check <- ifelse(times.weib.HRbased.placebo$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.weib.HRbased.placebo)){
  ifelse(times.weib.HRbased.placebo$check[i] == 0,
         times.weib.HRbased.placebo$est[i] <- times.weib.HRbased.placebo$est[i],
         times.weib.HRbased.placebo$est[i] <- times.weib.HRbased.placebo$est[i - 1] * (times.weib.HRbased.placebo.external$est[i]^HR)/(times.weib.HRbased.placebo.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.gomp.HRbased.placebo <- times.gomp.pivotal.original.placebo
# Identify external source
times.gomp.HRbased.placebo.external <- times.gomp.external.rebased.placebo
# Add a column to make a time check for when the cut point is reached
times.gomp.HRbased.placebo$check <- ifelse(times.gomp.HRbased.placebo$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.gomp.HRbased.placebo)){
  ifelse(times.gomp.HRbased.placebo$check[i] == 0,
         times.gomp.HRbased.placebo$est[i] <- times.gomp.HRbased.placebo$est[i],
         times.gomp.HRbased.placebo$est[i] <- times.gomp.HRbased.placebo$est[i - 1] * (times.gomp.HRbased.placebo.external$est[i]^HR)/(times.gomp.HRbased.placebo.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.ggam.HRbased.placebo <- times.ggam.pivotal.original.placebo
# Identify external source
times.ggam.HRbased.placebo.external <- times.ggam.external.rebased.placebo
# Add a column to make a time check for when the cut point is reached
times.ggam.HRbased.placebo$check <- ifelse(times.ggam.HRbased.placebo$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.ggam.HRbased.placebo)){
  ifelse(times.ggam.HRbased.placebo$check[i] == 0,
         times.ggam.HRbased.placebo$est[i] <- times.ggam.HRbased.placebo$est[i],
         times.ggam.HRbased.placebo$est[i] <- times.ggam.HRbased.placebo$est[i - 1] * (times.ggam.HRbased.placebo.external$est[i]^HR)/(times.ggam.HRbased.placebo.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.lnor.HRbased.placebo <- times.lnor.pivotal.original.placebo
# Identify external source
times.lnor.HRbased.placebo.external <- times.lnor.external.rebased.placebo
# Add a column to make a time check for when the cut point is reached
times.lnor.HRbased.placebo$check <- ifelse(times.lnor.HRbased.placebo$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.lnor.HRbased.placebo)){
  ifelse(times.lnor.HRbased.placebo$check[i] == 0,
         times.lnor.HRbased.placebo$est[i] <- times.lnor.HRbased.placebo$est[i],
         times.lnor.HRbased.placebo$est[i] <- times.lnor.HRbased.placebo$est[i - 1] * (times.lnor.HRbased.placebo.external$est[i]^HR)/(times.lnor.HRbased.placebo.external$est[i - 1]^HR))}

# Start by setting times to be same as basic model fitted to COU-AA-301 data
times.llog.HRbased.placebo <- times.llog.pivotal.original.placebo
# Identify external source
times.llog.HRbased.placebo.external <- times.llog.external.rebased.placebo
# Add a column to make a time check for when the cut point is reached
times.llog.HRbased.placebo$check <- ifelse(times.llog.HRbased.placebo$time <= piecewise_cutpoint,0,1)
# Set up as per piecewise model
for (i in 1:nrow(times.llog.HRbased.placebo)){
  ifelse(times.llog.HRbased.placebo$check[i] == 0,
         times.llog.HRbased.placebo$est[i] <- times.llog.HRbased.placebo$est[i],
         times.llog.HRbased.placebo$est[i] <- times.llog.HRbased.placebo$est[i - 1] * (times.llog.HRbased.placebo.external$est[i]^HR)/(times.llog.HRbased.placebo.external$est[i - 1]^HR))}

times.weib.HRbased.active$Abi <- "HR-based: Weibull"
times.gomp.HRbased.active$Abi <- "HR-based: Gompertz"
times.ggam.HRbased.active$Abi <- "HR-based: Generalised gamma"
times.expo.HRbased.active$Abi <- "HR-based: Exponential"
times.llog.HRbased.active$Abi <- "HR-based: Log-logistic"
times.lnor.HRbased.active$Abi <- "HR-based: Lognormal"

times.weib.HRbased.placebo$Abi <- "HR-based: Weibull"
times.gomp.HRbased.placebo$Abi <- "HR-based: Gompertz"
times.ggam.HRbased.placebo$Abi <- "HR-based: Generalised gamma"
times.expo.HRbased.placebo$Abi <- "HR-based: Exponential"
times.llog.HRbased.placebo$Abi <- "HR-based: Log-logistic"
times.lnor.HRbased.placebo$Abi <- "HR-based: Lognormal"

hazards.expo.HRbased.active <- hazards.expo.pivotal.original.active
hazards.expo.HRbased.active$check <- ifelse(hazards.expo.HRbased.active$time <= piecewise_cutpoint,0,1)

HR <- HR.active

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.expo.HRbased.active)){
  ifelse(hazards.expo.HRbased.active$check[i] == 0,
         hazards.expo.HRbased.active$est[i] <- hazards.expo.piecewise.active$est[i]^1,
         hazards.expo.HRbased.active$est[i] <- hazards.expo.piecewise.active$est[i]*HR)}

hazards.weib.HRbased.active <- hazards.weib.pivotal.original.active
hazards.weib.HRbased.active$check <- ifelse(hazards.weib.HRbased.active$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.weib.HRbased.active)){
  ifelse(hazards.weib.HRbased.active$check[i] == 0,
         hazards.weib.HRbased.active$est[i] <- hazards.weib.piecewise.active$est[i]^1,
         hazards.weib.HRbased.active$est[i] <- hazards.weib.piecewise.active$est[i]*HR)}

hazards.gomp.HRbased.active <- hazards.gomp.pivotal.original.active
hazards.gomp.HRbased.active$check <- ifelse(hazards.gomp.HRbased.active$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.gomp.HRbased.active)){
  ifelse(hazards.gomp.HRbased.active$check[i] == 0,
         hazards.gomp.HRbased.active$est[i] <- hazards.gomp.piecewise.active$est[i]^1,
         hazards.gomp.HRbased.active$est[i] <- hazards.gomp.piecewise.active$est[i]*HR)}

hazards.lnor.HRbased.active <- hazards.lnor.pivotal.original.active
hazards.lnor.HRbased.active$check <- ifelse(hazards.lnor.HRbased.active$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.lnor.HRbased.active)){
  ifelse(hazards.lnor.HRbased.active$check[i] == 0,
         hazards.lnor.HRbased.active$est[i] <- hazards.lnor.piecewise.active$est[i]^1,
         hazards.lnor.HRbased.active$est[i] <- hazards.lnor.piecewise.active$est[i]*HR)}

hazards.llog.HRbased.active <- hazards.llog.pivotal.original.active
hazards.llog.HRbased.active$check <- ifelse(hazards.llog.HRbased.active$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.llog.HRbased.active)){
  ifelse(hazards.llog.HRbased.active$check[i] == 0,
         hazards.llog.HRbased.active$est[i] <- hazards.llog.piecewise.active$est[i]^1,
         hazards.llog.HRbased.active$est[i] <- hazards.llog.piecewise.active$est[i]*HR)}

hazards.ggam.HRbased.active <- hazards.ggam.pivotal.original.active
hazards.ggam.HRbased.active$check <- ifelse(hazards.ggam.HRbased.active$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.ggam.HRbased.active)){
  ifelse(hazards.ggam.HRbased.active$check[i] == 0,
         hazards.ggam.HRbased.active$est[i] <- hazards.ggam.piecewise.active$est[i]^1,
         hazards.ggam.HRbased.active$est[i] <- hazards.ggam.piecewise.active$est[i]*HR)}

hazards.expo.HRbased.placebo <- hazards.expo.pivotal.original.placebo
hazards.expo.HRbased.placebo$check <- ifelse(hazards.expo.HRbased.placebo$time <= piecewise_cutpoint,0,1)
HR <- HR.placebo

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.expo.HRbased.placebo)){
  ifelse(hazards.expo.HRbased.placebo$check[i] == 0,
         hazards.expo.HRbased.placebo$est[i] <- hazards.expo.piecewise.placebo$est[i]^1,
         hazards.expo.HRbased.placebo$est[i] <- hazards.expo.piecewise.placebo$est[i]*HR)}

hazards.weib.HRbased.placebo <- hazards.weib.pivotal.original.placebo
hazards.weib.HRbased.placebo$check <- ifelse(hazards.weib.HRbased.placebo$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.weib.HRbased.placebo)){
  ifelse(hazards.weib.HRbased.placebo$check[i] == 0,
         hazards.weib.HRbased.placebo$est[i] <- hazards.weib.piecewise.placebo$est[i]^1,
         hazards.weib.HRbased.placebo$est[i] <- hazards.weib.piecewise.placebo$est[i]*HR)}

hazards.gomp.HRbased.placebo <- hazards.gomp.pivotal.original.placebo
hazards.gomp.HRbased.placebo$check <- ifelse(hazards.gomp.HRbased.placebo$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.gomp.HRbased.placebo)){
  ifelse(hazards.gomp.HRbased.placebo$check[i] == 0,
         hazards.gomp.HRbased.placebo$est[i] <- hazards.gomp.piecewise.placebo$est[i]^1,
         hazards.gomp.HRbased.placebo$est[i] <- hazards.gomp.piecewise.placebo$est[i]*HR)}

hazards.lnor.HRbased.placebo <- hazards.lnor.pivotal.original.placebo
hazards.lnor.HRbased.placebo$check <- ifelse(hazards.lnor.HRbased.placebo$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.lnor.HRbased.placebo)){
  ifelse(hazards.lnor.HRbased.placebo$check[i] == 0,
         hazards.lnor.HRbased.placebo$est[i] <- hazards.lnor.piecewise.placebo$est[i]^1,
         hazards.lnor.HRbased.placebo$est[i] <- hazards.lnor.piecewise.placebo$est[i]*HR)}

hazards.llog.HRbased.placebo <- hazards.llog.pivotal.original.placebo
hazards.llog.HRbased.placebo$check <- ifelse(hazards.llog.HRbased.placebo$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.llog.HRbased.placebo)){
  ifelse(hazards.llog.HRbased.placebo$check[i] == 0,
         hazards.llog.HRbased.placebo$est[i] <- hazards.llog.piecewise.placebo$est[i]^1,
         hazards.llog.HRbased.placebo$est[i] <- hazards.llog.piecewise.placebo$est[i]*HR)}

hazards.ggam.HRbased.placebo <- hazards.ggam.pivotal.original.placebo
hazards.ggam.HRbased.placebo$check <- ifelse(hazards.ggam.HRbased.placebo$time <= piecewise_cutpoint,0,1)

# Adjust extrapolated component based on estimated hazard ratio
for (i in 1:nrow(times.ggam.HRbased.placebo)){
  ifelse(hazards.ggam.HRbased.placebo$check[i] == 0,
         hazards.ggam.HRbased.placebo$est[i] <- hazards.ggam.piecewise.placebo$est[i]^1,
         hazards.ggam.HRbased.placebo$est[i] <- hazards.ggam.piecewise.placebo$est[i]*HR)}

hazards.weib.HRbased.active$Abi <- "HR-based: Weibull"
hazards.gomp.HRbased.active$Abi <- "HR-based: Gompertz"
hazards.ggam.HRbased.active$Abi <- "HR-based: Generalised gamma"
hazards.expo.HRbased.active$Abi <- "HR-based: Exponential"
hazards.llog.HRbased.active$Abi <- "HR-based: Log-logistic"
hazards.lnor.HRbased.active$Abi <- "HR-based: Lognormal"

hazards.weib.HRbased.placebo$Abi <- "HR-based: Weibull"
hazards.gomp.HRbased.placebo$Abi <- "HR-based: Gompertz"
hazards.ggam.HRbased.placebo$Abi <- "HR-based: Generalised gamma"
hazards.expo.HRbased.placebo$Abi <- "HR-based: Exponential"
hazards.llog.HRbased.placebo$Abi <- "HR-based: Log-logistic"
hazards.lnor.HRbased.placebo$Abi <- "HR-based: Lognormal"

# Plot HR-based models (rebased) # ----

HRbased_colours_list <- c("COU-AA-301: Original" = "black",
                            "COU-AA-301: Updated" = "grey40",
                            "External data (pooled)" = "grey80",
                            "HR-based: Weibull" = plot_pink,
                            "HR-based: Gompertz" = plot_forestgreen,
                            "HR-based: Generalised gamma" = plot_bronze,
                            "HR-based: Exponential" = plot_red,
                            "HR-based: Log-logistic" = plot_lightblue,
                            "HR-based: Lognormal" = plot_purple)

HRbased_lines_list <- c("COU-AA-301: Original" = "solid",
                          "COU-AA-301: Updated" = "solid",
                          "External data (pooled)" = "solid",
                          "HR-based: Weibull" = "dashed",
                          "HR-based: Gompertz" = "dashed",
                          "HR-based: Generalised gamma" = "dashed",
                          "HR-based: Exponential" = "dashed",
                          "HR-based: Log-logistic" = "dashed",
                          "HR-based: Lognormal" = "dashed")

HRbased_fill_list <- c("COU-AA-301: Original" = "black",
                       "COU-AA-301: Updated" = "grey40",
                       "External data (pooled)" = "grey80",
                        "HR-based: Weibull" = "white",
                        "HR-based: Gompertz" = "white",
                        "HR-based: Generalised gamma" = "white",
                        "HR-based: Exponential" = "white",
                        "HR-based: Log-logistic" = "white",
                        "HR-based: Lognormal" = "white")

KM_active_compare_HRbased <- KM_active_compare$plot +
  geom_line(data = times.expo.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.ggam.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.gomp.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.lnor.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.llog.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill = guide_legend(title="", nrow=3)) + 
  scale_colour_manual(values=HRbased_colours_list) +
  scale_linetype_manual(values=HRbased_lines_list)+
  scale_fill_manual(values=HRbased_fill_list)

KM_placebo_compare_HRbased <- KM_placebo_compare$plot +
  geom_line(data = times.expo.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.ggam.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.gomp.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.lnor.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.llog.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill = guide_legend(title="", nrow=3)) + 
  scale_colour_manual(values=HRbased_colours_list) +
  scale_linetype_manual(values=HRbased_lines_list)+
  scale_fill_manual(values=HRbased_fill_list)

KM_active_compareoldext_HRbased <- KM_active_oldexternal_ci_10$plot +
  geom_line(data = times.expo.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) +
  scale_colour_manual(values=HRbased_colours_list) +
  scale_linetype_manual(values=HRbased_lines_list)+
  scale_fill_manual(values=HRbased_fill_list)

KM_placebo_compareoldext_HRbased <- KM_placebo_oldexternal_ci_10$plot +
  geom_line(data = times.expo.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) + 
  geom_line(data = times.weib.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.ggam.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.gomp.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.lnor.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = times.llog.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour = guide_legend(title="", nrow=3), linetype = guide_legend(title="", nrow=3), fill= guide_legend(title="", nrow=3)) +
  scale_colour_manual(values=HRbased_colours_list) +
  scale_linetype_manual(values=HRbased_lines_list)+
  scale_fill_manual(values=HRbased_fill_list)

# Plot smoothed hazards versus HR-based models (rebased) # ----

curves_colours_list_hazards_HRbased <- c("COU-AA-301 (de Bono 2011) - AA, muhaz" = "black",
                                       "COU-AA-301 (de Bono 2011) - Placebo, muhaz" = "black",
                                       "HR-based: Exponential" = plot_red,
                                       "HR-based: Weibull" = plot_pink,
                                       "HR-based: Gompertz" = plot_forestgreen,
                                       "HR-based: Generalised gamma" = plot_bronze,
                                       "HR-based: Log-logistic" = plot_lightblue,
                                       "HR-based: Lognormal" = plot_purple)

curves_lines_list_hazards_HRbased <- c("COU-AA-301 (de Bono 2011) - AA, muhaz" = "solid",
                                     "COU-AA-301 (de Bono 2011) - Placebo, muhaz" = "solid",
                                     "HR-based: Exponential" = "dashed",
                                     "HR-based: Weibull" = "dashed",
                                     "HR-based: Gompertz" = "dashed",
                                     "HR-based: Generalised gamma" = "dashed",
                                     "HR-based: Log-logistic" = "dashed",
                                     "HR-based: Lognormal" = "dashed")

KM_active_compare_hazards_HRbased <- ggplot(hazards.muhaz.pivotal.original.active, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.weib.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.ggam.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.gomp.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.lnor.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.llog.HRbased.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards_HRbased) +
  scale_linetype_manual(values=curves_lines_list_hazards_HRbased)

KM_placebo_compare_hazards_HRbased <- ggplot(hazards.muhaz.pivotal.original.placebo, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.expo.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.weib.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.gomp.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.lnor.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.llog.HRbased.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_hazards_HRbased) +
  scale_linetype_manual(values=curves_lines_list_hazards_HRbased)

# Fitting Bayesian models (survHE) # ----

SD_prop_of_mean <- 0.05

# Active arm #

#Weibull - Gamma dist
#Mean: a/b
#SD: sqrt(a/b^2)
#a: target_mean*b
#thus, b: Mean/SD^2

input.hmc.weib <- exp(fit.weib.external.active$coefficients[1])

target_mean_weib <- as.numeric(input.hmc.weib)
target_sd_weib <- target_mean_weib * SD_prop_of_mean
#target_sd_weib <- as.numeric(fit.weib.external.active$res[1,4])
target_b_weib <- as.numeric(target_mean_weib/(target_sd_weib^2))
target_a_weib <- as.numeric(target_mean_weib * target_b_weib)

#Lognormal - uniform dist
#Mean: (a+b)/2
#SD: sqrt(((b-a)^2)/12)
#a: 2*Mean - b
#thus, b: 0.5*(sqrt(2*SD^2)+2*Mean)

input.hmc.lnor <- exp(fit.lnor.external.active$coefficients[2])

target_mean_lnor <- as.numeric(input.hmc.lnor)
target_sd_lnor <- target_mean_lnor * SD_prop_of_mean
target_b_lnor <- as.numeric(0.5*(sqrt(2*target_sd_lnor^2) + 2*target_mean_lnor))
target_a_lnor <- as.numeric(2*target_mean_lnor-target_b_lnor)

#Log-logistic - Gamma dist
#Mean: a/b
#SD: sqrt(a/b^2)
#a: target_mean*b
#thus, b: Mean/SD^2

input.hmc.llog <- exp(fit.llog.external.active$coefficients[1])

target_mean_llog <- as.numeric(input.hmc.llog)
target_sd_llog <- target_mean_llog * SD_prop_of_mean
target_b_llog <- as.numeric(target_mean_llog/(target_sd_llog^2))
target_a_llog <- as.numeric(target_mean_llog * target_b_llog)

#Gompertz - Gamma dist
#Mean: a/b
#SD: sqrt(a/b^2)
#a: target_mean*b
#thus, b: Mean/SD^2

input.hmc.gomp <- exp(fit.gomp.external.active$coefficients[1])
target_mean_gomp <- as.numeric(input.hmc.gomp)
target_sd_gomp <- target_mean_gomp * SD_prop_of_mean
target_b_gomp <- as.numeric(target_mean_gomp/(target_sd_gomp^2))
target_a_gomp <- as.numeric(target_mean_gomp * target_b_gomp)

#Gen gamma - two dists - 1 gamma, 1 normal
#gamma first

input.hmc.ggam <- exp(fit.ggam.external.active$coefficients[2])

target_mean_ggam <- as.numeric(input.hmc.ggam)
target_sd_ggam <- target_mean_ggam * SD_prop_of_mean
target_b_ggam <- as.numeric(target_mean_ggam/(target_sd_ggam^2))
target_a_ggam <- as.numeric(target_mean_ggam * target_b_ggam)

#then normal

input.hmc.ggam_2 <- fit.ggam.external.active$coefficients[3]
target_mu_Q_ggam <- input.hmc.ggam_2
target_sigma_Q_ggam <- target_mu_Q_ggam * SD_prop_of_mean

# all priors

hmc.priors <- list(weibull = list(a_alpha = target_a_weib,
                                  b_alpha = target_b_weib),
                   lognormal = list(a_alpha = target_a_lnor,
                                    b_alpha = target_b_lnor),
                   loglogistic = list(a_alpha = target_a_llog,
                                      b_alpha = target_b_llog),
                   gompertz = list(a_alpha = target_a_gomp,
                                   b_alpha = target_b_gomp),
                   gengamma = list(a_sigma = target_a_ggam,
                                   b_sigma = target_b_ggam,
                                   mu_Q = target_mu_Q_ggam,
                                   sigma_Q = target_sigma_Q_ggam))

hmc.weib.pivotal.original.active <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.active,
                                               dist="weibull", method = "hmc", priors=hmc.priors, seed=9999)

hmc.lnor.pivotal.original.active <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.active,
                                               dist="lognormal", method = "hmc", priors=hmc.priors, seed=9999)

hmc.llog.pivotal.original.active <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.active,
                                               dist="loglogistic", method = "hmc", priors=hmc.priors, seed=9999)

hmc.gomp.pivotal.original.active <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.active,
                                               dist="gompertz", method = "hmc", priors=hmc.priors, seed=9999) 

rownames(survival.data.pivotal.original.active) <- NULL #need this to avoid gen gamma code breaking, as relies on ID column

hmc.ggam.pivotal.original.active <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.active,
                                               dist="gengamma", method = "hmc", priors=hmc.priors, seed=9999)

prob.hmc.weib.pivotal.original.active <- make.surv(hmc.weib.pivotal.original.active, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.weib.pivotal.original.active <- as.data.frame(prob.hmc.weib.pivotal.original.active)
prob.hmc.weib.pivotal.original.active$Abi <- "Bayesian: Weibull"
prob.hmc.weib.pivotal.original.active$Abi <- as.factor(prob.hmc.weib.pivotal.original.active$Abi)

prob.hmc.lnor.pivotal.original.active <- make.surv(hmc.lnor.pivotal.original.active, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.lnor.pivotal.original.active <- as.data.frame(prob.hmc.lnor.pivotal.original.active)
prob.hmc.lnor.pivotal.original.active$Abi <- "Bayesian: Lognormal"
prob.hmc.lnor.pivotal.original.active$Abi <- as.factor(prob.hmc.lnor.pivotal.original.active$Abi)

prob.hmc.llog.pivotal.original.active <- make.surv(hmc.llog.pivotal.original.active, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.llog.pivotal.original.active <- as.data.frame(prob.hmc.llog.pivotal.original.active)
prob.hmc.llog.pivotal.original.active$Abi <- "Bayesian: Log-logistic"
prob.hmc.llog.pivotal.original.active$Abi <- as.factor(prob.hmc.llog.pivotal.original.active$Abi)

prob.hmc.gomp.pivotal.original.active <- make.surv(hmc.gomp.pivotal.original.active, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.gomp.pivotal.original.active <- as.data.frame(prob.hmc.gomp.pivotal.original.active)
prob.hmc.gomp.pivotal.original.active$Abi <- "Bayesian: Gompertz"
prob.hmc.gomp.pivotal.original.active$Abi <- as.factor(prob.hmc.gomp.pivotal.original.active$Abi)

prob.hmc.ggam.pivotal.original.active <- make.surv(hmc.ggam.pivotal.original.active, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.ggam.pivotal.original.active <- as.data.frame(prob.hmc.ggam.pivotal.original.active)
prob.hmc.ggam.pivotal.original.active$Abi <- "Bayesian: Generalised gamma"
prob.hmc.ggam.pivotal.original.active$Abi <- as.factor(prob.hmc.ggam.pivotal.original.active$Abi)

# Repeat for placebo arm #

#Weibull - Gamma dist
#Mean: a/b
#SD: sqrt(a/b^2)
#a: target_mean*b
#thus, b: Mean/SD^2

input.hmc.weib <- exp(fit.weib.external.placebo$coefficients[1])

target_mean_weib <- as.numeric(input.hmc.weib)
target_sd_weib <- target_mean_weib * SD_prop_of_mean
#target_sd_weib <- as.numeric(fit.weib.external.active$res[1,4])
target_b_weib <- as.numeric(target_mean_weib/(target_sd_weib^2))
target_a_weib <- as.numeric(target_mean_weib * target_b_weib)

#Lognormal - uniform dist
#Mean: (a+b)/2
#SD: sqrt(((b-a)^2)/12)
#a: 2*Mean - b
#thus, b: 0.5*(sqrt(2*SD^2)+2*Mean)

input.hmc.lnor <- exp(fit.lnor.external.placebo$coefficients[2])

target_mean_lnor <- as.numeric(input.hmc.lnor)
target_sd_lnor <- target_mean_lnor * SD_prop_of_mean
target_b_lnor <- as.numeric(0.5*(sqrt(2*target_sd_lnor^2) + 2*target_mean_lnor))
target_a_lnor <- as.numeric(2*target_mean_lnor-target_b_lnor)

#Log-logistic - Gamma dist
#Mean: a/b
#SD: sqrt(a/b^2)
#a: target_mean*b
#thus, b: Mean/SD^2

input.hmc.llog <- exp(fit.llog.external.placebo$coefficients[1])

target_mean_llog <- as.numeric(input.hmc.llog)
target_sd_llog <- target_mean_llog * SD_prop_of_mean
target_b_llog <- as.numeric(target_mean_llog/(target_sd_llog^2))
target_a_llog <- as.numeric(target_mean_llog * target_b_llog)

#Gompertz - Gamma dist
#Mean: a/b
#SD: sqrt(a/b^2)
#a: target_mean*b
#thus, b: Mean/SD^2

input.hmc.gomp <- exp(fit.gomp.external.placebo$coefficients[1])
target_mean_gomp <- as.numeric(input.hmc.gomp)
target_sd_gomp <- target_mean_gomp * SD_prop_of_mean
target_b_gomp <- as.numeric(target_mean_gomp/(target_sd_gomp^2))
target_a_gomp <- as.numeric(target_mean_gomp * target_b_gomp)

#Gen gamma - two dists - 1 gamma, 1 normal
#gamma first

input.hmc.ggam <- exp(fit.ggam.external.placebo$coefficients[2])

target_mean_ggam <- as.numeric(input.hmc.ggam)
target_sd_ggam <- target_mean_ggam * SD_prop_of_mean
target_b_ggam <- as.numeric(target_mean_ggam/(target_sd_ggam^2))
target_a_ggam <- as.numeric(target_mean_ggam * target_b_ggam)

#then normal

input.hmc.ggam_2 <- fit.ggam.external.placebo$coefficients[3]
target_mu_Q_ggam <- input.hmc.ggam_2
target_sigma_Q_ggam <- target_mu_Q_ggam * SD_prop_of_mean

# all priors

hmc.priors <- list(weibull = list(a_alpha = target_a_weib,
                                  b_alpha = target_b_weib),
                   lognormal = list(a_alpha = target_a_lnor,
                                    b_alpha = target_b_lnor),
                   loglogistic = list(a_alpha = target_a_llog,
                                      b_alpha = target_b_llog),
                   gompertz = list(a_alpha = target_a_gomp,
                                   b_alpha = target_b_gomp),
                   gengamma = list(a_sigma = target_a_ggam,
                                   b_sigma = target_b_ggam,
                                   mu_Q = target_mu_Q_ggam,
                                   sigma_Q = target_sigma_Q_ggam))

hmc.weib.pivotal.original.placebo <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.placebo,
                                               dist="weibull", method = "hmc", priors=hmc.priors, seed=9999)

hmc.lnor.pivotal.original.placebo <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.placebo,
                                               dist="lognormal", method = "hmc", priors=hmc.priors, seed=9999)

hmc.llog.pivotal.original.placebo <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.placebo,
                                               dist="loglogistic", method = "hmc", priors=hmc.priors, seed=9999)

hmc.gomp.pivotal.original.placebo <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.placebo,
                                               dist="gompertz", method = "hmc", priors=hmc.priors, seed=9999) 

rownames(survival.data.pivotal.original.placebo) <- NULL #need this to avoid gen gamma code breaking, as relies on ID column

hmc.ggam.pivotal.original.placebo <- fit.models(formula=Surv(Years, Event) ~ 1,
                                               data = survival.data.pivotal.original.placebo,
                                               dist="gengamma", method = "hmc", priors=hmc.priors, seed=9999)

prob.hmc.weib.pivotal.original.placebo <- make.surv(hmc.weib.pivotal.original.placebo, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.weib.pivotal.original.placebo <- as.data.frame(prob.hmc.weib.pivotal.original.placebo)
prob.hmc.weib.pivotal.original.placebo$Abi <- "Bayesian: Weibull"
prob.hmc.weib.pivotal.original.placebo$Abi <- as.factor(prob.hmc.weib.pivotal.original.placebo$Abi)

prob.hmc.lnor.pivotal.original.placebo <- make.surv(hmc.lnor.pivotal.original.placebo, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.lnor.pivotal.original.placebo <- as.data.frame(prob.hmc.lnor.pivotal.original.placebo)
prob.hmc.lnor.pivotal.original.placebo$Abi <- "Bayesian: Lognormal"
prob.hmc.lnor.pivotal.original.placebo$Abi <- as.factor(prob.hmc.lnor.pivotal.original.placebo$Abi)

prob.hmc.llog.pivotal.original.placebo <- make.surv(hmc.llog.pivotal.original.placebo, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.llog.pivotal.original.placebo <- as.data.frame(prob.hmc.llog.pivotal.original.placebo)
prob.hmc.llog.pivotal.original.placebo$Abi <- "Bayesian: Log-logistic"
prob.hmc.llog.pivotal.original.placebo$Abi <- as.factor(prob.hmc.llog.pivotal.original.placebo$Abi)

prob.hmc.gomp.pivotal.original.placebo <- make.surv(hmc.gomp.pivotal.original.placebo, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.gomp.pivotal.original.placebo <- as.data.frame(prob.hmc.gomp.pivotal.original.placebo)
prob.hmc.gomp.pivotal.original.placebo$Abi <- "Bayesian: Gompertz"
prob.hmc.gomp.pivotal.original.placebo$Abi <- as.factor(prob.hmc.gomp.pivotal.original.placebo$Abi)

prob.hmc.ggam.pivotal.original.placebo <- make.surv(hmc.ggam.pivotal.original.placebo, t = seq(lower_limit, upper_limit, increment)) 
prob.hmc.ggam.pivotal.original.placebo <- as.data.frame(prob.hmc.ggam.pivotal.original.placebo)
prob.hmc.ggam.pivotal.original.placebo$Abi <- "Bayesian: Generalised gamma"
prob.hmc.ggam.pivotal.original.placebo$Abi <- as.factor(prob.hmc.ggam.pivotal.original.placebo$Abi)

# Goodness of fit #

AICandBIC.pivotal.original.bayesian <- cbind(rbind(hmc.weib.pivotal.original.active$model.fitting$aic,
                                                   hmc.ggam.pivotal.original.active$model.fitting$aic,
                                                   hmc.gomp.pivotal.original.active$model.fitting$aic,
                                                   hmc.lnor.pivotal.original.active$model.fitting$aic,
                                                   hmc.llog.pivotal.original.active$model.fitting$aic),
                                             rbind(hmc.weib.pivotal.original.active$model.fitting$bic,
                                                   hmc.ggam.pivotal.original.active$model.fitting$bic,
                                                   hmc.gomp.pivotal.original.active$model.fitting$bic,
                                                   hmc.lnor.pivotal.original.active$model.fitting$bic,
                                                   hmc.llog.pivotal.original.active$model.fitting$bic),
                                             rbind(hmc.weib.pivotal.original.active$model.fitting$dic,
                                                   hmc.ggam.pivotal.original.active$model.fitting$dic,
                                                   hmc.gomp.pivotal.original.active$model.fitting$dic,
                                                   hmc.lnor.pivotal.original.active$model.fitting$dic,
                                                   hmc.llog.pivotal.original.active$model.fitting$dic),
                                             rbind(hmc.weib.pivotal.original.placebo$model.fitting$aic,
                                                   hmc.ggam.pivotal.original.placebo$model.fitting$aic,
                                                   hmc.gomp.pivotal.original.placebo$model.fitting$aic,
                                                   hmc.lnor.pivotal.original.placebo$model.fitting$aic,
                                                   hmc.llog.pivotal.original.placebo$model.fitting$aic),
                                             rbind(hmc.weib.pivotal.original.placebo$model.fitting$bic,
                                                   hmc.ggam.pivotal.original.placebo$model.fitting$bic,
                                                   hmc.gomp.pivotal.original.placebo$model.fitting$bic,
                                                   hmc.lnor.pivotal.original.placebo$model.fitting$bic,
                                                   hmc.llog.pivotal.original.placebo$model.fitting$bic),
                                             rbind(hmc.weib.pivotal.original.placebo$model.fitting$dic,
                                                   hmc.ggam.pivotal.original.placebo$model.fitting$dic,
                                                   hmc.gomp.pivotal.original.placebo$model.fitting$dic,
                                                   hmc.lnor.pivotal.original.placebo$model.fitting$dic,
                                                   hmc.llog.pivotal.original.placebo$model.fitting$dic))

rownames(AICandBIC.pivotal.original.bayesian) <- rbind("weib","ggam","gomp","lnor","llog")
colnames(AICandBIC.pivotal.original.bayesian) <- rbind("active, AIC","active, BIC","active, DIC","placebo, AIC","placebo, BIC","placebo, DIC")

# Putting in same format as models fitted in flexsurv #

times.weib.bayesian.active <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.weib.pivotal.original.active$S.S, Abi = prob.hmc.weib.pivotal.original.active$Abi)
times.ggam.bayesian.active <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.ggam.pivotal.original.active$S.S, Abi = prob.hmc.ggam.pivotal.original.active$Abi)
times.gomp.bayesian.active <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.gomp.pivotal.original.active$S.S, Abi = prob.hmc.gomp.pivotal.original.active$Abi)
times.lnor.bayesian.active <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.lnor.pivotal.original.active$S.S, Abi = prob.hmc.lnor.pivotal.original.active$Abi)
times.llog.bayesian.active <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.llog.pivotal.original.active$S.S, Abi = prob.hmc.llog.pivotal.original.active$Abi)

times.weib.bayesian.placebo <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.weib.pivotal.original.placebo$S.S, Abi = prob.hmc.weib.pivotal.original.placebo$Abi)
times.ggam.bayesian.placebo <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.ggam.pivotal.original.placebo$S.S, Abi = prob.hmc.ggam.pivotal.original.placebo$Abi)
times.gomp.bayesian.placebo <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.gomp.pivotal.original.placebo$S.S, Abi = prob.hmc.gomp.pivotal.original.placebo$Abi)
times.lnor.bayesian.placebo <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.lnor.pivotal.original.placebo$S.S, Abi = prob.hmc.lnor.pivotal.original.placebo$Abi)
times.llog.bayesian.placebo <- data.frame(time = seq(lower_limit, upper_limit, increment), est = prob.hmc.llog.pivotal.original.placebo$S.S, Abi = prob.hmc.llog.pivotal.original.placebo$Abi)

# Plot Bayesian models (survHE) # ----


curves_colours_list_Bayesian <- c("COU-AA-301: Original" = "black",
                                  "COU-AA-301: Updated" = "grey40",
                                  "External data (pooled)" = "grey80",
                                  "Bayesian: Weibull" = plot_pink,
                                  "Bayesian: Lognormal" = plot_purple,
                                  "Bayesian: Log-logistic" = plot_lightblue,
                                  "Bayesian: Gompertz" = plot_forestgreen,
                                  "Bayesian: Generalised gamma" = plot_bronze,
                                  "COU-AA-301 (de Bono 2011) - AA, muhaz" = "black",
                                  "COU-AA-301 (de Bono 2011) - Placebo, muhaz" = "black")

curves_lines_list_Bayesian <- c("COU-AA-301: Original" = "solid",
                                "COU-AA-301: Updated" = "solid",
                                "External data (pooled)" = "solid",
                                "Bayesian: Weibull" = "dashed",
                                "Bayesian: Lognormal" = "dashed",
                                "Bayesian: Log-logistic" = "dashed",
                                "Bayesian: Gompertz" = "dashed",
                                "Bayesian: Generalised gamma" = "dashed",
                                "COU-AA-301 (de Bono 2011) - AA, muhaz" = "solid",
                                "COU-AA-301 (de Bono 2011) - Placebo, muhaz" = "solid")

curves_fill_list_Bayesian <- c("COU-AA-301: Original" = "black",
                               "COU-AA-301: Updated" = "grey40",
                               "External data (pooled)" = "grey80",
                               "Bayesian: Weibull" = "white",
                               "Bayesian: Lognormal" = "white",
                               "Bayesian: Log-logistic" = "white",
                               "Bayesian: Gompertz" = "white",
                               "Bayesian: Generalised gamma" = "white")

legend_order <- c("COU-AA-301: Original",
                  "External data (pooled)",
                  "Bayesian: Generalised gamma",
                  "Bayesian: Gompertz",
                  "Bayesian: Log-logistic",
                  "Bayesian: Lognormal",
                  "Bayesian: Weibull")

legend_order_full <- c("COU-AA-301: Original",
                       "COU-AA-301: Updated",
                       "External data (pooled)",
                       "Bayesian: Generalised gamma",
                       "Bayesian: Gompertz",
                       "Bayesian: Log-logistic",
                       "Bayesian: Lognormal",
                       "Bayesian: Weibull")

KM_active_compareoldext_bayesian <- KM_active_oldexternal_ci_10$plot +
  geom_line(data = prob.hmc.weib.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.lnor.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.llog.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.gomp.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.ggam.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill = guide_legend(title="", nrow=3)) + 
  scale_linetype_manual(values=curves_lines_list_Bayesian, limits = legend_order) +
  scale_colour_manual(values=curves_colours_list_Bayesian, limits = legend_order) +
  scale_fill_manual(values=curves_fill_list_Bayesian, limits = legend_order)

KM_placebo_compareoldext_bayesian <- KM_placebo_oldexternal_ci_10$plot +
  geom_line(data = prob.hmc.weib.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.lnor.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.llog.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.gomp.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.ggam.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill = guide_legend(title="", nrow=3)) + 
  scale_linetype_manual(values=curves_lines_list_Bayesian, limits = legend_order) +
  scale_colour_manual(values=curves_colours_list_Bayesian, limits = legend_order) +
  scale_fill_manual(values=curves_fill_list_Bayesian, limits = legend_order)

KM_active_compare_bayesian <- KM_active_compare$plot +
  geom_line(data = prob.hmc.weib.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.lnor.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.llog.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.gomp.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.ggam.pivotal.original.active, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill = guide_legend(title="", nrow=3)) + 
  scale_linetype_manual(values=curves_lines_list_Bayesian, limits = legend_order_full) +
  scale_colour_manual(values=curves_colours_list_Bayesian, limits = legend_order_full) +
  scale_fill_manual(values=curves_fill_list_Bayesian, limits = legend_order_full)

KM_placebo_compare_bayesian <- KM_placebo_compare$plot +
  geom_line(data = prob.hmc.weib.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.lnor.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.llog.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.gomp.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  geom_line(data = prob.hmc.ggam.pivotal.original.placebo, aes(x = S.t, y = S.S, colour = Abi, linetype = Abi, fill = Abi), size=1) +
  guides(colour=guide_legend(title="", nrow=3),linetype = guide_legend(title="", nrow=3), fill = guide_legend(title="", nrow=3)) + 
  scale_linetype_manual(values=curves_lines_list_Bayesian, limits = legend_order_full) +
  scale_colour_manual(values=curves_colours_list_Bayesian, limits = legend_order_full) +
  scale_fill_manual(values=curves_fill_list_Bayesian, limits = legend_order_full)

# Plot smoothed hazards versus Bayesian models #----

#active arm #

#Weibull
weib_active_lambda_aft <- mean(as.data.frame(hmc.weib.pivotal.original.active$models$`Weibull (AFT)`)[,4]) #scale
weib_active_gamma <- mean(as.data.frame(hmc.weib.pivotal.original.active$models$`Weibull (AFT)`)[,3])  #shape
weib_active_lambda_ph <- weib_active_lambda_aft^(-weib_active_gamma) #shape_ph
hazards.weib.bayesian.active <- data.frame(time = hazards.weib.pivotal.original.active$time)
hazards.weib.bayesian.active$est <- weib_active_lambda_ph*weib_active_gamma*(hazards.weib.bayesian.active$time^(weib_active_gamma-1)) #hazard function
hazards.weib.bayesian.active$Abi <- "Bayesian: Weibull"

#lognormal
lnor_active_meanlog <- mean(as.data.frame(hmc.lnor.pivotal.original.active$models$`log-Normal`)[,4]) #meanlog
lnor_active_sdlog <- mean(as.data.frame(hmc.lnor.pivotal.original.active$models$`log-Normal`)[,3])  #sdlog
hazards.lnor.bayesian.active <- data.frame(time = hazards.lnor.pivotal.original.active$time)
hazards.lnor.bayesian.active$est <-(1/(lnor_active_sdlog*sqrt(2*pi)))*(1/hazards.lnor.bayesian.active$time)*exp(-(((log(hazards.lnor.bayesian.active$time)-lnor_active_meanlog)^2)/(2*(lnor_active_sdlog^2))))
hazards.lnor.bayesian.active$est <- hazards.lnor.bayesian.active$est/prob.hmc.lnor.pivotal.original.active$S.S
hazards.lnor.bayesian.active$Abi <- "Bayesian: Lognormal"

#loglogistic
llog_active_scale <- mean(as.data.frame(hmc.llog.pivotal.original.active$models$`log-Logistic`)[,4]) #scale
llog_active_shape <- mean(as.data.frame(hmc.llog.pivotal.original.active$models$`log-Logistic`)[,3])  #shape
hazards.llog.bayesian.active <- data.frame(time = hazards.llog.pivotal.original.active$time)
hazards.llog.bayesian.active$est <- ((llog_active_shape/llog_active_scale)*(hazards.llog.bayesian.active$time/llog_active_scale)^(llog_active_shape-1))/(1+(hazards.llog.bayesian.active$time/llog_active_scale)^llog_active_shape) #hazard function
hazards.llog.bayesian.active$Abi <- "Bayesian: Log-logistic"

#generalised gamma
ggam_active_mu <- mean(as.data.frame(hmc.ggam.pivotal.original.active$models$`Gen. Gamma`)[,3])
ggam_active_sigma <- mean(as.data.frame(hmc.ggam.pivotal.original.active$models$`Gen. Gamma`)[,2])
ggam_active_Q <- mean(as.data.frame(hmc.ggam.pivotal.original.active$models$`Gen. Gamma`)[,1])
hazards.ggam.bayesian.active <- data.frame(time = hazards.ggam.pivotal.original.active$time)
hazards.ggam.bayesian.active$est <- hgengamma(x = hazards.ggam.bayesian.active$time, mu=ggam_active_mu, sigma=ggam_active_sigma, Q=ggam_active_Q)
hazards.ggam.bayesian.active$Abi <- "Bayesian: Generalised gamma"

#Gompertz
#h(t) = scale*exp(shape*t)
gomp_active_scale <- mean(as.data.frame(hmc.gomp.pivotal.original.active$models$Gompertz)[,4]) #scale
gomp_active_shape <- mean(as.data.frame(hmc.gomp.pivotal.original.active$models$Gompertz)[,3])  #shape
hazards.gomp.bayesian.active <- data.frame(time = hazards.gomp.pivotal.original.active$time)
hazards.gomp.bayesian.active$est <- hgompertz(x = hazards.gomp.bayesian.active$time, shape = gomp_active_shape, rate = gomp_active_scale, log = FALSE)
hazards.gomp.bayesian.active$Abi <- "Bayesian: Gompertz"

#placebo arm #

#Weibull
weib_placebo_lambda_aft <- mean(as.data.frame(hmc.weib.pivotal.original.placebo$models$`Weibull (AFT)`)[,4]) #scale
weib_placebo_gamma <- mean(as.data.frame(hmc.weib.pivotal.original.placebo$models$`Weibull (AFT)`)[,3])  #shape
weib_placebo_lambda_ph <- weib_placebo_lambda_aft^(-weib_placebo_gamma) #shape_ph
hazards.weib.bayesian.placebo <- data.frame(time = hazards.weib.pivotal.original.placebo$time)
hazards.weib.bayesian.placebo$est <- weib_placebo_lambda_ph*weib_placebo_gamma*(hazards.weib.bayesian.placebo$time^(weib_placebo_gamma-1)) #hazard function
hazards.weib.bayesian.placebo$Abi <- "Bayesian: Weibull"

#lognormal
lnor_placebo_meanlog <- mean(as.data.frame(hmc.lnor.pivotal.original.placebo$models$`log-Normal`)[,4]) #meanlog
lnor_placebo_sdlog <- mean(as.data.frame(hmc.lnor.pivotal.original.placebo$models$`log-Normal`)[,3])  #sdlog
hazards.lnor.bayesian.placebo <- data.frame(time = hazards.lnor.pivotal.original.placebo$time)
hazards.lnor.bayesian.placebo$est <-(1/(lnor_placebo_sdlog*sqrt(2*pi)))*(1/hazards.lnor.bayesian.placebo$time)*exp(-(((log(hazards.lnor.bayesian.placebo$time)-lnor_placebo_meanlog)^2)/(2*(lnor_placebo_sdlog^2))))
hazards.lnor.bayesian.placebo$est <- hazards.lnor.bayesian.placebo$est/prob.hmc.lnor.pivotal.original.placebo$S.S
hazards.lnor.bayesian.placebo$Abi <- "Bayesian: Lognormal"

#loglogistic
#h(t) =  ((shape/scale)×(time/scale)^(shape-1))/(1+(time/scale)^shape)

llog_placebo_scale <- mean(as.data.frame(hmc.llog.pivotal.original.placebo$models$`log-Logistic`)[,4]) #scale
llog_placebo_shape <- mean(as.data.frame(hmc.llog.pivotal.original.placebo$models$`log-Logistic`)[,3])  #shape
hazards.llog.bayesian.placebo <- data.frame(time = hazards.llog.pivotal.original.placebo$time)
hazards.llog.bayesian.placebo$est <- ((llog_placebo_shape/llog_placebo_scale)*(hazards.llog.bayesian.placebo$time/llog_placebo_scale)^(llog_placebo_shape-1))/(1+(hazards.llog.bayesian.placebo$time/llog_placebo_scale)^llog_placebo_shape) #hazard function
hazards.llog.bayesian.placebo$Abi <- "Bayesian: Log-logistic"

#generalised gamma
ggam_placebo_mu <- mean(as.data.frame(hmc.ggam.pivotal.original.placebo$models$`Gen. Gamma`)[,3])
ggam_placebo_sigma <- mean(as.data.frame(hmc.ggam.pivotal.original.placebo$models$`Gen. Gamma`)[,2])
ggam_placebo_Q <- mean(as.data.frame(hmc.ggam.pivotal.original.placebo$models$`Gen. Gamma`)[,1])
hazards.ggam.bayesian.placebo <- data.frame(time = hazards.ggam.pivotal.original.placebo$time)
hazards.ggam.bayesian.placebo$est <- hgengamma(x = hazards.ggam.bayesian.placebo$time, mu=ggam_placebo_mu, sigma=ggam_placebo_sigma, Q=ggam_placebo_Q)
hazards.ggam.bayesian.placebo$Abi <- "Bayesian: Generalised gamma"

#Gompertz
#h(t) = scale*exp(shape*t)
gomp_placebo_scale <- mean(as.data.frame(hmc.gomp.pivotal.original.placebo$models$Gompertz)[,4]) #scale
gomp_placebo_shape <- mean(as.data.frame(hmc.gomp.pivotal.original.placebo$models$Gompertz)[,3])  #shape
hazards.gomp.bayesian.placebo <- data.frame(time = hazards.gomp.pivotal.original.placebo$time)
hazards.gomp.bayesian.placebo$est <- hgompertz(x = hazards.gomp.bayesian.placebo$time, shape = gomp_placebo_shape, rate = gomp_placebo_scale, log = FALSE)
hazards.gomp.bayesian.placebo$Abi <- "Bayesian: Gompertz"

legend_order_muhaz_active <- c("COU-AA-301 (de Bono 2011) - AA, muhaz",
                        "Bayesian: Generalised gamma",
                        "Bayesian: Gompertz",
                        "Bayesian: Log-logistic",
                        "Bayesian: Lognormal",
                        "Bayesian: Weibull")

KM_active_compare_hazards_bayesian <- ggplot(hazards.muhaz.pivotal.original.active, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.weib.bayesian.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.ggam.bayesian.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.gomp.bayesian.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.lnor.bayesian.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.llog.bayesian.active, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_Bayesian, limits = legend_order_muhaz_active) +
  scale_linetype_manual(values=curves_lines_list_Bayesian, limits = legend_order_muhaz_active)

legend_order_muhaz_placebo <- c("COU-AA-301 (de Bono 2011) - Placebo, muhaz",
                               "Bayesian: Generalised gamma",
                               "Bayesian: Gompertz",
                               "Bayesian: Log-logistic",
                               "Bayesian: Lognormal",
                               "Bayesian: Weibull")


KM_placebo_compare_hazards_bayesian <- ggplot(hazards.muhaz.pivotal.original.placebo, aes(x=time, y=est, colour= Abi, linetype = Abi)) +
  geom_line(size=1) +
  geom_line(data = hazards.weib.bayesian.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  geom_line(data = hazards.ggam.bayesian.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.gomp.bayesian.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.lnor.bayesian.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) +
  geom_line(data = hazards.llog.bayesian.placebo, aes(x = time, y = est, colour = Abi, linetype = Abi), size = 1) + 
  theme_classic() +
  theme(legend.position = "bottom") + 
  scale_y_continuous(limits = c(0,2), oob = rescale_none, breaks=seq(0,2,0.25)) +
  scale_x_continuous(limits = c(0,10), oob = rescale_none, breaks=seq(0,10,1)) +
  xlab("Time (years)") + 
  ylab("Hazard") + 
  theme(axis.title.y = element_text(size = 10, angle = 90, face = "bold")) + 
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "bold")) + 
  theme(axis.text.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(size = 10)) + 
  theme(legend.text = element_text(size = 8)) +
  guides(colour=guide_legend(title="", nrow=4),linetype = guide_legend(title="", nrow=4)) + 
  scale_colour_manual(values=curves_colours_list_Bayesian, limits = legend_order_muhaz_placebo) +
  scale_linetype_manual(values=curves_lines_list_Bayesian, limits = legend_order_muhaz_placebo)

# Extract survival times # ----

data_model_basic_active <- cbind(times.expo.pivotal.original.active,
                                 times.weib.pivotal.original.active,
                                 times.ggam.pivotal.original.active,
                                 times.gomp.pivotal.original.active,
                                 times.lnor.pivotal.original.active,
                                 times.llog.pivotal.original.active)
     
data_model_basic_placebo <- cbind(times.expo.pivotal.original.placebo,
                                  times.weib.pivotal.original.placebo,
                                  times.ggam.pivotal.original.placebo,
                                  times.gomp.pivotal.original.placebo,
                                  times.lnor.pivotal.original.placebo,
                                  times.llog.pivotal.original.placebo)

data_model_piecewise_active <- cbind(times.expo.piecewise.active,
                                 times.weib.piecewise.active,
                                 times.ggam.piecewise.active,
                                 times.gomp.piecewise.active,
                                 times.lnor.piecewise.active,
                                 times.llog.piecewise.active)

data_model_piecewise_placebo <- cbind(times.expo.piecewise.placebo,
                                  times.weib.piecewise.placebo,
                                  times.ggam.piecewise.placebo,
                                  times.gomp.piecewise.placebo,
                                  times.lnor.piecewise.placebo,
                                  times.llog.piecewise.placebo)

data_model_HRbased_active <- cbind(times.expo.HRbased.active,
                                     times.weib.HRbased.active,
                                     times.ggam.HRbased.active,
                                     times.gomp.HRbased.active,
                                     times.lnor.HRbased.active,
                                     times.llog.HRbased.active)

data_model_HRbased_placebo <- cbind(times.expo.HRbased.placebo,
                                      times.weib.HRbased.placebo,
                                      times.ggam.HRbased.placebo,
                                      times.gomp.HRbased.placebo,
                                      times.lnor.HRbased.placebo,
                                      times.llog.HRbased.placebo)

data_model_bayesian_active <- cbind(times.weib.bayesian.active,
                                 times.ggam.bayesian.active,
                                 times.gomp.bayesian.active,
                                 times.lnor.bayesian.active,
                                 times.llog.bayesian.active)

data_model_bayesian_placebo <- cbind(times.weib.bayesian.placebo,
                                  times.ggam.bayesian.placebo,
                                  times.gomp.bayesian.placebo,
                                  times.lnor.bayesian.placebo,
                                  times.llog.bayesian.placebo)

write.csv(data_model_basic_active, "active, basic.csv")
write.csv(data_model_basic_placebo, "placebo, basic.csv")
write.csv(data_model_piecewise_active, "active, piecewise.csv")
write.csv(data_model_piecewise_placebo, "placebo, piecewise.csv")
write.csv(data_model_HRbased_active, "active, HRbased.csv")
write.csv(data_model_HRbased_placebo, "placebo, HRbased.csv")
write.csv(data_model_bayesian_active, "active, bayesian.csv")
write.csv(data_model_bayesian_placebo, "placebo, bayesian.csv")

summary(KM.est.active.pivotal.original, times = c(1,12.8/12,1.25,1.5,20.2/12,1.75,2))
summary(KM.est.placebo.pivotal.original, times = c(1,12.8/12,1.25,1.5,20.2/12,1.75,2))

summary(KM.est.active.pivotal.updated, times = c(1,12.8/12,1.25,1.5,20.2/12,1.75,2))
summary(KM.est.placebo.pivotal.updated, times = c(1,12.8/12,1.25,1.5,20.2/12,1.75,2))

# RMST estimates (KM) ----

# Note: the rmsth function does not account for the end of follow-up, so  enter timepoints after follow-up with caution!

rmst_current <- survival.data.pivotal.original.active
rmst_times = c(1,12.8/12,1.25,1.5,20.2/12,1.75,2)
rmst_one <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[1], eps = 1.0e-08)
rmst_two <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[2], eps = 1.0e-08)
rmst_thr <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[3], eps = 1.0e-08)
rmst_fou <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[4], eps = 1.0e-08)
rmst_fiv <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[5], eps = 1.0e-08)
rmst_six <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[6], eps = 1.0e-08)
rmst_sev <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[7], eps = 1.0e-08)

rmst_one$rmst
rmst_two$rmst
rmst_thr$rmst
rmst_fou$rmst
rmst_fiv$rmst
rmst_six$rmst
rmst_sev$rmst

rmst_current <- survival.data.pivotal.updated.active
rmst_times = c(1,12.8/12,1.25,1.5,20.2/12,1.75,2)
rmst_one <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[1], eps = 1.0e-08)
rmst_two <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[2], eps = 1.0e-08)
rmst_thr <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[3], eps = 1.0e-08)
rmst_fou <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[4], eps = 1.0e-08)
rmst_fiv <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[5], eps = 1.0e-08)
rmst_six <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[6], eps = 1.0e-08)
rmst_sev <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[7], eps = 1.0e-08)

rmst_one$rmst
rmst_two$rmst
rmst_thr$rmst
rmst_fou$rmst
rmst_fiv$rmst
rmst_six$rmst
rmst_sev$rmst

rmst_current <- survival.data.pivotal.original.placebo
rmst_times = c(1,12.8/12,1.25,1.5,20.2/12,1.75,2)
rmst_one <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[1], eps = 1.0e-08)
rmst_two <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[2], eps = 1.0e-08)
rmst_thr <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[3], eps = 1.0e-08)
rmst_fou <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[4], eps = 1.0e-08)
rmst_fiv <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[5], eps = 1.0e-08)
rmst_six <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[6], eps = 1.0e-08)
rmst_sev <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[7], eps = 1.0e-08)

rmst_one$rmst
rmst_two$rmst
rmst_thr$rmst
rmst_fou$rmst
rmst_fiv$rmst
rmst_six$rmst
rmst_sev$rmst

rmst_current <- survival.data.pivotal.updated.placebo
rmst_times = c(1,12.8/12,1.25,1.5,20.2/12,1.75,2)
rmst_one <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[1], eps = 1.0e-08)
rmst_two <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[2], eps = 1.0e-08)
rmst_thr <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[3], eps = 1.0e-08)
rmst_fou <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[4], eps = 1.0e-08)
rmst_fiv <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[5], eps = 1.0e-08)
rmst_six <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[6], eps = 1.0e-08)
rmst_sev <- rmsth(y = rmst_current$Years, d = rmst_current$Event, tcut = rmst_times[7], eps = 1.0e-08)

rmst_one$rmst
rmst_two$rmst
rmst_thr$rmst
rmst_fou$rmst
rmst_fiv$rmst
rmst_six$rmst
rmst_sev$rmst

#### ----

# Plots included in report # ----

adjustment <- scale_y_continuous(
  labels = scales::number_format(scale = 100,
                                 accuracy = 1,
                                 decimal.mark = '.',
                                 suffix = "%"), breaks = c(seq(0,1,0.1)))
KM_original_AA #Figure 2

KM_active #Figure 3
KM_placebo #Figure 3

muhaz_active #Figure 3
muhaz_placebo #Figure 3
muhaz_active_rebased #SA Figure 1
muhaz_placebo_rebased #SA Figure 1

muhaz_pehaz_plot #Figure 2

muhaz_active_pooled
muhaz_active_pooled_rebased
muhaz_placebo_pooled
muhaz_placebo_pooled_rebased

KM_active_compare
KM_placebo_compare

KM_active_newold_ci
KM_placebo_newold_ci

KM_active_compare_basicmodels 
KM_placebo_compare_basicmodels

KM_active_external_basicmodels
KM_placebo_external_basicmodels
KM_active_external_rebased_basicmodels + adjustment #SA Figure 3
KM_placebo_external_rebased_basicmodels + adjustment #SA Figure 3

KM_active_original_compare_basicmodels
KM_placebo_original_compare_basicmodels
KM_active_compareoldext_basicmodels + adjustment #SA Figure 2
KM_placebo_compareoldext_basicmodels + adjustment #SA Figure 2
KM_active_compare_hazards #SA Figure 2
KM_placebo_compare_hazards #SA Figure 2

KM_active_compare_basicmodels_reduced
KM_placebo_compare_basicmodels_reduced
KM_active_compare_basicmodels_hazards_reduced
KM_placebo_compare_basicmodels_hazards_reduced

KM_active_compare_piecewisemodels
KM_placebo_compare_piecewisemodels
KM_active_compareoldext_piecewisemodels + adjustment #SA Figure 4
KM_placebo_compareoldext_piecewisemodels + adjustment #SA Figure 4
KM_active_compare_hazards_piecewise #SA Figure 4
KM_placebo_compare_hazards_piecewise #SA Figure 4

PH.test.active
PH.test.placebo
HR.active
HR.placebo

KM_active_compare_HRbased
KM_placebo_compare_HRbased
KM_active_compareoldext_HRbased + adjustment #SA Figure 5
KM_placebo_compareoldext_HRbased + adjustment #SA Figure 5
KM_active_compare_hazards_HRbased #SA Figure 5
KM_placebo_compare_hazards_HRbased #SA Figure 5

KM_active_compare_bayesian
KM_placebo_compare_bayesian
KM_active_compareoldext_bayesian + adjustment #SA Figure 6
KM_placebo_compareoldext_bayesian + adjustment #SA Figure 6 
KM_active_compare_hazards_bayesian #SA Figure 6
KM_placebo_compare_hazards_bayesian #SA Figure 6

AICandBIC.pivotal.original
AICandBIC.external.rebased
AICandBIC.pivotal.original.bayesian #also includes DIC

# Confidence intervals for the Cox PH model # ----

PH.test.active
PH.test.placebo

exp(confint(PH.test.active))
exp(confint(PH.test.placebo))
