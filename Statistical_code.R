#--------------------------------
#libraries
#--------------------------------
library(haven)
library(lavaan)
library(knitr)
library(semPlot)
library(dplyr)
library(equaltestMI)
library(summarytools)
library(readxl)

url <- "https://raw.githubusercontent.com/haaszevelin/KPT/main/KPT_2022.xlsx"
download.file(url, destfile = "KPT_2022.xlsx", mode = "wb")

data <- read_excel("KPT_2022.xlsx")

data$gender <- factor(data$gender,
                     levels = c(1,2),
                     labels = c("Men", "Women"))

table(data$gender)

data <- subset(data,
                 !is.na(gender))
table(data$gender)

colnames(data)

#-----------------------
# CFA
#-----------------------

overall.model <- 'Gf =~ GFF_Raven + GFF_Figures; Ga =~ GAF_Firstsounds + GAF_Lastsounds; 
PSMF =~ PSMF_Bodyparts + PSMF_Fingers + PSMF_Hands; Gwm =~ WMF_FDS + WMF_BDS
+ WMF_Corsi; Balance_Open =~ BF_BOR + BF_BOL; Balance_Closed =~ BF_BCR + BF_BCL '


overall.fit <- cfa(model = overall.model,
                   data = data,
                   meanstructure = TRUE)

summary(overall.fit,
        standardized = TRUE,
        rsquare = TRUE,
        fit.measure = TRUE)

table_fit <- matrix(NA, nrow=7, ncol=9)
colnames(table_fit) = c("Model", "X2", "df", "CFI","TLI", "RMSEA", "SRMR", "AIC", "BIC")
table_fit[1, ] <- c("Overall Model", round(fitmeasures(overall.fit,
                                                       c("chisq", "df", "cfi", "tli",
                                                         "rmsea", "srmr", "aic", "bic")),3))

kable(table_fit)

#------------------------------
#plot
#-----------------------------

semPaths(overall.fit,
         what = "std",
         whatLabels = "std",
         style = "ram",
         layout = "circle",
         residuals = FALSE,
         intercepts = FALSE,
         sizeLat = 6,
         sizeMan = 7,
         edge.label.cex = 1,
         edge.color = "black",
         edge.width = 0.5,
         fade = FALSE,
         asize = 2,
         esize =5,
         label.cex = 1.2,
         nCharNodes = 0,
         curvePivot = TRUE,
         mar = c(2,2,2,2),
         color = list(lat = "#B3D1F0", man = "#E6F0FA", arrow = "#08306B"),
         nodeLabels = c(
           "GFF_Raven" = "Raven",
           "GFF_Figurak" = "Figures",
           "GAF_NY1" = "First sound",
           "GAF_NY2" = "Last sound",
           "PSMF_Testreszek" = "Body parts",
           "PSMF_Ujjak" = "Fingers",
           "PSMF_Kezek"= "Hands",
           "WMF_FDS" = "FDS",
           "WMF_BDS" = "BDS",
           "WMF_Corsi" = "Corsi",
           "BF_NYJ" = "Right",
           "BF_NYB" = "Left",
           "BF_CSJ" = "Right",
           "BF_CSB" = "Left",
           "Fluid" = "Gf",
           "Phono" = "Ga",
           "Senso" = "PSMF",
           "Wmemory" = "Gwm",
           "Balance_Open" = "BOF",
           "Balance_Closed" = "BCF")
)


summary(overall.fit, standardized = TRUE)
inspect(overall.fit, "cor.lv")

#---------------------------------------------------------
#bootstrap for standardized ci-s (for the factor loadings)
#---------------------------------------------------------
overall.fit <- cfa(
  model = overall.model,  
  data = data,
  std.lv = TRUE,           
  se = "bootstrap",        
  bootstrap = 1000         
)

loadings_ci <- standardizedSolution(overall.fit, type = "std.all", ci = TRUE) %>%
  filter(op == "=~") %>%   
  select(lhs, rhs, est.std, ci.lower, ci.upper) %>%
  rename(
    Factor = lhs,
    Variable = rhs,
    Loading = est.std,
    CI_lower = ci.lower,
    CI_upper = ci.upper
  )

print(loadings_ci)

#-------------------------------------------------
# CI for latent factor covariances
#-------------------------------------------------

factors <- c("Gf", "PSMF", "Ga", "Gwm", "BCF", "BOF")

pe <- parameterEstimates(
  overall.fit,
  boot.ci.type = "perc",
  level = 0.95
)

latent_cov_ci <- pe %>%
  dplyr::filter(
    op == "~~",
    lhs %in% factors,
    rhs %in% factors,
    lhs != rhs
  ) %>%
  rowwise() %>%
  mutate(pair = paste0(sort(c(lhs, rhs)), collapse = "~~")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  dplyr::select(lhs, rhs, est, ci.lower, ci.upper)

kable(latent_cov_ci, digits = 3, caption = "CFA Standardized Factor Loadings with 95% Bootstrap CI")

#--------------------------------------------------
# Standardized latent correlations with CI
#--------------------------------------------------

cor_ci <- standardizedSolution(
  overall.fit,
  type = "std.all",
  ci = TRUE
) %>%
  filter(
    op == "~~",
    lhs %in% factors,
    rhs %in% factors,
    lhs != rhs
  ) %>%
  rowwise() %>%
  mutate(pair = paste0(sort(c(lhs, rhs)), collapse = "~~")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(lhs, rhs, est.std, ci.lower, ci.upper) %>%
  rename(
    Factor1 = lhs,
    Factor2 = rhs,
    Correlation = est.std,
    CI_lower = ci.lower,
    CI_upper = ci.upper
  )

kable(
  cor_ci,
  digits = 3,
  caption = "Standardized Latent Factor Correlations with 95% Bootstrap CI"
)

#-------------------------------
# MGCFA
#-------------------------------
variables <- c("GFF_Raven", "GFF_Figures", "GAF_Firstsounds", "GAF_Lastsounds", 
               "PSMF_Bodyparts", "PSMF_Fingers" , "PSMF_Hands", "WMF_FDS", "WMF_BDS", 
               "WMF_Corsi", "BF_BOR", "BF_BOL", "BF_BCR", "BF_BCL")
data_cleaned <- data %>%
  filter(complete.cases(select(., all_of(variables))))
nrow(data_cleaned)

variables <- c("GFF_Raven", "GFF_Figures", "GAF_Firstsounds", "GAF_Lastsounds", 
               "PSMF_Bodyparts", "PSMF_Fingers" , "PSMF_Hands", "WMF_FDS", "WMF_BDS", 
               "WMF_Corsi", "BF_BOR", "BF_BOL", "BF_BCR", "BF_BCL")

data_cleaned <- data %>%
  filter(complete.cases(select(., all_of(variables))))

MG.model <- eqMI.main(model = overall.model,
                      data = data_cleaned,
                      group = "gender",
                      meanstructure = TRUE,
                      output = "both",
                      equivalence.test = TRUE,
                      adjRMSEA = TRUE,
                      projection = TRUE,
                      bootstrap = FALSE)
ls()

## MEN

summary(MG.model$convention.sem$LavaanOut$fit.configural.g1,
        standardized = TRUE,
        rsquare = TRUE,
        fit.measure = TRUE)

## WOMEN

summary(MG.model$convention.sem$LavaanOut$fit.configural.g2,
        standardized = TRUE,
        rsquare = TRUE,
        fit.measure = TRUE)

## overall tables

table_fit[2, ] <- c("Men model",
                    round(fitmeasures(MG.model$convention.sem$LavaanOut$fit.configural.g1,
                                      c("chisq", "df", "cfi","tli", "rmsea", "srmr", "aic", "bic")),3))

table_fit[3, ] <- c("Women model",
                    round(fitmeasures(MG.model$convention.sem$LavaanOut$fit.configural.g2,
                                      c("chisq", "df", "cfi", "tli", "rmsea", "srmr", "aic", "bic")),3))

kable(table_fit)

## Configural invariance
summary(MG.model$convention.sem$LavaanOut$fit.combine.groups,
        standardized = TRUE,
        rsquare = TRUE,
        fit.measure = TRUE)

table_fit[4, ] <- c("configural model", round(fitmeasures(MG.model$convention.sem$LavaanOut$fit.combine.groups,
                                                          c("chisq","df", "cfi", "tli", "rmsea", "srmr", "aic", "bic")),3))

kable(table_fit)

## Metric invariance
summary(MG.model$convention.sem$LavaanOut$fit.metric,
        standardized = TRUE,
        rsquare = TRUE,
        fit.measure =TRUE)

table_fit[5, ] <- c("Metric model", round(fitmeasures(MG.model$convention.sem$LavaanOut$fit.combine.groups,
                                                      c("chisq","df", "cfi", "tli", "rmsea", "srmr", "aic", "bic")),3))

kable(table_fit)

## Scalar invariance
summary(MG.model$convention.sem$LavaanOut$fit.scalar,
        standardized = TRUE,
        rsquare = TRUE,
        fit.measure = TRUE)

table_fit[6, ] <- c("Scalar model", round(fitmeasures(MG.model$convention.sem$LavaanOut$fit.scalar,
                                                      c("chisq","df", "cfi", "tli", "rmsea", "srmr", "aic", "bic")),3))

kable(table_fit)

## Strict (error) invariance (are the item residuals the same?)
summary(MG.model$convention.sem$LavaanOut$fit.strict.residuals,
        standardized = TRUE,
        rsquare = TRUE,
        fit.measure = TRUE)

table_fit[7, ] <- c("Strict model", round(fitmeasures(MG.model$convention.sem$LavaanOut$fit.strict.residuals,
                                                      c("chisq", "df", "cfi", "tli", "rmsea", "srmr", "aic", "bic")),3))

kable(table_fit)

##EQ testing
EQ.model <- eqMI.main(model = overall.model,
                      data = data_cleaned,
                      group = "gender",
                      meanstructure = TRUE,
                      output = "both", ## mean and covariance
                      equivalence.test = TRUE,
                      adjRMSEA = TRUE,
                      projection= TRUE,
                      bootstrap = FALSE,
                      quiet = T)

EQ.model$eqMI.stat

summary(data_cleaned)
descr(data_cleaned, stats = c("n", "mean", "sd", "min", "max", "skewness", "kurtosis"))



