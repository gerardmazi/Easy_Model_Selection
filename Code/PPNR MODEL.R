################################################################################################################################
#
#
#                                                             PPNR MODELING
#                                                             
#
#
################################################################################################################################

# Optimization of model selection for time series regression
'
Gerard Mazi
08/21/2018
gerard.mazi@gmail.com
862.221.2477
'
library(data.table)                 # For processing data
library(dplyr)                      # For processing data
library(plyr)                       # For processing data
library(zoo)                        # Working with time
library(leaps)                      # Perform model selection (best subsets, stepwise, etc.)
library(car)                        # Voad vif() function for multicollinearity
library(lmtest)                     # Breusch-Pagani and White test for heteroskedasticity
library(sandwich)                   # Performs robust coefficients
library(gvlma)                      # Global validation of liner model assumptions
library(caret)                      # Box-Cox transformation for heteroskedasticity
library(e1071)                      # Required to run Box-Cox transformation
library(fUnitRoots)                 # adfTest function to test for stationarity
library(nortest)                    # Normality test


wd = 'C:/Users/gmazi/Desktop/PPNR/Model'
input = paste0(wd, '/Input')
output = paste0(wd, '/Output')
code = paste0(wd, '/Code')
setwd(wd)

#===============================================================================================================================
# BANKING DATA IMPORT
#===============================================================================================================================
# Load aggregated banking data
bank = fread(input = paste0(input, '/data_industry.csv'), stringsAsFactors = F, 
             drop = c('loan','SF_orig_total','total_assets'))

# Add index
bank = cbind(index = seq_along(1:nrow(bank)),bank)

# Convert quarter to date format
bank$quarter = as.yearqtr(bank$quarter, format = "%Y Q%q")



#===============================================================================================================================
# DFAST MACROECONOMIC DATA
#===============================================================================================================================
# Load DFAST variables which include historical macroeconomic data and stressed scenarios
macro = fread(paste0(input, '/DFAST_base.csv'), stringsAsFactors = F)
macro_adv = fread(paste0(input, '/DFAST_adverse.csv'), stringsAsFactors = F)
macro_sev = fread(paste0(input, '/DFAST_severe.csv'), stringsAsFactors = F)

#lapply(list(macro,macro_adv,macro_sev), function(x) mutate(x, BBB_SP = BBB - T10Y))
#purrr::map(list(macro,macro_adv,macro_sev), function(x) mutate(x, BBB_SP = BBB - T10Y))

# Calc BBB spread
macro$BBB_SP = macro$BBB - macro$T10Y
macro_adv$BBB_SP = macro_adv$BBB - macro_adv$T10Y
macro_sev$BBB_SP = macro_sev$BBB - macro_sev$T10Y

# Calc Treasury Yield Spread (10 Year Treasury less 3 Month Treasury)
macro$T_SP = macro$T10Y - macro$T3M
macro_adv$T_SP = macro_adv$T10Y - macro_adv$T3M
macro_sev$T_SP = macro_sev$T10Y - macro_sev$T3M

# Convert quarter to date format for consitency and merging with bank data
macro$quarter = as.yearqtr(macro$quarter, format = "Q%q %Y")
macro_adv$quarter = as.yearqtr(macro_adv$quarter, format = "Q%q %Y")
macro_sev$quarter = as.yearqtr(macro_sev$quarter, format = "Q%q %Y")

# Define growth vs. differencing for each variable
growth = c("GDPR","GDPN","DPIR","DPIN","CPI","DJ","HPI","CRE_PI","VIX")
diff = c("UR","T3M","T5Y","T10Y","BBB","MTG","PRIM","BBB_SP", "T_SP")
all_vars = c(growth, diff)

# Transform all variables where appropriate: QoQ diff or growth; YoY diff or growth; Lag 1 through 4
#### Base
macro = mutate_at(macro, setNames(diff, paste0(diff, "_QoQ_D")), function(x) x - dplyr::lag(x))
macro = mutate_at(macro, setNames(diff, paste0(diff, "_YoY_D")), function(x) x - dplyr::lag(x, 4))
macro = mutate_at(macro, setNames(growth, paste0(growth, "_QoQ_G")), function(x) x / dplyr::lag(x) - 1)
macro = mutate_at(macro, setNames(growth, paste0(growth, "_YoY_G")), function(x) x / dplyr::lag(x, 4) - 1)

pre_lag = names(macro)[-1]

for (i in 1:4){
  macro = mutate_at(macro, setNames(pre_lag, paste0(pre_lag, "_L", i)), function(x) dplyr::lag(x, i))
}

### Adverse
macro_adv = mutate_at(macro_adv, setNames(diff, paste0(diff, "_QoQ_D")), function(x) x - dplyr::lag(x))
macro_adv = mutate_at(macro_adv, setNames(diff, paste0(diff, "_YoY_D")), function(x) x - dplyr::lag(x, 4))
macro_adv = mutate_at(macro_adv, setNames(growth, paste0(growth, "_QoQ_G")), function(x) x / dplyr::lag(x) - 1)
macro_adv = mutate_at(macro_adv, setNames(growth, paste0(growth, "_YoY_G")), function(x) x / dplyr::lag(x, 4) - 1)

for (i in 1:4){
  macro_adv = mutate_at(macro_adv, setNames(pre_lag, paste0(pre_lag, "_L", i)), function(x) dplyr::lag(x, i))
}

### Severe
macro_sev = mutate_at(macro_sev, setNames(diff, paste0(diff, "_QoQ_D")), function(x) x - dplyr::lag(x))
macro_sev = mutate_at(macro_sev, setNames(diff, paste0(diff, "_YoY_D")), function(x) x - dplyr::lag(x, 4))
macro_sev = mutate_at(macro_sev, setNames(growth, paste0(growth, "_QoQ_G")), function(x) x / dplyr::lag(x) - 1)
macro_sev = mutate_at(macro_sev, setNames(growth, paste0(growth, "_YoY_G")), function(x) x / dplyr::lag(x, 4) - 1)

for (i in 1:4){
  macro_sev = mutate_at(macro_sev, setNames(pre_lag, paste0(pre_lag, "_L", i)), function(x) dplyr::lag(x, i))
}

# Merge bank and macro data
bank = merge.data.frame(bank, macro, by = 'quarter', all.x = T, all.y = F)

# Exclude base (non-transformend) variables along with thir lags due to lack of stationarity
all_vars_lag = paste0(all_vars[!all_vars %in% c("BBB_SP", "T_SP")], "_L")           # Keep BBB and Treasury spreads

bank = bank[!names(bank) %in% all_vars[!all_vars %in% c("BBB_SP", "T_SP")]]
bank = bank[, !grepl(paste(all_vars_lag, collapse = "|"), colnames(bank))]

# output forecast quarters from each DFAST dataset (i.e., > 2018 Q2)
fcst_base = macro[macro$quarter > "2018 Q2",]
fcst_adv = macro_adv[macro_adv$quarter > "2018 Q2",]
fcst_sev = macro_sev[macro_sev$quarter > "2018 Q2",]

calc_4qtr_ma = function(var_data){
  new = c()
  new[1] = var_data[1]
  new[2] = mean(var_data[1:2], na.rm = T)
  new[3] = mean(var_data[1:3], na.rm = T)
  for(i in 4:length(var_data)){
    new[i] = mean(var_data[(i-3):i], na.rm=T)
  }
  return(new)
}

#===============================================================================================================================
# MODEL SELECTION
#===============================================================================================================================
# Generate dataset for modeling
loan_data = select(bank, quarter, index, SF_loans,
                   starts_with("GDP"),
                   starts_with("DPI"),
                   starts_with("UR"),
                   starts_with("CPI"),
                   starts_with("T"),
                   starts_with("BBB"),
                   starts_with("MTG"),
                   starts_with("PRIM"),
                   starts_with("HPI"),
                   starts_with("CRE_PI"))


# Loess transform
#loan_data = mutate(loan_data, SF_loans = predict(loess(SF_loans ~ index, data = loan_data, span = 0.6)))

# 4 Quarter moving average transform
#loan_data$SF_loans = calc_4qtr_ma(loan_data$SF_loans)

# Percent of total assets transform
#loan_data = mutate(loan_data, SF_loans = (SF_loans) / (total_assets))

# Quarterly growth transform
loan_data = mutate(loan_data, SF_loans = SF_loans / dplyr::lag(SF_loans) - 1)
loan_data = na.omit(loan_data)    # Remove NAs

# Perform a best subsets model selection
best = regsubsets(SF_loans ~ ., data = loan_data[,-(1:2)], nbest = 10, nvmax = 2, method = 'exhaustive', really.big = T)

# Output all factors in the best subsets output in a vectorized matrix
regMat = summary(best)$which[,-1]

# List out all regressors
regressors = colnames(regMat)

# Build regression formulas with list or regressors and the vectorized matrix
allModelsList = apply(regMat, 1, function(x) as.formula(paste(c("SF_loans ~ 1", regressors[x]), collapse=" + ")))

# Run all models individually and save the results to this object for future diagnostics analysis
allModelsResults = lapply(allModelsList, function(x) lm(x, data=loan_data))

# Output various diagnostics
dfpValues = as.data.frame(do.call(rbind, lapply(allModelsResults, function(x) (coef(summary(x))[,"Pr(>|t|)"]))))
names(dfpValues) = c("intcpt_p-val", paste0("p-val_", seq_along(1:(ncol(dfpValues)-1))))
NoOfCoef = unlist(apply(regMat, 1, sum))
adjR2    = unlist(lapply(allModelsResults, function(x) summary(x)$adj.r.squared))
BIC = unlist(lapply(allModelsResults, function(x) BIC(x)))

# Output the coefficients
coef_out = as.data.frame(do.call(rbind, lapply(allModelsResults, function(x) summary(x)$coef[,"Estimate"])))
colnames(coef_out)[-1] = paste0("coef_", seq_along(1:(ncol(coef_out)-1)))

# Output the Durbin-Watson statistic for autocorrelation
DW = ldply(allModelsResults, function(x) as.data.frame(t(dwtest(x)$statistic)))[-1]

# Output the VIF statistic for multicollinearity
vif = rbind(
  data.frame(vif_1 = rep(NA, 10), 
             vif_2 = rep(NA, 10), 
             vif_3 = rep(NA, 10)),
  data.frame(vif_1 = unlist(lapply(allModelsResults[-(1:10)], function(x) vif(x)[1])),
             vif_2 = unlist(lapply(allModelsResults[-(1:10)], function(x) vif(x)[2])),
             vif_3 = unlist(lapply(allModelsResults[-(1:10)], function(x) vif(x)[3])))
  )

# Output the Breusch-Pagan p-value for heteroskedasticity
BP = data.frame(
  BP = unlist(lapply(allModelsResults, function(x) bptest(x)[[1]])),
  BP_p = unlist(lapply(allModelsResults, function(x) bptest(x)[[4]]))
  )

# Output normality tests
normality = cbind(
  shapiro = unlist(lapply(allModelsResults, function(x) shapiro.test(resid(x))[[2]])),
  ks_test = unlist(lapply(allModelsResults, function(x) ks.test(resid(x), 'pnorm')[[2]])),
  cvm_test = unlist(lapply(allModelsResults, function(x) cvm.test(resid(x))[[2]])),
  ad_test = unlist(lapply(allModelsResults, function(x) ad.test(resid(x))[[2]]))
  )

# Output predictions for all three DFAST scenarios
pred = cbind(
  ldply(allModelsResults, function(x) as.data.frame(t(predict(x, fcst_base)))),
  ldply(allModelsResults, function(x) as.data.frame(t(predict(x, fcst_adv)))),
  ldply(allModelsResults, function(x) as.data.frame(t(predict(x, fcst_sev))))
  )

# Aggregate all diagnostics, sort, and output to csv
results = data.frame(model = as.character(allModelsList),
                     NoOfCoef = NoOfCoef,
                     coef_out,
                     dfpValues,
                     adjR2 = adjR2,
                     BIC = BIC,
                     DW,
                     vif,
                     BP,
                     normality,
                     pred)
results = results[order(results$NoOfCoef, -results$adjR2),]
fwrite(results, "Output/results.csv")

