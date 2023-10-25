################################################################################
##Compare baseline characteristics of Anchor vs Main
################################################################################

# Display baseline characteristics of ANCHOR trial (binary variables reformatted as percentages)
target_pop_standard %>%
  mutate(
    SEX=sprintf("%1.0f%%", 100*SEX),
    ECOG=sprintf("%1.0f%%", 100*ECOG),
    ISS=sprintf("%1.0f%%", 100*ISS)
  )

# Summarise baseline characteristics of MAIN IPD trial at aggregate level
intervention_input %>%
  select(
    AGE,
    SEX,
    ECOG,
    ISS
  ) %>%
  drop_na() %>%
  summarise(
    AGE=round(mean(AGE),2),
    SEX=paste0(
      as.character(round(sum(intervention_input$SEX)/length(intervention_input$SEX)*100),4),'%'),
    ECOG=paste0(
      as.character(round(sum(intervention_input$ECOG)/length(intervention_input$ECOG)*100),4),'%'),
    ISS=paste0(
      as.character(round(sum(intervention_input$ISS)/length(intervention_input$ISS)*100),4),'%')
  )

################################################################################
# Check effect modifiers in MAIN trial
################################################################################
cox_all_covariates_main <-
  coxph(
    Surv(Time, Event)~
      ARM + AGE + SEX + ECOG + ISS + 
      ARM*AGE + ARM*SEX + ARM*ECOG + ARM*ISS, 
    data=intervention_input  
  )


##Try differemt significance levels
summary(cox_all_covariates_main, conf.int=0.9)



################################################################################
###Matching Procedure
################################################################################

# Select vars to include in MAIC
maic_vars <- c(
  "AGE",
  "SEX",
  "ECOG",
  "ISS"
)



##center baseline characteristics
main_centered <- intervention_input %>%
  mutate(AGE = AGE - target_pop$age,
         SEX = SEX - target_pop$prop.male,
         ECOG = ECOG - target_pop$prop.ecog1,
         ISS = ISS - target_pop$prop.iss3
  ) %>%
  drop_na() # Standard MAIC can't handle missing values


# Complicated looking function to calculate weights for MAIC by Newton-Raphson
objfn <- function(a1, X){
  sum(exp(X %*% a1))
} # This is a function for Q(b) defined above.


# Gradient function => Derivative of Q(b).
gradfn <- function(a1, X){
  colSums(sweep(X, 1, exp(X %*% a1), '*'))
}


## The following is the in built R function used to optimise Q(b) ##
## using Newton-Raphson techniques ##
print(
  opt1 <- optim(
    par = rep(0,dim(main_centered %>% select(maic_vars))[2]) ,
    fn = objfn, gr = gradfn, 
    X = as.matrix(main_centered %>% select(maic_vars)),
    method = 'BFGS'
  )
)

a1 = opt1$par


# Create dataset for analysis by making sure weights are attached to dataset
df_analysis <- intervention_input %>%
  drop_na() %>%
  mutate(wt = as.vector(exp(as.matrix(main_centered %>% select(maic_vars)) %*% a1)))

wt = as.vector(exp(as.matrix(main_centered %>% select(maic_vars)) %*% a1))



################################################################################
#### Weight Diagnostics ####
################################################################################
# Check estimation of weights has worked- miniumum requirement is that the
# baseline characteristics match at the aggregate level
reweighted_anchor_baseline <- as.data.frame(
  df_analysis %>%
    drop_na() %>%
    bind_cols(as.data.frame(wt)) %>%
    summarise(
      AGE = weighted.mean(AGE, wt),
      SEX = 100*weighted.mean(SEX,wt),
      ECOG = 100*weighted.mean(ECOG,wt),
      ISS = 100*weighted.mean(ISS,wt)
    )
)

reweighted_anchor_baseline
target_pop_standard

# Effective sample size
ESS = sum(wt)^2/sum(wt^2)
ESS

# Re-scaled weights
wt.rs <- (wt/sum(wt))*nrow(df_analysis)
qplot(wt.rs, geom='histogram', xlab='rescaled weight', binwidth=1.5)


################################################################################
###INDIRECT COMPARISON###
################################################################################
# Adjusted and unadjusted hazard ratios from MAIN
cox_unadjusted <- coxph(
  Surv(Time, Event)~ ARM ,
  data=intervention_input
)

cox_adjusted <- coxph(
  Surv(Time, Event)~ ARM ,
  data=df_analysis,
  weights=wt
)


################################################################################
##SPECIFY DATA FOR ANCHOR TRIAL
# Anchor log hazard ratio and standard error
anchor_log_hr_pe <- log(0.712)
anchor_log_hr_se <- (log(0.906) - log(0.56))/(2*1.96)


# Main log hazard ratio and standard error- adjusted and unadjusted
main_unadjusted_log_hr_pe <- log(summary(cox_unadjusted)$conf.int[1])
main_unadjusted_log_hr_se <- (log(summary(cox_unadjusted)$conf.int[4])- log(summary(cox_unadjusted)$conf.int[3]))/(2*1.96)
main_adjusted_log_hr_pe <- log(summary(cox_adjusted)$conf.int[1])
main_adjusted_log_hr_se <- (log(summary(cox_adjusted)$conf.int[4])- log(summary(cox_adjusted)$conf.int[3]))/(2*1.96)

#ITC- adjusted and unadjusted
adjusted_itc_pe <- exp(main_adjusted_log_hr_pe - anchor_log_hr_pe)
adjusted_itc_se <- sqrt((anchor_log_hr_se^2) + (main_adjusted_log_hr_se^2))
adjusted_itc_lci <- exp(log(adjusted_itc_pe) - adjusted_itc_se)
adjusted_itc_uci <- exp(log(adjusted_itc_pe) + adjusted_itc_se)

unadjusted_itc_pe <- exp(main_unadjusted_log_hr_pe - anchor_log_hr_pe)
unadjusted_itc_se <- sqrt((anchor_log_hr_se^2) + (main_unadjusted_log_hr_se^2))
unadjusted_itc_lci <- exp(log(unadjusted_itc_pe) - unadjusted_itc_se)
unadjusted_itc_uci <- exp(log(unadjusted_itc_pe) + unadjusted_itc_se)


##display results###

HR_CI <- as.data.frame(cbind("HR" = unadjusted_itc_pe, "HR_low_CI" = unadjusted_itc_lci, "HR_upp_CI" = unadjusted_itc_uci))
HR_CI

HR_CI_adj <- as.data.frame(cbind("HR" = adjusted_itc_pe, "HR_low_CI" = adjusted_itc_lci, "HR_upp_CI" = adjusted_itc_uci))
HR_CI_adj




################################################################################
## SUMMARY ##
################################################################################

## Summary

# Produce a summary of HRs and CIs
HR_summ <-  rbind(HR_CI, HR_CI_adj) %>% # Adjusted and unadjusted
  mutate(Method = c("Unadjusted HR (95% CI)",
                    "Adjusted HR (95% CI)")) %>%
  
  #apply rounding for numeric columns
  mutate_if(.predicate = is.numeric, sprintf, fmt = "%.3f") %>%
  #format for output
  transmute(Method, HR_95_CI = paste0(HR, " (", HR_low_CI, " to ", HR_upp_CI, ")"))


# turns the results to a table suitable for word/ powerpoint
HR_table <- HR_summ %>%
  regulartable() %>% #make it a flextable object
  set_header_labels(Method = "Method",  HR_95_CI = "Hazard ratio (95% CI)")  %>%
  font(font = 'Arial', part = 'all') %>%
  fontsize(size = 14, part = 'all') %>%
  bold(part = 'header') %>%
  align(align = 'center', part = 'all') %>%
  align(j = 1, align = 'left', part = 'all') %>%
  border_outer(border = fp_border()) %>%
  border_inner_h(border = fp_border()) %>%
  border_inner_v(border = fp_border()) %>%
  autofit(add_w = 0.2, add_h = 2)

HR_table

####FOREST PLOT###
# Forest plot of results
plotdat <- data_frame(
  id = 1:5,
  Comparison = factor(c(1, 1, 2, 2, 3),
                      labels = c("DRd vs VRd", "DRd vs Rd", "VRd vs Rd")),
  LogHR = c(log(adjusted_itc_pe), log(unadjusted_itc_pe), 
            main_adjusted_log_hr_pe, main_unadjusted_log_hr_pe, 
            anchor_log_hr_pe),
  HR = c(adjusted_itc_pe, unadjusted_itc_pe, 
         exp(main_adjusted_log_hr_pe), exp(main_unadjusted_log_hr_pe), 
         exp(anchor_log_hr_pe)),
  SE = c(adjusted_itc_se, unadjusted_itc_se,
         main_adjusted_log_hr_se, main_unadjusted_log_hr_se,
         anchor_log_hr_se),
  lo = exp(LogHR - SE),
  hi = exp(LogHR + SE),
  type = c("Adjusted", "Unadjusted",
           "Adjusted", "Unadjusted",
           "Unadjusted")
)



ggplot(aes(x = HR, y = id, col = type, shape = type), data = plotdat) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_point(size = 3.5) +
  geom_segment(aes(y = id, yend = id, x = lo, xend = hi), na.rm = TRUE) +
  xlab("Hazards Ratio") +
  facet_grid(Comparison~., switch = "y", scales = "free_y", space = "free_y") +
  scale_y_reverse(name = "", breaks = NULL, expand = c(0, 0.6))

