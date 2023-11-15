################################################################################
## DEFINE PARAMETERS ##
################################################################################
v.n            <- c("Pre-Progression","Progressed","Dead") # state names
n.s            <- length(v.n)    # number of states
n.t            <- 30*52.143      # number of cycles to run / time horizon in weeks
c.l            <- 1              # cycle length
times          <- seq(from = 0, to = n.t, by = c.l)  # sequence of times to be considered in the model

################################################################################
## GENERATE KM AND PLOTS ##
################################################################################
##km
km.drdpfs <- survfit(Surv(Time, Event)~1, data=drdpfs, type = "kaplan-meier")
km.drdttd <- survfit(Surv(Time, Event)~1, data=drdttd, type = "kaplan-meier")
km.drdos <- survfit(Surv(Time, Event)~1, data=drdos, type = "kaplan-meier")
km.vrdpfs <- survfit(Surv(Time, Event)~1, data=vrdpfs, type = "kaplan-meier")
km.vrdos <- survfit(Surv(Time, Event)~1, data=vrdos, type = "kaplan-meier")

#plot PFS curves together
plot(km.vrdpfs, col = 2, main = "Progression-free Survival", xlab = "Time (weeks)", ylab = "Survival (%)", lwd = 1, conf.int= F)
lines(km.drdpfs, col = 1, conf.int= F)
legend("bottomleft", col = c(1,2), bty = "n", lwd = rep(1), c("DRd", "VRd"))

#plot OS curves together
plot(km.vrdos, col = 2, main = "Overall Survival", xlab = "Time (weeks)", ylab = "Survival (%)", lwd = 1, conf.int= F)
lines(km.drdos, col = 1, conf.int= F)
legend("bottomleft", col = c(1,2), bty = "n", lwd = rep(1), c("DRd", "VRd"))


################################################################################
## CREATE FUNCTION FOR FITTING PARAMETRIC DISTRIBUTIONS ##
################################################################################

fit.pm <- function(time = "Time", event = "Event", data = data, label = label)  
{
  data$time  = data[, time]  
  data$event = data[, event]
  
  # Progression free survival  
  KM.fit     <- survfit(Surv(time, event) ~ 1, data = data, type = "kaplan-meier")# fit Kaplan-Meier curve 
  fit.exp    <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "exp"    ) # fit model with exponential distribution 
  fit.weib   <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "weibull") # fit model with Weibull distribution
  fit.gg     <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "gengamma")# fit model with gamma distribution 
  fit.lnorm  <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "lnorm"  ) # fit model with lognormal distribution
  fit.llogis <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "llogis" ) # fit model with loglogistic distribution
  fit.gomp   <- flexsurvreg(Surv(time, event) ~ 1, data = data, dist = "gompertz" ) # fit model with gompertz distribution
  
  # extarapolate all models beyond the KM curve
  plot(KM.fit, ylab = "Survival Probability", xlab = "Time (weeks)", ylim = c(0,1), xlim = c(0,n.t), conf.int= F, main = label)
  lines(fit.exp ,     t = times, col = cols[1], ci = F)
  lines(fit.weib,     t = times, col = cols[2], ci = F)
  lines(fit.gg,       t = times, col = cols[3], ci = F)
  lines(fit.lnorm,    t = times, col = cols[4], ci = F)
  lines(fit.llogis,   t = times, col = cols[5], ci = F)
  lines(fit.gomp,     t = times, col = cols[6], ci = F)
  lines(KM.fit, conf.int = F)
  legend("topright", cex = 0.9, col = c("black", cols[1], cols[2], cols[3], cols[4], cols[5], cols[6]),
         bty = "n", lwd = rep(2, 4),
         c("Kaplan-Meier", "Exponential", "Weibull", "Gen. Gamma", "Log-normal", "Log-logistic", "Gompertz"))
  
  # compare AIC values
  AIC <- c(    Exponential = AIC(fit.exp),
               Weibull     = AIC(fit.weib), 
               GenGamma    = AIC(fit.gg),
               Lognormal   = AIC(fit.lnorm), 
               Loglogistic = AIC(fit.llogis),
               Gompertz    = AIC(fit.gomp))
  
  AIC= round(AIC,3)
  
  # compare BIC values
  BIC <- c(    Exponential = BIC(fit.exp),
               Weibull     = BIC(fit.weib), 
               GenGamma    = BIC(fit.gg),
               Lognormal   = BIC(fit.lnorm), 
               Loglogistic = BIC(fit.llogis),
               Gompertz    = BIC(fit.gomp))
  
  BIC= round(BIC,3)
  
  gof <- data.frame(AIC, BIC)
  
  res <- list(Exponential = fit.exp,
              Weibull     = fit.weib,
              GenGamma    = fit.gg,
              Lognormal   = fit.lnorm, 
              Gamma       = fit.gg,
              Loglogistic = fit.llogis,
              Gompertz    = fit.gomp,
              AIC         = AIC,
              BIC         = BIC,
              GOF         = gof)
  res
}

################################################################################
## FITTING PARAMETRIC DISTRIBUTIONS 
################################################################################

fit.drdpfs <- fit.pm(time = "Time", event = "Event", data = drdpfs, label = "DRd PFS") #Fit drd pfs
fit.drdttd <- fit.pm(time = "Time", event = "Event", data = drdttd, label = "DRd TTD") #Fit drd ttd
fit.vrdpfs <- fit.pm(time = "Time", event = "Event", data = vrdpfs, label = "VRd PFS") #Fit drd pfs
fit.drdos <- fit.pm(time = "Time", event = "Event", data = drdos, label = "DRd OS") #Fit drd os
fit.vrdos <- fit.pm(time = "Time", event = "Event", data = vrdos, label = "VRd OS") #Fit vrd os

##display AIC and BIC results
fit.drdpfs[[10]]
fit.drdttd[[10]]
fit.vrdpfs[[10]]
fit.drdos[[10]]
fit.vrdos[[10]]

#select best fitting models
b.drdpfs <- fit.drdpfs$Exponential
plot(km.drdpfs, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted PFS for DRd, Exponential")) 
lines(b.drdpfs,  col = 2, t = times, lty = 2)

b.drdttd <- fit.drdttd$Exponential
plot(km.drdttd, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted TTD for DRd, Exponential")) 
lines(b.drdttd,  col = 2, t = times, lty = 2)

b.vrdpfs <- fit.vrdpfs$Exponential
plot(km.vrdpfs, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted PFS for VRd, Exponential")) 
lines(b.vrdpfs,  col = 2, t = times, lty = 2)

b.drdos <- fit.drdos$Exponential
plot(km.drdos, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted OS for DRd, Exponential")) 
lines(b.drdos,  col = 2, t = times, lty = 2)   

b.vrdos <- fit.vrdos$Exponential
plot(km.vrdos, ylab = "Survival Probability", xlab = "Time (months)",  main = paste ("True vs Fitted OS for VRd, Exponential")) 
lines(b.vrdos,  col = 2, t = times, lty = 2)   


#select most conservative models
c.drdpfs <- fit.drdpfs$Exponential
plot(km.drdpfs, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted PFS for DRd, Exponential")) 
lines(c.drdpfs,  col = 4, t = times, lty = 2)

c.drdttd <- fit.drdttd$GenGamma
plot(km.drdttd, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted TTD for DRd, Gen. Gamma")) 
lines(c.drdttd,  col = 4, t = times, lty = 2)

c.drdos <- fit.drdos$GenGamma
plot(km.drdos, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted OS for DRd, Gen. Gamma")) 
lines(c.drdos,  col = 4, t = times, lty = 2)   

c.vrdos <- fit.vrdos$Gompertz
plot(km.vrdos, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted OS for VRd, Gompertz")) 
lines(c.vrdos,  col = 4, t = times, lty = 2)  

#Weibull models for all
w.drdpfs <- fit.drdpfs$Weibull
plot(km.drdpfs, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted PFS for DRd, Weibull")) 
lines(w.drdpfs,  col = 4, t = times, lty = 2)

w.drdttd <- fit.drdttd$Weibull
plot(km.drdttd, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted TTD for DRd, Weibull")) 
lines(w.drdttd,  col = 4, t = times, lty = 2)

w.drdos <- fit.drdos$Weibull
plot(km.drdos, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted OS for DRd, Weibull")) 
lines(w.drdos,  col = 4, t = times, lty = 2)   

w.vrdos <- fit.vrdos$Weibull
plot(km.vrdos, ylab = "Survival Probability", xlab = "Time (months)",  main = paste ("True vs Fitted OS for VRd, Weibull")) 
lines(w.vrdos,  col = 4, t = times, lty = 2)   


################################################################################
## CREATE FUNCTION FOR PARTITIONED SURVIVAL ANALYSIS ##
################################################################################

partsa <- function(fit.pfs, fit.os, title = "treatment", time = times)
{
  
  pfs.surv             <- summary(fit.pfs, t = time, ci = F)[[1]]$est
  os.surv              <- summary(fit.os, t = time, ci = F)[[1]]$est
  cycle                <- time
  year                 <- cycle/52.143
  pfs                  <- pfs.surv
  os                   <- os.surv
  preprog              <- pfs.surv                    # probability of remaining progression-free
  prog                 <- os.surv - pfs.surv          # estimate the probability of remaining in the progressed state
  prog[prog < 0]       <- 0                           # in cases where the probability is negative replace with zero
  dead                 <- 1 - prog - preprog          # probability of being dead
  check                <- prog + preprog + dead
  trace1               <- data.frame(year = year, cycle = cycle, pfs = pfs, os = os, preprog = preprog, prog = prog, dead = dead, check = check)
  trace2               <- data.frame(preprog = preprog, prog = prog, dead = dead)
  res                  <- list(complete = trace1, matrix = trace2)
  return(res)
}

################################################################################
## CREATE PARTITIONED SURVIVAL MATRIX
################################################################################

## apply PSM function to selected curve ## CHANGE TO BEST FIT OR MOST CONSERVATIVE
psm.drd <- partsa(b.drdttd, b.drdos, time = times)
psm.vrd <- partsa(b.vrdpfs, b.vrdos, time = times)

#plot matrices
matrix.drd <- as.matrix(psm.drd$matrix)    
matplot(matrix.drd, type = 'l', lty=1, xlab = "Time (weeks)", ylab = "Proportion of patients", main = "DRd")                      
legend("right", v.n, col=1:n.s, lty=rep(1,n.s), bty='n')

matrix.vrd <- as.matrix(psm.vrd$matrix)    
matplot(matrix.vrd, type = 'l', lty=1, xlab = "Time (weeks)", ylab = "Proportion of patients", main = "VRd")                      
legend("right", v.n, col=1:n.s, lty=rep(1,n.s), bty='n')

trace.drd <- psm.drd$complete
trace.vrd <- psm.vrd$complete 


################################################################################
## FIT PIECEWISE EXPONENTIAL MODEL FOR VRD OS
################################################################################

c <- c(13,162,385)

os.pw <- survSplit(vrdos,cut=c,end="Time",event="Event",start="start") 
os.pw$gt1 <- (os.pw$start==0)+0
os.pw$gt2 <- (os.pw$start==c[1])+0
os.pw$gt3 <- (os.pw$start==c[2])+0
os.pw$gt4 <- (os.pw$start==c[3])+0

# Fit piecewise exponential model 

pw.exp <- flexsurvreg(Surv(start, Time, Event)~gt1 + gt2 + gt3,
                      dist="exponential", data=os.pw)
h1 <- exp(pw.exp$res.t[1]+pw.exp$res.t[2])
h2 <- exp(pw.exp$res.t[1]+pw.exp$res.t[3])
h3 <- exp(pw.exp$res.t[1]+pw.exp$res.t[4])
h4 <- exp(pw.exp$res.t[1])

pw.exp$AIC


plot(pw.exp)

################################################################################
## GENERATE CURVES USING HAZARDS RATIO
################################################################################

##define HRs from ITC
##DRd versus VRd HR
hr.pfs <- 0.7211
ll.pfs <- 0.5048
ul.pfs <- 1.0300

hr.itc <- 1.1309
ll.os <- 0.8942
ul.os <- 1.4303
hr.os <- 1


##PFS CALCULATIONS
scale.pfs <- b.drdpfs$coefficients #get coefficient / scale parameter / lambda
invhr.pfs = 1+(1-hr.pfs) #get HR for VRd to DRd
loghr.pfs = log(invhr.pfs)
haz.pfs = exp(scale.pfs + loghr.pfs) #get new haz for comparator arm
riskprob.pfs = 1-exp(-haz.pfs) #get per cycle prob of event


#Set up for new trace
trace.vrd$pfs <- 1


#calculate cumulative survival
for(i in 2:nrow(trace.vrd)) {
  trace.vrd$pfs[i] <- trace.vrd$pfs[i-1]*(1-riskprob.pfs)
}

p1 = 1-exp(-h1)
p2 = 1-exp(-h2)
p3 = 1-exp(-h3)
p4 = 1-exp(-h4)

for(i in 2:nrow(trace.vrd)) {
  trace.vrd$os[i] <- trace.vrd$os[i-1]*(1-p1)
}

for(i in 14:nrow(trace.vrd)) {
  trace.vrd$os[i] <- trace.vrd$os[i-1]*(1-p2)
}

for(i in 163:nrow(trace.vrd)) {
  trace.vrd$os[i] <- trace.vrd$os[i-1]*(1-p3)
}



#recompute trace for DRD
trace.drd$preprog                     <- trace.drd$pfs  # probability of remaining progression-free
trace.drd$prog                        <- trace.drd$os - trace.drd$pfs   # estimate the probability of remaining in the progressed state
trace.drd$prog[trace.drd$prog < 0]    <- 0  # in cases where the probability is negative replace with zero
trace.drd$dead                        <- 1 - trace.drd$prog - trace.drd$preprog   # probability of being dead
trace.drd$check                       <- trace.drd$prog + trace.drd$preprog + trace.drd$dead

#recompute trace for VRD
trace.vrd$preprog                     <- trace.vrd$pfs  # probability of remaining progression-free
trace.vrd$prog                        <- trace.vrd$os - trace.vrd$pfs   # estimate the probability of remaining in the progressed state
trace.vrd$prog[trace.vrd$prog < 0]    <- 0  # in cases where the probability is negative replace with zero
trace.vrd$dead                        <- 1 - trace.vrd$prog - trace.vrd$preprog   # probability of being dead
trace.vrd$check                       <- trace.vrd$prog + trace.vrd$preprog + trace.vrd$dead


################################################################################
## PARAMETERS FOR COST-EFFECTIVENESS ANALYSIS ##
################################################################################
##### COST INPUTS
## DRUG COSTS
#Unit costs from hospital (minimum value)
dc_dara           <- #cost#
dc_bort           <- #cost#
dc_len            <- #cost#
dc_dex            <- #cost#
dc_drdcycle       <- ((dc_len*21)+(dc_dex*4*4))/4 #Rd cost of DRd per week
dc_vrdcycle       <- ((dc_len*21)+(dc_dex*2*8))/4 #Rd cost of VRd induction per week

#Chemo administration and disease monitoring
dc_admin          <- 184.02 #per visit
dc_mgt            <- 152.64/4 #per visit, divided by 4 to get weekly cost

##DRd (With PAP incorporated)
dc_drdc1_2        <- ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        <- ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          <- ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

##VRd
dc_vrdind         <- ((((dc_bort*1)+dc_admin)*4)/4)+dc_vrdcycle+dc_mgt
dc_vrdmaint       <- ((dc_len*21)+(dc_dex*2*4))/4 #Rd maintenance


#Cost of subsequent treatments
dc_pdrd           <- 1661.69/4
dc_pvrd           <- 4288.35/4

## HOSPICE
c_hospice         <- 8750/4

## HOSPITALIZATION FOR AE
cae_neut          <- 13266.3
cae_anae          <- 2567.49
cae_pneu          <- 13760.94

##### UTILITIES
## STATES per cycle (divided by weeks)
u_pf              <- 0.725/52.143
u_pp              <- 0.65/52.143

## AE DISUTILITY (NEGATIVE VALUES) per event
u_neut            <- (0.145*9.79)/365.25
u_anae            <- (0.25*1.94)/365.25
u_pneu            <- (0.14*11.14)/365.25
u_diar            <- (0.103*4.87)/365.25
u_neuro           <- (0.065*2.5)/365.25

#### AE RISK
rae_drdneut       <- 1-exp(-(-log(1-((136+61)/364)))/(56.2*4.345))
rae_drdleuk       <- 1-exp(-(-log(1-((41+19+37+5)/364)))/(56.2*4.345))
rae_drdpneu       <- 1-exp(-(-log(1-((62+5+3)/364)))/(56.2*4.345))
rae_drdanae       <- 1-exp(-(-log(1-((1+60)/364)))/(56.2*4.345))
rae_drdthromb     <- 1-exp(-(-log(1-((23+9)/364)))/(56.2*4.345))

rae_vrdneut       <- 1-exp(-(-log(1-((29+5+1)/241)))/(84*4.345))
rae_vrdlymph      <- 1-exp(-(-log(1-(5/241)))/(84*4.345))
rae_vrdanae       <- 1-exp(-(-log(1-((73+41)/241)))/(84*4.345))
rae_vrdthromb     <- 1-exp(-(-log(1-(7/241)))/(84*4.345))
rae_vrddiar       <- 1-exp(-(-log(1-((49+3+1)/241)))/(84*4.345))
rae_vrdneuro      <- 1-exp(-(-log(1-((76+4)/241)))/(84*4.345))


##### DISCOUNT RATE (ANNUAL)
dr_cost           <- 0.03
dr_qaly           <- 0.03


################################################################################
## ESTIMATE COSTS AND LIFE YEARS ##
################################################################################
## DRD; calculate costs, life years, and qalys per cycle
trace.drd$drugcosts    <- trace.drd$preprog*(ifelse(trace.drd$cycle<=24, ifelse(trace.drd$cycle<=8, dc_drdc1_2, dc_drdc3_6), dc_drdc7))
trace.drd$costs        <- trace.drd$drugcosts+(trace.drd$prog*(dc_pdrd))+(trace.drd$dead*(c_hospice))+(trace.drd$preprog*cae_neut*rae_drdneut)+(trace.drd$preprog*cae_anae*rae_drdanae)+(trace.drd$preprog*cae_pneu*rae_drdpneu)
trace.drd$disccosts    <- trace.drd$costs/(1 + dr_cost) ^ trace.drd$year
trace.drd$lifeyears    <- (trace.drd$preprog + trace.drd$prog)/(1 + dr_qaly) ^ trace.drd$year
trace.drd$qalys        <- (trace.drd$preprog*u_pf)+(trace.drd$prog*u_pp)-(trace.drd$preprog*rae_drdneut*u_neut)-(trace.drd$preprog*rae_drdpneu*u_pneu)-(trace.drd$preprog*rae_drdanae*u_anae)
trace.drd$discqalys    <- trace.drd$qalys/(1 + dr_qaly) ^ trace.drd$year

## VRD; calculate costs, life years, and qalys per cycle
trace.vrd$drugcosts    <- trace.vrd$preprog*(ifelse(trace.vrd$cycle<=24, dc_vrdind, dc_vrdmaint))
trace.vrd$costs        <- trace.vrd$drugcosts+(trace.vrd$prog*(dc_pvrd))+(trace.vrd$dead*(c_hospice))+(trace.vrd$preprog*cae_neut*rae_vrdneut)+(trace.vrd$preprog*cae_anae*rae_vrdanae)
trace.vrd$disccosts    <- trace.vrd$costs/(1 + dr_cost) ^ trace.vrd$year
trace.vrd$lifeyears    <- (trace.vrd$preprog + trace.vrd$prog)/(1 + dr_qaly) ^ trace.vrd$year
trace.vrd$qalys        <- (trace.vrd$preprog*u_pf)+(trace.vrd$prog*u_pp)-(trace.vrd$preprog*rae_vrdneut*u_neut)-(trace.vrd$preprog*rae_vrdanae*u_anae)-(trace.vrd$preprog*rae_vrddiar*u_diar)-(trace.vrd$preprog*rae_vrdneuro*u_neuro)
trace.vrd$discqalys    <- trace.vrd$qalys/(1 + dr_qaly) ^ trace.vrd$year


## make top row 0 values
trace.drd[1, c(9:14)] <- 0
trace.vrd[1, c(9:14)] <- 0

## SUMMARY OUTPUTS
## DRD
tpreprog.drd           <- sum(trace.drd$preprog) #get sum of those pre-progression
tc.drd                 <- sum(trace.drd$disccosts) #total cost
tly.drd                <- sum(trace.drd$lifeyears)/52.143 #get total life years
tqaly.drd              <- sum(trace.drd$discqalys) #total qalys

## VRD
tpreprog.vrd           <- sum(trace.vrd$preprog) #get sum of those pre-progression
tc.vrd                 <- sum(trace.vrd$disccosts) #total cost 
tly.vrd                <- sum(trace.vrd$lifeyears)/52.143 #get total life years
tqaly.vrd              <- sum(trace.vrd$discqalys) #total qalys

## create summary table
icer.n <- tc.drd - tc.vrd
icer.d <- tqaly.drd - tqaly.vrd
icer <- icer.n/icer.d
result <- matrix(c(tc.drd, tly.drd, tqaly.drd, icer, tc.vrd, tly.vrd, tqaly.vrd, NA), ncol = 4, byrow = TRUE)
colnames(result) <- c("Total Costs", "Total Life Years", "Total QALYs", "ICER")
rownames(result) <- c("DRd", "VRd")
result <- as.data.frame(result) %>% mutate(across(where(is.numeric), ~ round(., 2)))

#print result in readable format
format.data.frame(result, nsmall = 2, big.mark = ",") 


################################################################################
## COMPARISON BY TREATMENT, BEST FITTING ##
################################################################################
#PFS COMPARISON
drdpfs <- drdpfs %>% mutate(Arm = 1)
vrdpfs <- vrdpfs %>% mutate(Arm = 0)
pfs.km <- rbind(drdpfs, vrdpfs)
kmfit.pfs <- survfit(Surv(Time, Event) ~ Arm, data = pfs.km, type = "kaplan-meier")

#PLOTS
plot(km.drdpfs, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = cols[1],
     ylab = "Survival Probability", xlab = "Time (weeks)",  main = "Progression-free survival", conf.int = F)
lines(b.drdpfs,  col = cols[1], t = times, lty = 4, ci = F)
lines(trace.vrd$pfs,  col = cols[4], lty = 4, lwd = 3)
legend("topright", cex = 0.9, col = c(cols[1], cols[4], cols[1]),
       bty = "n", lwd = rep(3, 4), lty = c(4,4,1),
       c("DRd PFS Exponential (best fit)", "VRd PFS (based on ITC HR)", "KM DRd PFS"))

#OS COMPARISON
drdos <- drdos %>% mutate(Arm = 1)
vrdos <- vrdos %>% mutate(Arm = 0)
os.km <- rbind(drdos, vrdos)
kmfit.os <- survfit(Surv(Time, Event) ~ Arm, data = os.km, type = "kaplan-meier")

plot(kmfit.os, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = c(cols[4], cols[1]),
     ylab = "Survival Probability", xlab = "Time (weeks)",  main = "Overall survival", conf.int = F)
lines(trace.drd$os,  col = cols[1], lty = 4, lwd = 2)
lines(trace.vrd$os,  col = cols[4], lty = 4, lwd = 2)
lines(kmfit.os, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = c(cols[4], cols[1]))
legend("topright", cex = 0.9, col = c(cols[1], cols[4], cols[1], cols[4]),
       bty = "n", lwd = rep(3, 4), lty = c(1,1,4,4),
       c("KM DRd OS", "KM VRd OS", "DRd OS Exponential", "VRd OS Piecewise Exponential"))

