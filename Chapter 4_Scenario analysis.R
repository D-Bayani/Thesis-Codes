################################################################################
## SCENARIO ANALYSIS AND VALIDATION ##
################################################################################

################################################################################
## WANING OF TREATMENT EFFECT ##
################################################################################

scen.wte <- trace.drd

#incorporate waning of treatment effect
for(i in 523:nrow(scen.wte)) {
  scen.wte$os[i] <- scen.wte$os[i-1]*(1-p3)
}


scen.wte$preprog                     <- scen.wte$pfs  # probability of remaining progression-free
scen.wte$prog                        <- scen.wte$os - scen.wte$pfs   # estimate the probability of remaining in the progressed state
scen.wte$prog[scen.wte$prog < 0]     <- 0  # in cases where the probability is negative replace with zero
scen.wte$dead                        <- 1 - scen.wte$prog - scen.wte$preprog   # probability of being dead
scen.wte$check                       <- scen.wte$prog + scen.wte$preprog + scen.wte$dead

plot(kmfit.os, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = c(cols[4], cols[1]),
     ylab = "Survival Probability", xlab = "Time (weeks)",  main = "Overall survival", conf.int = F)
lines(scen.wte$os,  col = cols[1], lty = 4, lwd = 2)
lines(trace.vrd$os,  col = cols[4], lty = 4, lwd = 2)
lines(kmfit.os, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = c(cols[4], cols[1]))
legend("topright", cex = 0.9, col = c(cols[1], cols[4], cols[1], cols[4]),
       bty = "n", lwd = rep(3, 4), lty = c(1,1,4,4),
       c("KM DRd OS", "KM VRd OS", "DRd OS Exponential (with waning effect after 10 years)", "VRd OS Piecewise Exponential"))

## DRD; calculate costs, life years, and qalys per cycle
scen.wte$drugcosts    <- scen.wte$preprog*(ifelse(scen.wte$cycle<=24, ifelse(scen.wte$cycle<=8, dc_drdc1_2, dc_drdc3_6), dc_drdc7))
scen.wte$costs        <- scen.wte$drugcosts+(scen.wte$prog*(dc_pdrd))+(scen.wte$dead*(c_hospice))+(scen.wte$preprog*cae_neut*rae_drdneut)+(scen.wte$preprog*cae_anae*rae_drdanae)+(scen.wte$preprog*cae_pneu*rae_drdpneu)
scen.wte$disccosts    <- scen.wte$costs/(1 + dr_cost) ^ scen.wte$year
scen.wte$lifeyears    <- (scen.wte$preprog + scen.wte$prog)/(1 + dr_qaly) ^ scen.wte$year
scen.wte$qalys        <- (scen.wte$preprog*u_pf)+(scen.wte$prog*u_pp)-(scen.wte$preprog*rae_drdneut*u_neut)-(scen.wte$preprog*rae_drdpneu*u_pneu)-(scen.wte$preprog*rae_drdanae*u_anae)
scen.wte$discqalys    <- scen.wte$qalys/(1 + dr_qaly) ^ scen.wte$year


## make top row 0 values
scen.wte[1, c(9:14)] <- 0


## SUMMARY OUTPUTS
## DRD
sce1.tpreprog.drd           <- sum(scen.wte$preprog) #get sum of those pre-progression
sce1.tc.drd                 <- sum(scen.wte$disccosts) #total cost
sce1.tly.drd                <- sum(scen.wte$lifeyears)/52.143 #get total life years
sce1.tqaly.drd              <- sum(scen.wte$discqalys) #total qalys

## VRD
tpreprog.vrd           <- sum(trace.vrd$preprog) #get sum of those pre-progression
tc.vrd                 <- sum(trace.vrd$disccosts) #total cost 
tly.vrd                <- sum(trace.vrd$lifeyears)/52.143 #get total life years
tqaly.vrd              <- sum(trace.vrd$discqalys) #total qalys

## create summary table
icer.n.sce1 <- sce1.tc.drd - tc.vrd
icer.d.sce1 <- sce1.tqaly.drd - tqaly.vrd
icer.sce1 <- icer.n.sce1/icer.d.sce1
result <- matrix(c(tc.drd, tly.drd, tqaly.drd, icer, sce1.tc.drd, sce1.tly.drd, sce1.tqaly.drd, icer.sce1, tc.vrd, tly.vrd, tqaly.vrd, NA), ncol = 4, byrow = TRUE)
colnames(result) <- c("Total Costs", "Total Life Years", "Total QALYs", "ICER")
rownames(result) <- c("DRd", "DRd (with TE waning)", "VRd")
result <- as.data.frame(result) %>% mutate(across(where(is.numeric), ~ round(., 2)))

#print result in readable format
format.data.frame(result, nsmall = 2, big.mark = ",") 




################################################################################
## VALIDATION OF EXTRAPOLATION BASED ON UPDATED DATA CUT ##
################################################################################

##Load dataset
drdpfs60 <- read.csv("filename")
drdos60 <- read.csv("filename")

##convert to weekly cycles
drdpfs60$Time = drdpfs60$Time*4.345
drdos60$Time = drdos60$Time*4.345

##generate KM functions
km.drdpfs60 <- survfit(Surv(Time, Event)~1, data=drdpfs60, type = "kaplan-meier")
km.drdos60 <- survfit(Surv(Time, Event)~1, data=drdos60, type = "kaplan-meier")

#plot PFS curves together
plot(km.drdpfs, col = 2, main = "Progression-free Survival", xlab = "Time (weeks)", ylab = "Survival (%)", lwd = 1, conf.int= F)
lines(km.drdpfs60, col = 1, conf.int= F)
legend("bottomleft", col = c(1,2), bty = "n", lwd = rep(1), c("DRd updated data cut", "DRd baseline"))

#plot OS curves together
plot(km.drdos, col = 2, main = "Overall Survival", xlab = "Time (weeks)", ylab = "Survival (%)", lwd = 1, conf.int= F)
lines(km.drdos60, col = 1, conf.int= F)
legend("bottomleft", col = c(1,2), bty = "n", lwd = rep(1), c("DRd updated data cut", "DRd baseline"))

##Fit functions
fit.drdpfs60 <- fit.pm(time = "Time", event = "Event", data = drdpfs60, label = "DRd PFS updated data cut") 
fit.drdos60 <- fit.pm(time = "Time", event = "Event", data = drdos60, label = "DRd OS updated data cut") 

##display AIC and BIC results
fit.drdpfs60[[10]]
fit.drdos60[[10]]

#select best fitting models
b.drdpfs60 <- fit.drdpfs60$Exponential
plot(km.drdpfs60, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted PFS (UDC) for DRd, Exponential")) 
lines(b.drdpfs60,  col = 2, t = times, lty = 2)

b.drdos60 <- fit.drdos60$Exponential
plot(km.drdos60, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted OS (UDC) for DRd, Exponential")) 
lines(b.drdos60,  col = 2, t = times, lty = 2)


#Create PSM based on new functions
## apply PSM function to selected curve ## CHANGE TO BEST FIT OR MOST CONSERVATIVE
psm.drd60 <- partsa(b.drdttd, b.drdos60, time = times)

#plot matrices
matrix.drd60 <- as.matrix(psm.drd60$matrix)    
matplot(matrix.drd60, type = 'l', lty=1, xlab = "Time (weeks)", ylab = "Proportion of patients", main = "DRd (UDC)")                      
legend("right", v.n, col=1:n.s, lty=rep(1,n.s), bty='n')


################################################################################
## COMPARISON BY DATA CUT, BEST FITTING ##
################################################################################
#PFS COMPARISON
drdpfs <- drdpfs %>% mutate(Arm = 1)
drdpfs60 <- drdpfs60 %>% mutate(Arm = 0)
pfs.km60 <- rbind(drdpfs, drdpfs60)
kmfit.pfs60 <- survfit(Surv(Time, Event) ~ Arm, data = pfs.km60, type = "kaplan-meier")

plot(kmfit.pfs60, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = c(cols[4], cols[1]),
     ylab = "Survival Probability", xlab = "Time (weeks)",  main = "Progression-free survival Comparison", conf.int = F)
lines(b.drdpfs,  col = cols[1], t = times, lty = 4, lwd = 2, ci = F)
lines(b.drdpfs60,  col = cols[4], t = times, lty = 4, lwd = 2, ci = F)
lines(kmfit.pfs60, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = c(cols[4], cols[1]))
legend("topright", cex = 0.9, col = c(cols[1], cols[4], cols[1], cols[4]),
       bty = "n", lwd = rep(3, 4), lty = c(1,1,4,4),
       c("KM DRd PFS baseline", "KM DRd PFS (UDC)", "DRd PFS baseline", "DRd PFS (UDC)"))


#OS COMPARISON
drdos <- drdos %>% mutate(Arm = 1)
drdos60 <- drdos60 %>% mutate(Arm = 0)
os.km60 <- rbind(drdos, drdos60)
kmfit.os60 <- survfit(Surv(Time, Event) ~ Arm, data = os.km60, type = "kaplan-meier")

plot(kmfit.os60, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = c(cols[4], cols[1]),
     ylab = "Survival Probability", xlab = "Time (weeks)",  main = "Overall survival Comparison", conf.int = F)
lines(b.drdos,  col = cols[1], t = times, lty = 4, lwd = 2, ci = F)
lines(b.drdos60,  col = cols[4], t = times, lty = 4, lwd = 2, ci = F)
lines(kmfit.os60, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = c(cols[4], cols[1]))
legend("topright", cex = 0.9, col = c(cols[1], cols[4], cols[1], cols[4]),
       bty = "n", lwd = rep(3, 4), lty = c(1,1,4,4),
       c("KM DRd OS baseline", "KM DRd OS (UDC)", "DRd OS baseline", "DRd OS (UDC)"))



##Redo analysis using UDC
trace.drd <- psm.drd60$complete
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
dc_dara           <- 7242.10
dc_bort           <- 246.67
dc_len            <- 47.61
dc_dex            <- 1.50
dc_drdcycle       <- ((dc_len*21)+(dc_dex*4*4))/4 #Rd cost of DRd per week
dc_vrdcycle       <- ((dc_len*14)+(dc_dex*2*8))/3 #Rd cost of VRd induction per week

#Chemo administration and disease monitoring
dc_admin          <- 184.02 #per visit
dc_mgt            <- 152.64/4 #per visit, divided by 4 to get weekly cost

##DRd (With PAP incorporated)
dc_drdc1_2        <- ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        <- ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          <- ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

##VRd
dc_vrdind         <- ((((dc_bort*1.3*1.6)+dc_admin)*4)/3)+dc_vrdcycle+dc_mgt
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
## USE OF FLATIRON DATA ##
################################################################################

##Load dataset
vrdpfsfi <- read.csv("filename")

##convert to weekly cycles
vrdpfsfi$Time = vrdpfsfi$Time*4.345

##generate KM functions
km.vrdpfsfi <- survfit(Surv(Time, Event)~1, data=vrdpfsfi, type = "kaplan-meier")

#plot PFS curves together
plot(km.vrdpfs, col = 2, main = "Progression-free Survival Comparison", xlab = "Time (weeks)", ylab = "Survival (%)", lwd = 1, conf.int= F)
lines(km.vrdpfsfi, col = 1, conf.int= F)
legend("bottomleft", col = c(1,2), bty = "n", lwd = rep(1), c("VRd from Flatiron", "VRd from SWOG-0777"))

##Fit functions
fit.vrdpfsfi <- fit.pm(time = "Time", event = "Event", data = vrdpfsfi, label = "VRd PFS Flatiron") 

##display AIC and BIC results
fit.vrdpfsfi[[10]]


#select best fitting models
b.vrdpfsfi <- fit.vrdpfsfi$Lognormal
plot(km.vrdpfsfi, ylab = "Survival Probability", xlab = "Time (weeks)",  main = paste ("True vs Fitted PFS (UDC) for DRd, Exponential")) 
lines(b.vrdpfsfi,  col = 2, t = times, lty = 2)


#FLATIRON VS SWOG COMPARISON
vrdpfs <- vrdpfs %>% mutate(Arm = 1)
vrdpfsfi <- vrdpfsfi %>% mutate(Arm = 0)
pfs.kmfi <- rbind(vrdpfs, vrdpfsfi)
kmfit.pfsfi <- survfit(Surv(Time, Event) ~ Arm, data = pfs.kmfi, type = "kaplan-meier")

plot(km.vrdpfsfi, lwd = 3, ylim = c(0,1), xlim = c(0,n.t), col = c(cols[4], cols[1]),
     ylab = "Survival Probability", xlab = "Time (weeks)",  main = "Progression-free survival compariosn", conf.int = F)
lines(trace.vrd$pfs,  col = cols[1], lty = 4, lwd = 2)
lines(b.vrdpfsfi,  col = cols[4], t = times, lty = 4, lwd = 2, ci = F)
legend("topright", cex = 0.9, col = c(cols[4], cols[1], cols[4]),
       bty = "n", lwd = rep(3, 4), lty = c(1,4,4),
       c("KM VRd PFS (Flatiron)", "VRd PFS (based on ITC HR)", "VRd PFS (based on Flatiron data)"))


#Create PSM based on new functions
## apply PSM function to selected curve ## CHANGE TO BEST FIT OR MOST CONSERVATIVE
psm.vrdfi <- partsa(b.vrdpfsfi, b.vrdos, time = times)

#plot matrices
matrix.vrdfi <- as.matrix(psm.vrdfi$matrix)    
matplot(matrix.vrdfi, type = 'l', lty=1, xlab = "Time (weeks)", ylab = "Proportion of patients", main = "VRd Flatiron PFS")                      
legend("right", v.n, col=1:n.s, lty=rep(1,n.s), bty='n')

trace.vrd <- psm.vrdfi$complete
trace.drd <- psm.drd$complete


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
dc_dara           <- 7242.10
dc_bort           <- 246.67
dc_len            <- 47.61
dc_dex            <- 1.50
dc_drdcycle       <- ((dc_len*21)+(dc_dex*4*4))/4 #Rd cost of DRd per week
dc_vrdcycle       <- ((dc_len*14)+(dc_dex*2*8))/3 #Rd cost of VRd induction per week

#Chemo administration and disease monitoring
dc_admin          <- 184.02 #per visit
dc_mgt            <- 152.64/4 #per visit, divided by 4 to get weekly cost

##DRd (With PAP incorporated)
dc_drdc1_2        <- ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        <- ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          <- ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

##VRd
dc_vrdind         <- ((((dc_bort*1.3*1.6)+dc_admin)*4)/3)+dc_vrdcycle+dc_mgt
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
trace.drd$lifeyears    <- (trace.drd$preprog + trace.drd$prog)/(1 + dr_qaly)^ trace.drd$year
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
