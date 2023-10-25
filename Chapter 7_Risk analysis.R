options(scipen=999)

################################################################################
## FUNCTION TO RUN MODEL ##
################################################################################

psm <- function(l.params){
  with(as.list(l.params), {
    
    ##set -up
    cycle                 <- times
    year                  <- cycle/52.143
    pfs                   <- 1
    os                    <- 1
    preprog               <- pfs
    prog                  <- os - pfs          # estimate the probability of remaining in the progressed state
    prog[prog < 0]        <- 0                           # in cases where the probability is negative replace with zero
    dead                  <- 1 - prog - preprog          # probability of being dead
    check                 <- prog + preprog + dead
    trace1                <- data.frame(year = year, cycle = cycle, pfs = pfs, os = os, preprog = preprog, prog = prog, dead = dead, check = check)
    trace2                <- data.frame(year = year, cycle = cycle, pfs = pfs, os = os, preprog = preprog, prog = prog, dead = dead, check = check)
    
    ##efficacy parameters
    #trace 1 (DRd)
    pfs.rate1 = 1-exp(-exp(drd_ttd))
    
    for(i in 2:nrow(trace1)) {
      trace1$pfs[i] <- trace1$pfs[i-1]*(1-pfs.rate1)
    }
    
    os.rate1 = 1-exp(-exp(drd_os))
    
    for(i in 2:nrow(trace1)) {
      trace1$os[i] <- trace1$os[i-1]*(1-os.rate1)
    }
    
    
    #trace 2 (VRd)
    pfs.rate2 = 1-exp(-exp(drd_pfs + log(1+(1-hr))))
    
    for(i in 2:nrow(trace2)) {
      trace2$pfs[i] <- trace2$pfs[i-1]*(1-pfs.rate2)
    }
    
    p.1 = 1-exp(-exp(vrd_os1+vrd_os2))
    p.2 = 1-exp(-exp(vrd_os1+vrd_os3))
    p.3 = 1-exp(-exp(vrd_os1+vrd_os4))
    
    for(i in 2:nrow(trace2)) {
      trace2$os[i] <- trace2$os[i-1]*(1-p.1)
    }
    
    for(i in 14:nrow(trace2)) {
      trace2$os[i] <- trace2$os[i-1]*(1-p.2)
    }
    
    for(i in 163:nrow(trace2)) {
      trace2$os[i] <- trace2$os[i-1]*(1-p.3)
    }
    
    
    
    trace1$preprog                     <- trace1$pfs  # probability of remaining progression-free
    trace1$prog                        <- trace1$os - trace1$pfs   # estimate the probability of remaining in the progressed state
    trace1$prog[trace1$prog < 0]    <- 0  # in cases where the probability is negative replace with zero
    trace1$dead                        <- 1 - trace1$prog - trace1$preprog   # probability of being dead
    trace1$check                       <- trace1$prog + trace1$preprog + trace1$dead
    
    trace2$preprog                     <- trace2$pfs  # probability of remaining progression-free
    trace2$prog                        <- trace2$os - trace2$pfs   # estimate the probability of remaining in the progressed state
    trace2$prog[trace2$prog < 0]    <- 0  # in cases where the probability is negative replace with zero
    trace2$dead                        <- 1 - trace2$prog - trace2$preprog   # probability of being dead
    trace2$check                       <- trace2$prog + trace2$preprog + trace2$dead
    
    
    trace1$drugcosts    <- trace1$preprog*(ifelse(trace1$cycle<=24, ifelse(trace1$cycle<=8, dc_drdc1_2, dc_drdc3_6), dc_drdc7))
    trace1$costs        <- trace1$drugcosts+(trace1$prog*(dc_pdrd))+(trace1$dead*(c_hospice))+(trace1$preprog*cae_neut*rae_drdneut)+(trace1$preprog*cae_anae*rae_drdanae)+(trace1$preprog*cae_pneu*rae_drdpneu)
    trace1$disccosts    <- trace1$costs/(1 + dr_cost) ^ trace1$year
    trace1$lifeyears    <- trace1$preprog + trace1$prog 
    trace1$qalys        <- (trace1$preprog*u_pf)+(trace1$prog*u_pp)-(trace1$preprog*rae_drdneut*u_neut)-(trace1$preprog*rae_drdpneu*u_pneu)-(trace1$preprog*rae_drdanae*u_anae)
    trace1$discqalys    <- trace1$qalys/(1 + dr_qaly) ^ trace1$year
    
    
    trace2$drugcosts    <- trace2$preprog*(ifelse(trace2$cycle<=24, dc_vrdind, dc_vrdmaint))
    trace2$costs        <- trace2$drugcosts+(trace2$prog*(dc_pvrd))+(trace2$dead*(c_hospice))+(trace2$preprog*cae_neut*rae_vrdneut)+(trace2$preprog*cae_anae*rae_vrdanae)
    trace2$disccosts    <- trace2$costs/(1 + dr_cost) ^ trace2$year
    trace2$lifeyears    <- trace2$preprog + trace2$prog 
    trace2$qalys        <- (trace2$preprog*u_pf)+(trace2$prog*u_pp)-(trace2$preprog*rae_vrdneut*u_neut)-(trace2$preprog*rae_vrdanae*u_anae)-(trace2$preprog*rae_vrddiar*u_diar)-(trace2$preprog*rae_vrdneuro*u_neuro)
    trace2$discqalys    <- trace2$qalys/(1 + dr_qaly) ^ trace2$year
    
    trace1[1, c(9:14)] <- 0
    trace2[1, c(9:14)] <- 0
    
    t1preprog           <- sum(trace1$preprog) #get sum of those pre-progression
    t1c                 <- sum(trace1$disccosts) #total cost
    t1ly                <- sum(trace1$lifeyears)/52.143 #get total life years
    t1qaly              <- sum(trace1$discqalys) #total qalys
    
    t2preprog           <- sum(trace2$preprog) #get sum of those pre-progression
    t2c                 <- sum(trace2$disccosts) #total cost
    t2ly                <- sum(trace2$lifeyears)/52.143 #get total life years
    t2qaly              <- sum(trace2$discqalys) #total qalys
    
    ## result
    icer.n <- t1c - t2c
    icer.d <- t1qaly - t2qaly
    icer <- icer.n / icer.d
    result <- matrix(c(t1c, t1qaly, icer, t2c, t2qaly, NA), ncol = 3, byrow = TRUE)
    colnames(result) <- c("Total Costs", "Total QALYs", "ICER")
    rownames(result) <- c("DRd", "VRd")
    result <- as.data.frame(result) %>% mutate(across(where(is.numeric), ~ round(., 2)))
    
    #print result in readable format
    format.data.frame(result, nsmall = 2, big.mark = ",") 
    
    data <- data.frame(l.params, t1c, t2c, t1qaly, t2qaly)
    out  <- list(drd.trace = trace1, vrd.trace = trace2, result, data)
    
    
    return(out)
  }
  )
}



################################################################################
## FUNCTION TO DRAW RANDOM VALUES ##
################################################################################
draw <- function() {
  p.params <-  list(
    
    #efficacy 
    drd_ttd           = rnorm(1, mean = b.drdttd$res.t[1], sd = b.drdttd$res.t[4]),
    drd_pfs           = rnorm(1, mean = b.drdpfs$res.t[1], sd = b.drdttd$res.t[4]),
    drd_os            = rnorm(1, mean = b.drdos$res.t[1], sd = b.drdos$res.t[4]),
    
    vrd_os1            = pw.exp$res.t[1],
    vrd_os2            = pw.exp$res.t[2],
    vrd_os3            = pw.exp$res.t[3],
    vrd_os4            = pw.exp$res.t[4],
    
    #hazards ratio
    hr                = rlnorm(1, meanlog = -0.3270, sdlog = 0.1819),
    
    #Unit costs from hospital (minimum value)
    dc_dara           = dc_dara,
    dc_bort           = dc_bort,
    dc_len            = dc_len,
    dc_dex            = dc_dex,
    dc_drdcycle       = dc_drdcycle, #Rd cost of DRd per week
    dc_vrdcycle       = dc_vrdcycle, #Rd cost of VRd induction per week
    
    #Chemo administration and disease monitoring
    dc_admin          = dc_admin, #per visit
    dc_mgt            = dc_mgt, #per visit, divided by 4 to get weekly cost
    
    ##DRd
    dc_drdc1_2        = dc_drdc1_2,
    dc_drdc3_6        = dc_drdc3_6,
    dc_drdc7          = dc_drdc7,
    
    ##VRd
    dc_vrdind         = dc_vrdind,
    dc_vrdmaint       = dc_vrdmaint, #Rd maintenance
    
    
    #Cost of subsequent treatments
    dc_pdrd           = (rgamma(1, shape = 1.681341814, scale = 988.3118272))/4,
    dc_pvrd           = (rgamma(1, shape = 1.538700343, scale = 2786.994894))/4,
    
    ## HOSPICE
    c_hospice         = (rgamma(1, shape = 56.57942829, scale = 154.6498483))/4,
    
    ## HOSPITALIZATION FOR AE
    cae_neut          = rgamma(1, shape = 0.287901009, scale = 46079.37995),
    cae_anae          = rgamma(1, shape = 0.924930001, scale = 2775.874929),
    cae_pneu          = rgamma(1, shape = 0.348637828, scale = 39470.58774),
    
    ##### UTILITIES
    ## STATES per cycle (divided by weeks)
    u_pf              = (rbeta(1, shape1 = 2.394768629, shape2 = 0.908360515)/52.143),
    u_pp              = (rbeta(1, shape1 = 157.1490556, shape2 = 84.61872222)/52.143),
    
    ## AE DISUTILITY (NEGATIVE VALUES) per event
    u_neut            = 0.145*rgamma(1, shape = 0.706176942, scale = 13.863381)/365.25,
    u_anae            = 0.25*rgamma(1, shape = 0.483498413, scale = 4.01242268)/365.25,
    u_pneu            = 0.14*rgamma(1, shape = 0.553766629, scale = 20.11677738)/365.25,
    u_diar            = 0.103*rgamma(1, shape = 1.106358662, scale = 4.401827515)/365.25,
    u_neuro           = 0.065*rgamma(1, shape = 12.39833366, scale = 0.20164)/365.25,
    
    #### AE RISK
    rae_drdneut       = (rbeta(1, shape1 =  0.319, shape2 =  99.68)),
    rae_drdpneu       = (rbeta(1, shape1 =  0.087, shape2 =  99.91)),
    rae_drdanae       = (rbeta(1, shape1 =  0.075, shape2 =  99.92)),
    
    rae_vrdneut       = (rbeta(1, shape1 =  0.043, shape2 =  99.96)),
    rae_vrdanae       = (rbeta(1, shape1 =  0.176, shape2 =  99.82)),
    rae_vrddiar       = (rbeta(1, shape1 =  0.068, shape2 =  99.93)),
    rae_vrdneuro      = (rbeta(1, shape1 =  0.111, shape2 =  99.89))
  )
  
  params <- data.frame(p.params)
  return(params)
}


draw()

################################################################################
## FUNCTION TO LOOP BASED ON NUMBER OF TRIALS ##
################################################################################

psa <- function(trials = trials) {
  trials = trials
  psa.result <- vector("list", trials)
  for (i in 1:trials) {
    psa.result[[i]] <- psm(draw())[[4]]
  }
  psa.result <- do.call(rbind, psa.result)
  psa.result <- psa.result[, -c(9:21)]
  return(psa.result)
}


###RUN PSA FOR SCENARIOS (S1: Base case)
psa.s1 <- psa(trials = 1000)

s1.par <- psa.s1[, c(1:3, 8:28)]
s1.cos <- psa.s1[, c(29, 30)]
s1.qal <- psa.s1[, c(31, 32)]

colnames(s1.cos) <- c("DRd", "VRd")
colnames(s1.qal) <- c("DRd", "VRd")

write.csv(s1.par, "s1.par.csv", row.names = FALSE)
write.csv(s1.cos, "s1.costs.csv", row.names = FALSE)
write.csv(s1.qal, "s1.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (S2) 30%
dc_dara           = 7242.10*0.7
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

psa.s2 <- psa(trials = 1000)

s2.par <- psa.s2[, c(1:3, 8:28)]
s2.cos <- psa.s2[, c(29, 30)]
s2.qal <- psa.s2[, c(31, 32)]

colnames(s2.cos) <- c("DRd", "VRd")
colnames(s2.qal) <- c("DRd", "VRd")

write.csv(s2.par, "s2.par.csv", row.names = FALSE)
write.csv(s2.cos, "s2.costs.csv", row.names = FALSE)
write.csv(s2.qal, "s2.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (D1) - 5%
dc_dara           = 7242.10*0.95
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

psa.d1 <- psa(trials = 1000)

d1.par <- psa.d1[, c(1:3, 8:28)]
d1.cos <- psa.d1[, c(29, 30)]
d1.qal <- psa.d1[, c(31, 32)]

colnames(d1.cos) <- c("DRd", "VRd")
colnames(d1.qal) <- c("DRd", "VRd")

write.csv(d1.par, "d1.par.csv", row.names = FALSE)
write.csv(d1.cos, "d1.costs.csv", row.names = FALSE)
write.csv(d1.qal, "d1.qalys.csv", row.names = FALSE)

###RUN PSA FOR PRICE DISCOUNT (D2) - 10%
dc_dara           = 7242.10*0.90
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

psa.d2 <- psa(trials = 1000)

d2.par <- psa.d2[, c(1:3, 8:28)]
d2.cos <- psa.d2[, c(29, 30)]
d2.qal <- psa.d2[, c(31, 32)]

colnames(d2.cos) <- c("DRd", "VRd")
colnames(d2.qal) <- c("DRd", "VRd")

write.csv(d2.par, "d2.par.csv", row.names = FALSE)
write.csv(d2.cos, "d2.costs.csv", row.names = FALSE)
write.csv(d2.qal, "d2.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (D3) - 15%
dc_dara           = 7242.10*0.85
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt


psa.d3 <- psa(trials = 1000)

d3.par <- psa.d3[, c(1:3, 8:28)]
d3.cos <- psa.d3[, c(29, 30)]
d3.qal <- psa.d3[, c(31, 32)]

colnames(d3.cos) <- c("DRd", "VRd")
colnames(d3.qal) <- c("DRd", "VRd")

write.csv(d3.par, "d3.par.csv", row.names = FALSE)
write.csv(d3.cos, "d3.costs.csv", row.names = FALSE)
write.csv(d3.qal, "d3.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (D4) - 20%
dc_dara           = 7242.10*0.80
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt


psa.d4 <- psa(trials = 1000)

d4.par <- psa.d4[, c(1:3, 8:28)]
d4.cos <- psa.d4[, c(29, 30)]
d4.qal <- psa.d4[, c(31, 32)]

colnames(d4.cos) <- c("DRd", "VRd")
colnames(d4.qal) <- c("DRd", "VRd")

write.csv(d4.par, "d4.par.csv", row.names = FALSE)
write.csv(d4.cos, "d4.costs.csv", row.names = FALSE)
write.csv(d4.qal, "d4.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (D5) - 25%
dc_dara           = 7242.10*0.75
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

psa.d5 <- psa(trials = 1000)

d5.par <- psa.d5[, c(1:3, 8:28)]
d5.cos <- psa.d5[, c(29, 30)]
d5.qal <- psa.d5[, c(31, 32)]

colnames(d5.cos) <- c("DRd", "VRd")
colnames(d5.qal) <- c("DRd", "VRd")

write.csv(d5.par, "d5.par.csv", row.names = FALSE)
write.csv(d5.cos, "d5.costs.csv", row.names = FALSE)
write.csv(d5.qal, "d5.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (D7) - 35%
dc_dara           = 7242.10*0.65
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

psa.d7 <- psa(trials = 1000)

d7.par <- psa.d7[, c(1:3, 8:28)]
d7.cos <- psa.d7[, c(29, 30)]
d7.qal <- psa.d7[, c(31, 32)]

colnames(d7.cos) <- c("DRd", "VRd")
colnames(d7.qal) <- c("DRd", "VRd")

write.csv(d7.par, "d7.par.csv", row.names = FALSE)
write.csv(d7.cos, "d7.costs.csv", row.names = FALSE)
write.csv(d7.qal, "d7.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (D8) - 40%
dc_dara           = 7242.10*0.6
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

psa.d8 <- psa(trials = 1000)

d8.par <- psa.d8[, c(1:3, 8:28)]
d8.cos <- psa.d8[, c(29, 30)]
d8.qal <- psa.d8[, c(31, 32)]

colnames(d8.cos) <- c("DRd", "VRd")
colnames(d8.qal) <- c("DRd", "VRd")

write.csv(d8.par, "d8.par.csv", row.names = FALSE)
write.csv(d8.cos, "d8.costs.csv", row.names = FALSE)
write.csv(d8.qal, "d8.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (D9) - 45%
dc_dara           = 7242.10*0.55
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

psa.d9 <- psa(trials = 1000)

d9.par <- psa.d9[, c(1:3, 8:28)]
d9.cos <- psa.d9[, c(29, 30)]
d9.qal <- psa.d9[, c(31, 32)]

colnames(d9.cos) <- c("DRd", "VRd")
colnames(d9.qal) <- c("DRd", "VRd")

write.csv(d9.par, "d9.par.csv", row.names = FALSE)
write.csv(d9.cos, "d9.costs.csv", row.names = FALSE)
write.csv(d9.qal, "d9.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (D10) - 50%
dc_dara           = 7242.10*0.5
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

psa.d10 <- psa(trials = 1000)

d10.par <- psa.d10[, c(1:3, 8:28)]
d10.cos <- psa.d10[, c(29, 30)]
d10.qal <- psa.d10[, c(31, 32)]

colnames(d10.cos) <- c("DRd", "VRd")
colnames(d10.qal) <- c("DRd", "VRd")

write.csv(d10.par, "d10.par.csv", row.names = FALSE)
write.csv(d10.cos, "d10.costs.csv", row.names = FALSE)
write.csv(d10.qal, "d10.qalys.csv", row.names = FALSE)


###RUN PSA FOR PRICE DISCOUNT (D11) - 55%
dc_dara           = 7242.10*0.45
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt

psa.d11 <- psa(trials = 1000)

d11.par <- psa.d11[, c(1:3, 8:28)]
d11.cos <- psa.d11[, c(29, 30)]
d11.qal <- psa.d11[, c(31, 32)]

colnames(d11.cos) <- c("DRd", "VRd")
colnames(d11.qal) <- c("DRd", "VRd")

write.csv(d11.par, "d11.par.csv", row.names = FALSE)
write.csv(d11.cos, "d11.costs.csv", row.names = FALSE)
write.csv(d11.qal, "d11.qalys.csv", row.names = FALSE)


###RUN PSA FOR PERFECT EVIDENCE (MAIA) (S3)
dc_dara           = 7242.10
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt


#create function to draw scenario 3 parameters only
draw.s3 <- function() {
  p.params <-  list(
    
    #efficacy 
    drd_ttd           = b.drdttd$res.t[1],
    drd_pfs           = b.drdpfs$res.t[1],
    drd_os            = b.drdos$res.t[1],
    
    vrd_os1            = pw.exp$res.t[1],
    vrd_os2            = pw.exp$res.t[2],
    vrd_os3            = pw.exp$res.t[3],
    vrd_os4            = pw.exp$res.t[4],
    
    #hazards ratio
    hr                = rlnorm(1, meanlog = -0.3270, sdlog = 0.1819),
    
    #Unit costs from hospital (minimum value)
    dc_dara           = dc_dara,
    dc_bort           = dc_bort,
    dc_len            = dc_len,
    dc_dex            = dc_dex,
    dc_drdcycle       = dc_drdcycle, #Rd cost of DRd per week
    dc_vrdcycle       = dc_vrdcycle, #Rd cost of VRd induction per week
    
    #Chemo administration and disease monitoring
    dc_admin          = dc_admin, #per visit
    dc_mgt            = dc_mgt, #per visit, divided by 4 to get weekly cost
    
    ##DRd
    dc_drdc1_2        = dc_drdc1_2,
    dc_drdc3_6        = dc_drdc3_6,
    dc_drdc7          = dc_drdc7,
    
    ##VRd
    dc_vrdind         = dc_vrdind,
    dc_vrdmaint       = dc_vrdmaint, #Rd maintenance
    
    
    #Cost of subsequent treatments
    dc_pdrd           = (rgamma(1, shape = 1.681341814, scale = 988.3118272))/4,
    dc_pvrd           = (rgamma(1, shape = 1.538700343, scale = 2786.994894))/4,
    
    ## HOSPICE
    c_hospice         = (rgamma(1, shape = 56.57942829, scale = 154.6498483))/4,
    
    ## HOSPITALIZATION FOR AE
    cae_neut          = rgamma(1, shape = 0.287901009, scale = 46079.37995),
    cae_anae          = rgamma(1, shape = 0.924930001, scale = 2775.874929),
    cae_pneu          = rgamma(1, shape = 0.348637828, scale = 39470.58774),
    
    ##### UTILITIES
    ## STATES per cycle (divided by weeks)
    u_pf              = (rbeta(1, shape1 = 2.394768629, shape2 = 0.908360515)/52.143),
    u_pp              = (rbeta(1, shape1 = 157.1490556, shape2 = 84.61872222)/52.143),
    
    ## AE DISUTILITY (NEGATIVE VALUES) per event
    u_neut            = 0.145*rgamma(1, shape = 0.706176942, scale = 13.863381)/365.25,
    u_anae            = 0.25*rgamma(1, shape = 0.483498413, scale = 4.01242268)/365.25,
    u_pneu            = 0.14*rgamma(1, shape = 0.553766629, scale = 20.11677738)/365.25,
    u_diar            = 0.103*rgamma(1, shape = 1.106358662, scale = 4.401827515)/365.25,
    u_neuro           = 0.065*rgamma(1, shape = 12.39833366, scale = 0.20164)/365.25,
    
    #### AE RISK
    rae_drdneut       = (rbeta(1, shape1 =  0.319, shape2 =  99.68)),
    rae_drdpneu       = (rbeta(1, shape1 =  0.087, shape2 =  99.91)),
    rae_drdanae       = (rbeta(1, shape1 =  0.075, shape2 =  99.92)),
    
    rae_vrdneut       = (rbeta(1, shape1 =  0.043, shape2 =  99.96)),
    rae_vrdanae       = (rbeta(1, shape1 =  0.176, shape2 =  99.82)),
    rae_vrddiar       = (rbeta(1, shape1 =  0.068, shape2 =  99.93)),
    rae_vrdneuro      = (rbeta(1, shape1 =  0.111, shape2 =  99.89))
  )
  
  params <- data.frame(p.params)
  return(params)
}

draw.s3()

##create function for psa for s3 only
psa.s3 <- function(trials = trials) {
  trials = trials
  psa.result <- vector("list", trials)
  for (i in 1:trials) {
    psa.result[[i]] <- psm(draw.s3())[[4]]
  }
  psa.result <- do.call(rbind, psa.result)
  psa.result <- psa.result[, -c(1:7,9:21)]
  return(psa.result)
}

psa.s3 <- psa.s3(trials = 1000)

s3.par <- psa.s3[, c(1:21)]
s3.cos <- psa.s3[, c(22, 23)]
s3.qal <- psa.s3[, c(24, 25)]

colnames(s3.cos) <- c("DRd", "VRd")
colnames(s3.qal) <- c("DRd", "VRd")

write.csv(s3.par, "s3.par.csv", row.names = FALSE)
write.csv(s3.cos, "s3.costs.csv", row.names = FALSE)
write.csv(s3.qal, "s3.qalys.csv", row.names = FALSE)


###RUN PSA ASSUMING NO INDIRECT EVIDENCE, NO ITC (S4)
dc_dara           = 7242.10
dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt
dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt


#create function to draw scenario 4 parameters only
draw.s4 <- function() {
  p.params <-  list(
    
    #efficacy 
    drd_ttd           = rnorm(1, mean = b.drdttd$res.t[1], sd = b.drdttd$res.t[4]),
    drd_pfs           = rnorm(1, mean = b.drdpfs$res.t[1], sd = b.drdttd$res.t[4]),
    drd_os            = rnorm(1, mean = b.drdos$res.t[1], sd = b.drdos$res.t[4]),
    
    vrd_os1            = pw.exp$res.t[1],
    vrd_os2            = pw.exp$res.t[2],
    vrd_os3            = pw.exp$res.t[3],
    vrd_os4            = pw.exp$res.t[4],
    
    #hazards ratio
    hr                = 0.7211,
    
    #Unit costs from hospital (minimum value)
    dc_dara           = dc_dara,
    dc_bort           = dc_bort,
    dc_len            = dc_len,
    dc_dex            = dc_dex,
    dc_drdcycle       = dc_drdcycle, #Rd cost of DRd per week
    dc_vrdcycle       = dc_vrdcycle, #Rd cost of VRd induction per week
    
    #Chemo administration and disease monitoring
    dc_admin          = dc_admin, #per visit
    dc_mgt            = dc_mgt, #per visit, divided by 4 to get weekly cost
    
    ##DRd
    dc_drdc1_2        = dc_drdc1_2,
    dc_drdc3_6        = dc_drdc3_6,
    dc_drdc7          = dc_drdc7,
    
    ##VRd
    dc_vrdind         = dc_vrdind,
    dc_vrdmaint       = dc_vrdmaint, #Rd maintenance
    
    
    #Cost of subsequent treatments
    dc_pdrd           = (rgamma(1, shape = 1.681341814, scale = 988.3118272))/4,
    dc_pvrd           = (rgamma(1, shape = 1.538700343, scale = 2786.994894))/4,
    
    ## HOSPICE
    c_hospice         = (rgamma(1, shape = 56.57942829, scale = 154.6498483))/4,
    
    ## HOSPITALIZATION FOR AE
    cae_neut          = rgamma(1, shape = 0.287901009, scale = 46079.37995),
    cae_anae          = rgamma(1, shape = 0.924930001, scale = 2775.874929),
    cae_pneu          = rgamma(1, shape = 0.348637828, scale = 39470.58774),
    
    ##### UTILITIES
    ## STATES per cycle (divided by weeks)
    u_pf              = (rbeta(1, shape1 = 2.394768629, shape2 = 0.908360515)/52.143),
    u_pp              = (rbeta(1, shape1 = 157.1490556, shape2 = 84.61872222)/52.143),
    
    ## AE DISUTILITY (NEGATIVE VALUES) per event
    u_neut            = 0.145*rgamma(1, shape = 0.706176942, scale = 13.863381)/365.25,
    u_anae            = 0.25*rgamma(1, shape = 0.483498413, scale = 4.01242268)/365.25,
    u_pneu            = 0.14*rgamma(1, shape = 0.553766629, scale = 20.11677738)/365.25,
    u_diar            = 0.103*rgamma(1, shape = 1.106358662, scale = 4.401827515)/365.25,
    u_neuro           = 0.065*rgamma(1, shape = 12.39833366, scale = 0.20164)/365.25,
    
    #### AE RISK
    rae_drdneut       = (rbeta(1, shape1 =  0.319, shape2 =  99.68)),
    rae_drdpneu       = (rbeta(1, shape1 =  0.087, shape2 =  99.91)),
    rae_drdanae       = (rbeta(1, shape1 =  0.075, shape2 =  99.92)),
    
    rae_vrdneut       = (rbeta(1, shape1 =  0.043, shape2 =  99.96)),
    rae_vrdanae       = (rbeta(1, shape1 =  0.176, shape2 =  99.82)),
    rae_vrddiar       = (rbeta(1, shape1 =  0.068, shape2 =  99.93)),
    rae_vrdneuro      = (rbeta(1, shape1 =  0.111, shape2 =  99.89))
  )
  
  params <- data.frame(p.params)
  return(params)
}

draw.s4()

##create function for psa for s3 only
psa.s4 <- function(trials = trials) {
  trials = trials
  psa.result <- vector("list", trials)
  for (i in 1:trials) {
    psa.result[[i]] <- psm(draw.s4())[[4]]
  }
  psa.result <- do.call(rbind, psa.result)
  psa.result <- psa.result[, -c(4:21)]
  return(psa.result)
}

psa.s4 <- psa.s4(trials = 1000)

s4.par <- psa.s4[, c(1:23)]
s4.cos <- psa.s4[, c(24, 25)]
s4.qal <- psa.s4[, c(26, 27)]

colnames(s4.cos) <- c("DRd", "VRd")
colnames(s4.qal) <- c("DRd", "VRd")

write.csv(s4.par, "s4.par.csv", row.names = FALSE)
write.csv(s4.cos, "s4.costs.csv", row.names = FALSE)
write.csv(s4.qal, "s4.qalys.csv", row.names = FALSE)



###RUN PSA WITH CONDITIONAL TREATMENT CONTINUATION (S5)

psm.s5 <- function(l.params){
  with(as.list(l.params), {
    
    ##set -up
    cycle                 <- times
    year                  <- cycle/52.143
    pfs                   <- 1
    os                    <- 1
    preprog               <- pfs
    prog                  <- os - pfs          # estimate the probability of remaining in the progressed state
    prog[prog < 0]        <- 0                           # in cases where the probability is negative replace with zero
    dead                  <- 1 - prog - preprog          # probability of being dead
    check                 <- prog + preprog + dead
    trace1                <- data.frame(year = year, cycle = cycle, pfs = pfs, os = os, preprog = preprog, prog = prog, dead = dead, check = check)
    trace2                <- data.frame(year = year, cycle = cycle, pfs = pfs, os = os, preprog = preprog, prog = prog, dead = dead, check = check)
    
    ##efficacy parameters
    #trace 1 (DRd)
    pfs.rate1 = 1-exp(-exp(drd_ttd))
    
    for(i in 2:nrow(trace1)) {
      trace1$pfs[i] <- trace1$pfs[i-1]*(1-pfs.rate1)
    }
    
    os.rate1 = 1-exp(-exp(drd_os))
    
    for(i in 2:nrow(trace1)) {
      trace1$os[i] <- trace1$os[i-1]*(1-os.rate1)
    }
    
    
    #trace 2 (VRd)
    pfs.rate2 = 1-exp(-exp(drd_pfs + log(1+(1-hr))))
    
    for(i in 2:nrow(trace2)) {
      trace2$pfs[i] <- trace2$pfs[i-1]*(1-pfs.rate2)
    }
    
    p.1 = 1-exp(-exp(vrd_os1+vrd_os2))
    p.2 = 1-exp(-exp(vrd_os1+vrd_os3))
    p.3 = 1-exp(-exp(vrd_os1+vrd_os4))
    
    for(i in 2:nrow(trace2)) {
      trace2$os[i] <- trace2$os[i-1]*(1-p.1)
    }
    
    for(i in 14:nrow(trace2)) {
      trace2$os[i] <- trace2$os[i-1]*(1-p.2)
    }
    
    for(i in 163:nrow(trace2)) {
      trace2$os[i] <- trace2$os[i-1]*(1-p.3)
    }
    
    
    
    trace1$preprog                     <- trace1$pfs  # probability of remaining progression-free
    trace1$prog                        <- trace1$os - trace1$pfs   # estimate the probability of remaining in the progressed state
    trace1$prog[trace1$prog < 0]    <- 0  # in cases where the probability is negative replace with zero
    trace1$dead                        <- 1 - trace1$prog - trace1$preprog   # probability of being dead
    trace1$check                       <- trace1$prog + trace1$preprog + trace1$dead
    
    trace2$preprog                     <- trace2$pfs  # probability of remaining progression-free
    trace2$prog                        <- trace2$os - trace2$pfs   # estimate the probability of remaining in the progressed state
    trace2$prog[trace2$prog < 0]    <- 0  # in cases where the probability is negative replace with zero
    trace2$dead                        <- 1 - trace2$prog - trace2$preprog   # probability of being dead
    trace2$check                       <- trace2$prog + trace2$preprog + trace2$dead
    
    
    trace1$drugcosts    <- trace1$preprog*(ifelse(trace1$cycle<=13, dc_drdc1_2, dc_drdc1_2*0.932))
    trace1$costs        <- trace1$drugcosts+(trace1$prog*(dc_pdrd))+(trace1$dead*(c_hospice))+(trace1$preprog*cae_neut*rae_drdneut)+(trace1$preprog*cae_anae*rae_drdanae)+(trace1$preprog*cae_pneu*rae_drdpneu)
    trace1$disccosts    <- trace1$costs/(1 + dr_cost) ^ trace1$year
    trace1$lifeyears    <- trace1$preprog + trace1$prog 
    trace1$qalys        <- (trace1$preprog*u_pf)+(trace1$prog*u_pp)-(trace1$preprog*rae_drdneut*u_neut)-(trace1$preprog*rae_drdpneu*u_pneu)-(trace1$preprog*rae_drdanae*u_anae)
    trace1$discqalys    <- trace1$qalys/(1 + dr_qaly) ^ trace1$year
    
    
    trace2$drugcosts    <- trace2$preprog*(ifelse(trace2$cycle<=24, dc_vrdind, dc_vrdmaint))
    trace2$costs        <- trace2$drugcosts+(trace2$prog*(dc_pvrd))+(trace2$dead*(c_hospice))+(trace2$preprog*cae_neut*rae_vrdneut)+(trace2$preprog*cae_anae*rae_vrdanae)
    trace2$disccosts    <- trace2$costs/(1 + dr_cost) ^ trace2$year
    trace2$lifeyears    <- trace2$preprog + trace2$prog 
    trace2$qalys        <- (trace2$preprog*u_pf)+(trace2$prog*u_pp)-(trace2$preprog*rae_vrdneut*u_neut)-(trace2$preprog*rae_vrdanae*u_anae)-(trace2$preprog*rae_vrddiar*u_diar)-(trace2$preprog*rae_vrdneuro*u_neuro)
    trace2$discqalys    <- trace2$qalys/(1 + dr_qaly) ^ trace2$year
    
    trace1[1, c(9:14)] <- 0
    trace2[1, c(9:14)] <- 0
    
    t1preprog           <- sum(trace1$preprog) #get sum of those pre-progression
    t1c                 <- sum(trace1$disccosts) #total cost
    t1ly                <- sum(trace1$lifeyears)/52.143 #get total life years
    t1qaly              <- sum(trace1$discqalys) #total qalys
    
    t2preprog           <- sum(trace2$preprog) #get sum of those pre-progression
    t2c                 <- sum(trace2$disccosts) #total cost
    t2ly                <- sum(trace2$lifeyears)/52.143 #get total life years
    t2qaly              <- sum(trace2$discqalys) #total qalys
    
    ## result
    icer.n <- t1c - t2c
    icer.d <- t1qaly - t2qaly
    icer <- icer.n / icer.d
    result <- matrix(c(t1c, t1qaly, icer, t2c, t2qaly, NA), ncol = 3, byrow = TRUE)
    colnames(result) <- c("Total Costs", "Total QALYs", "ICER")
    rownames(result) <- c("DRd", "VRd")
    result <- as.data.frame(result) %>% mutate(across(where(is.numeric), ~ round(., 2)))
    
    #print result in readable format
    format.data.frame(result, nsmall = 2, big.mark = ",") 
    
    data <- data.frame(l.params, t1c, t2c, t1qaly, t2qaly)
    out  <- list(drd.trace = trace1, vrd.trace = trace2, result, data)
    
    
    return(out)
  }
  )
}

psa.s5 <- function(trials = trials) {
  trials = trials
  psa.result <- vector("list", trials)
  for (i in 1:trials) {
    psa.result[[i]] <- psm.s5(draw())[[4]]
  }
  psa.result <- do.call(rbind, psa.result)
  psa.result <- psa.result[, -c(9:21)]
  return(psa.result)
}


psa.s5 <- psa.s5(trials = 1000)

s5.par <- psa.s5[, c(1:3, 8:28)]
s5.cos <- psa.s5[, c(29, 30)]
s5.qal <- psa.s5[, c(31, 32)]

colnames(s5.cos) <- c("DRd", "VRd")
colnames(s5.qal) <- c("DRd", "VRd")

write.csv(s5.par, "s5.par.csv", row.names = FALSE)
write.csv(s5.cos, "s5.costs.csv", row.names = FALSE)
write.csv(s5.qal, "s5.qalys.csv", row.names = FALSE)

################################################################################
## LIST DETERMINISTIC PARAMETERS ##
################################################################################
d.params <- list(
  #rate
  drd_ttd           = b.drdttd$res.t[1],
  drd_pfs           = b.drdpfs$res.t[1],
  drd_os            = b.drdos$res.t[1],
  
  vrd_os1            = pw.exp$res.t[1],
  vrd_os2            = pw.exp$res.t[2],
  vrd_os3            = pw.exp$res.t[3],
  vrd_os4            = pw.exp$res.t[4],
  
  
  #hazards ratio
  hr                = 0.7211,
  
  dc_dara           = 7242.10,
  dc_bort           = 246.67,
  dc_len            = 47.61,
  dc_dex            = 1.50,
  dc_drdcycle       = ((dc_len*21)+(dc_dex*4*4))/4, #Rd cost of DRd per week
  dc_vrdcycle       = ((dc_len*14)+(dc_dex*2*8))/3, #Rd cost of VRd induction per week
  
  #Chemo administration and disease monitoring
  dc_admin          = 184.02, #per visit
  dc_mgt            = 152.64/4, #per visit, divided by 4 to get weekly cost
  
  ##DRd
  dc_drdc1_2        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt,
  dc_drdc3_6        = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt,
  dc_drdc7          = ((dc_dara+dc_admin)/4)+dc_drdcycle+dc_mgt,
  
  ##VRd
  dc_vrdind         = ((((dc_bort*1.3*1.6)+dc_admin)*4)/3)+dc_vrdcycle+dc_mgt,
  dc_vrdmaint       = ((dc_len*21)+(dc_dex*2*4))/4, #Rd maintenance
  
  
  #Cost of subsequent treatments
  dc_pdrd           = 1661.69/4,
  dc_pvrd           = 4288.35/4,
  
  ## HOSPICE
  c_hospice         = 8750/4,
  
  ## HOSPITALIZATION FOR AE
  cae_neut          = 13266.3,
  cae_anae          = 2567.49,
  cae_pneu          = 13760.94,
  
  ##### UTILITIES
  ## STATES per cycle (divided by weeks)
  u_pf              = 0.725/52.143,
  u_pp              = 0.65/52.143,
  
  ## AE DISUTILITY (NEGATIVE VALUES) per event
  u_neut            = (0.145*9.79)/365.25,
  u_anae            = (0.25*1.94)/365.25,
  u_pneu            = (0.14*11.14)/365.25,
  u_diar            = (0.103*4.87)/365.25,
  u_neuro           = (0.065*2.5)/365.25,
  
  #### AE RISK
  rae_drdneut       = 1-exp(-(-log(1-((136+61)/364)))/(56.2*4.345)),
  rae_drdleuk       = 1-exp(-(-log(1-((41+19+37+5)/364)))/(56.2*4.345)),
  rae_drdpneu       = 1-exp(-(-log(1-((62+5+3)/364)))/(56.2*4.345)),
  rae_drdanae       = 1-exp(-(-log(1-((1+60)/364)))/(56.2*4.345)),
  rae_drdthromb     = 1-exp(-(-log(1-((23+9)/364)))/(56.2*4.345)),
  
  rae_vrdneut       = 1-exp(-(-log(1-((29+5+1)/241)))/(84*4.345)),
  rae_vrdlymph      = 1-exp(-(-log(1-(5/241)))/(84*4.345)),
  rae_vrdanae       = 1-exp(-(-log(1-((73+41)/241)))/(84*4.345)),
  rae_vrdthromb     = 1-exp(-(-log(1-(7/241)))/(84*4.345)),
  rae_vrddiar       = 1-exp(-(-log(1-((49+3+1)/241)))/(84*4.345)),
  rae_vrdneuro      = 1-exp(-(-log(1-((76+4)/241)))/(84*4.345)),
  
  
  ##### DISCOUNT RATE (ANNUAL)
  dr_cost           = 0.03,
  dr_qaly           = 0.03)


basecase <- psm(d.params)
basecase[[3]]

