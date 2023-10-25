################################################################################
## DEMOGRAPHICS ##
################################################################################

#check SG and PRs
billing <- billing[, c(1, 8)]
billing <- billing %>%
  filter(`Billing Class Category`!="NA")
billing <- unique(billing)
summary(billing)

demog <- full_join(demog, billing, by = c("Subject Code"))

##exclude SG and PR##

sgpr <- billing %>%
  filter(`Billing Class Category`== "SINGAPOREAN" | `Billing Class Category`== "PR")

#exclude individuals from datasets
demog <- semi_join(demog, sgpr, by = "Subject Code")
service <- semi_join(service, sgpr, by = "Subject Code")
casemove <- semi_join(casemove, sgpr, by = "Subject Code")
diagnosis <- semi_join(diagnosis, sgpr, by = "Subject Code")
ipharm <- semi_join(ipharm, sgpr, by = "Subject Code")

summary(demog)

################################################################################
## GET SURVIVAL ##
################################################################################

survival <- demog[, c(1,3)]
survival$event = 1
survival <- survival %>%
  mutate(event = ifelse(is.na(`Death Date`), 0,1))

#get earliest date
survstart <- casemove %>%
  group_by(`Subject Code`) %>%
  slice(which.min(`Movement Start Date`))

survstart <- survstart[, c(1,23)] #select cols

survival <- survival %>%
  full_join(survstart, by=c("Subject Code"))

survival$lastdate <- survival$`Death Date`
survival$lastdate[is.na(survival$lastdate)] <- "2022-04-30"
survival$time = as.numeric(difftime(survival$lastdate, survival$`Movement Start Date`, units = "days"))

survival <- survival[, c(1,3,6)]
km <- survfit(Surv(time, event)~1, data=survival, type = "kaplan-meier")

################################################################################
## SEPARATE INPATIENT WITH OUTPATIENT ##
################################################################################

inpatient <- casemove %>%
  filter(`Case Type`=="1")

inpatient <- unique(inpatient[, 3])

outpatient <- casemove %>%
  filter(`Case Type`=="2")

outpatient <- outpatient[, c(3,5)]


################################################################################
## COST PER CASE (ALL) ##
################################################################################

#set up for total cost per case
totalcostpercase <- service %>%
  select(c(1,2,21,23,25))

#get total amount per case
totalamount <- totalcostpercase %>%
  group_by(`Masked Case No`) %>%
  summarise(`Total Amount` = sum(`Gross Amount`, `Tax Amount`, na.rm=T))

#get total subsidy per case
totalsubsidy <- totalcostpercase %>%
  group_by(`Masked Case No`) %>%
  summarise(`Total Subsidy` = sum(`Subsidy Amount`, na.rm=T))

#add total amount as a variable
totalcostpercase <- totalcostpercase %>%
  full_join(totalamount, by=c("Masked Case No")) %>%
  full_join(totalsubsidy, by=c("Masked Case No"))

#omit other variables
totalcostpercase <- totalcostpercase %>% 
  group_by(`Masked Case No`) %>%
  slice(which.max(`Total Amount`)) %>%
  slice(which.max(`Total Subsidy`)) %>%
  select(c(`Subject Code`,`Masked Case No`,`Total Amount`, `Total Subsidy`))

totalcostpercase$percentsubsidy = (totalcostpercase$`Total Subsidy`/totalcostpercase$`Total Amount`) #new var for % subsidy


#mix of inpatient and outpatient cases!

################################################################################
## COST PER PATIENT (ALL) ##
################################################################################

#set up for total cost per patient
totalcostperpx <- service %>%
  select(c(1,2,21,23,25))

#get total amount per case
totalamountpx <- totalcostperpx %>%
  group_by(`Subject Code`) %>%
  summarise(`Total Amount` = sum(`Gross Amount`, `Tax Amount`, na.rm=T))

#get total subsidy per case
totalsubsidypx <- totalcostperpx %>%
  group_by(`Subject Code`) %>%
  summarise(`Total Subsidy` = sum(`Subsidy Amount`, na.rm=T))

#add total amount as a variable
totalcostperpx <- totalcostperpx %>%
  full_join(totalamountpx, by=c("Subject Code")) %>%
  full_join(totalsubsidypx, by=c("Subject Code"))

#omit other variables
totalcostperpx <- totalcostperpx %>% 
  group_by(`Subject Code`) %>%
  slice(which.max(`Total Amount`)) %>%
  slice(which.max(`Total Subsidy`)) %>%
  select(c(`Subject Code`,`Masked Case No`,`Total Amount`, `Total Subsidy`))

totalcostperpx$percentsubsidy = (totalcostperpx$`Total Subsidy`/totalcostperpx$`Total Amount`) #new var for % subsidy


#mix of inpatient and outpatient!


################################################################################
## INPATIENT LOS AND COST PER CASE ##
################################################################################

#get service data of inpatient cases only
inptsvc<- semi_join(service, inpatient, by=c("Masked Case No")) 

#get admission date
start <- inptsvc %>%
  group_by(`Masked Case No`) %>%
  slice(which.min(`Service Start Date`))
start <- start[, c(2,13)] #select cols

#get discharge date
end <- inptsvc %>%
  group_by(`Masked Case No`) %>%
  slice(which.max(`Service End Date`))
end <- end[, c(2,14)] #select cols

#combine and calculate duration of stay
los <- start %>%
  full_join(end, by=c("Masked Case No"))

los$dur = as.numeric(difftime(los$`Service End Date`, los$`Service Start Date`, units = "days"))

#get cost per inpatient case
inptcpc <- semi_join(totalcostpercase, inpatient, by=c("Masked Case No")) %>%
  full_join(los, by=c("Masked Case No"))

rm(start, end) #remove dfs


################################################################################
## IPHARM COST PER CASE ## INPATIENT
################################################################################
##inpatient ipharm
ipharminpt <- semi_join(ipharm, inpatient, by=c("Masked Case No")) 

#get total amount per case
phinptamt <- ipharminpt %>%
  group_by(`Masked Case No`) %>%
  summarise(`Total Amount` = sum(`Gross Amount`, na.rm=T))

#get total subsidy per case
phinptsub <- ipharminpt %>%
  group_by(`Masked Case No`) %>%
  summarise(`Total Subsidy` = sum(`Subsidy Amount`, na.rm=T))

#add total amount as a variable
ipharminptcpc <- ipharminpt %>%
  full_join(phinptamt, by=c("Masked Case No")) %>%
  full_join(phinptsub, by=c("Masked Case No"))

#omit other variables
ipharminptcpc <- ipharminptcpc %>% 
  group_by(`Masked Case No`) %>%
  slice(which.max(`Total Amount`)) %>%
  slice(which.max(`Total Subsidy`)) %>%
  select(c(`Subject Code`,`Masked Case No`,`Total Amount`, `Total Subsidy`))

ipharminptcpc$percentsubsidy = (ipharminptcpc$`Total Subsidy`/ipharminptcpc$`Total Amount`) #new var for % subsidy

summary(ipharminptcpc)


################################################################################
## IPHARM COST PER CASE ## OUTPATIENT
################################################################################
##outpatient ipharm
ipharmopt <- semi_join(ipharm, outpatient, by=c("Masked Case No", "Visit No")) 

#get total amount per case
phoptamt <- ipharmopt %>%
  group_by(`Masked Case No`, `Visit No`) %>%
  summarise(`Total Amount` = sum(`Gross Amount`, na.rm=T))

#get total subsidy per case
phoptsub <- ipharmopt %>%
  group_by(`Masked Case No`, `Visit No`) %>%
  summarise(`Total Subsidy` = sum(`Subsidy Amount`, na.rm=T))

#add total amount as a variable
ipharmoptcpv <- ipharmopt %>%
  full_join(phoptamt, by=c("Masked Case No", "Visit No")) %>%
  full_join(phoptsub, by=c("Masked Case No", "Visit No"))

#omit other variables
ipharmoptcpv <- ipharmoptcpv %>% 
  group_by(`Masked Case No`, `Visit No`) %>%
  slice(which.max(`Total Amount`)) %>%
  slice(which.max(`Total Subsidy`)) %>%
  select(c(`Subject Code`,`Masked Case No`,`Visit No`,`Total Amount`, `Total Subsidy`))

ipharmoptcpv$percentsubsidy = (ipharmoptcpv$`Total Subsidy`/ipharmoptcpv$`Total Amount`) #new var for % subsidy

summary(ipharmoptcpv)

################################################################################
## DIAGNOSIS BASED COST PER CASE FOR INPATIENT ##
################################################################################

#diagnosis based inpatient data; combine inpatient cases with diagnosis (primary and secondary)
inptdiag <- semi_join(diagnosis, inpatient, by=c("Masked Case No")) 

#get primary diagnoses for inpatient
inptdiagpri <- inptdiag %>%
  filter(`Diagnosis Type`=="P", `Diagnosis Category`=="D2")

inptdiagpri <- full_join(inptdiagpri, inptcpc, by=c("Masked Case No"))
inptdiag <- full_join(inptdiag, inptcpc, by=c("Masked Case No")) 

topinptdiag <- as.data.frame(summary(inptdiagpri$`Diagnosis Description`)) #Get most common causes of hospitalization
summary(inptdiagpri$`Total Amount`) #get summary costs of hospitalization

######NEUTROPAENIA##### 
##get only those with primary diagnosis and main treatment(P/D2)
neut <- inptdiag %>%
  filter(`Diagnosis Code`=="D70",`Diagnosis Type`=="P", `Diagnosis Category`=="D2")
summary(neut)


######PNEUMONIA##### 
##get only those with primary diagnosis and main treatment
pneu <- inptdiag %>%
  filter(`Diagnosis Code`=="J189",`Diagnosis Type`=="P", `Diagnosis Category`=="D2")
summary(pneu)


######ANAEMIA##### 
##get all those with diagnosis and main treatment (D1/D2)
anae <- inptdiag %>%
  filter(`Diagnosis Code`=="D630"| `Diagnosis Code`=="D649",`Diagnosis Type`=="P",`Diagnosis Category`=="D2")
summary(anae)


######LYMPHOPAENIA/LEUKOPAENIA#####  FOR VALIDATION W CLINICIANS
#get all those with diagnosis 
lymph <- inptdiag %>%
  filter(`Diagnosis Code`=="D728")
summary(lymph)

######THROMBOCYTOPAENIA#####  
##get only those with primary diagnosis and main treatment
thromb <- inptdiag %>%
  filter(`Diagnosis Code`=="D696" | `Diagnosis Code`=="D695",`Diagnosis Type`=="P", `Diagnosis Category`=="D2")
summary(thromb)

##get all those with diagnosis and main treatment (D1/D2)
thromb <- inptdiag %>%
  filter(`Diagnosis Code`=="D696" | `Diagnosis Code`=="D695",`Diagnosis Category`=="D2")
summary(thromb)


######Gastroenteritis#####
diar <- inptdiag %>%
  filter(`Diagnosis Code`=="A099" | `Diagnosis Code`=="K521",`Diagnosis Type`=="P", `Diagnosis Category`=="D2")
summary(diar)

######neuropathy#####
neur <- inptdiag %>%
  filter(`Diagnosis Code`=="G629" | `Diagnosis Code`=="G620",`Diagnosis Type`=="P")
summary(neur)

######OTHER COMPLICATIONS#####
#get all those with diagnosis 
othcomp <- inptdiag %>%
  filter(`Diagnosis Code`=="T808" | `Diagnosis Code`=="T801")
summary(othcomp)


######PALLIATIVE CARE#####
#get all those with diagnosis 
palliative <- inptdiag %>%
  filter(`Diagnosis Code`=="Z515")
summary(palliative)


################################################################################
## PHARMACOTHERAPY INPATIENT ##
################################################################################
##get all those with diagnosis and main treatment (D1/D2)
pmct <- inptdiag %>%
  filter(`Diagnosis Code`=="Z511",`Diagnosis Category`=="D2")
summary(pmct)


################################################################################
##COMBINE WITH IPHARM AND SERVICE DATA
################################################################################

pmctdetail <- pmct %>%
  full_join(ipharminpt, by=c("Masked Case No"))

pmctdetail <- pmctdetail[!is.na(pmctdetail$`Subject Code.x`), ] #remove rows with NA, those not for pmct cases
pmctdetail <- pmctdetail[, -c(3:6, 11, 18:20)]


pmctsvc <- pmct %>%
  full_join(inptsvc, by=c("Masked Case No"))

pmctsvc <- pmctsvc[!is.na(pmctsvc$`Subject Code.x`), ] #remove rows with NA, those not for pmct cases


###GET COSTS PER REGIMEN###
#get all dara patients
dara <- pmctdetail %>%
  filter(`Item Code`=="0010-10-180-W" | `Item Code`=="0010-10-179-E" | `Item Code`=="0310-10-017-T" | `Item Code`=="0310-10-018-T")
dara <- semi_join(pmctdetail,dara, by=c("Masked Case No")) #get full details of those on dara

#get paid dara patients
darapaid <- pmctdetail %>%
  filter(`Item Code`=="0010-10-180-W" | `Item Code`=="0010-10-179-E")
darapaid <- semi_join(pmctdetail,darapaid, by=c("Masked Case No")) #get full details of those on paid dara

#get len patients
len <- pmctdetail %>%
  filter(`Item Code`=="0004-92-128-I" | `Item Code`=="0004-92-127-J"| `Item Code`=="0304-92-007-T")
len <- semi_join(pmctdetail,len, by=c("Masked Case No")) #get full details of those on len

#get thal patients
thal <- pmctdetail %>%
  filter(`Item Code`=="0004-92-062-I")
thal <- semi_join(pmctdetail,thal, by=c("Masked Case No")) #get full details of those on thal

#get bort patients
bort <- pmctdetail %>%
  filter(`Item Code`=="0010-10-073-E")
bort <- semi_join(pmctdetail,bort, by=c("Masked Case No")) #get full details of those on bort

#get carf patients
carf <- pmctdetail %>%
  filter(`Item Code`=="0010-92-164-L")
carf <- semi_join(pmctdetail,carf, by=c("Masked Case No")) #get full details of those on carf

#get cyclo patients
cyclo <- pmctdetail %>%
  filter(`Item Code`=="0004-10-005-W")
cyclo <- semi_join(pmctdetail,cyclo, by=c("Masked Case No")) #get full details of those on cyclo

#get dexa patients
dexa <- pmctdetail %>%
  filter(`Item Code`=="0004-68-083-K" | `Item Code`=="0010-68-006-E")
dexa <- semi_join(pmctdetail,dexa, by=c("Masked Case No")) #get full details of those on dexa

#get ixa patients
ixa <- pmctdetail %>%
  filter(`Item Code`=="2125-54-96X-3" | `Item Code`=="4792-94-81X-0")
ixa <- semi_join(pmctdetail,ixa, by=c("Masked Case No")) #get full details of those on ixa

#get poma patients
poma <- pmctdetail %>%
  filter(`Item Code`=="1843-35-43X-1" | `Item Code`=="0004-92-172-J")
poma <- semi_join(pmctdetail,poma, by=c("Masked Case No")) #get full details of those on poma

#get rituximab patients
rit <- pmctdetail %>%
  filter(`Item Code`=="0010-10-053-K")
rit <- semi_join(pmctdetail,rit, by=c("Masked Case No")) #get full details of those on ritux


################################################################################
## OUTPATIENT COST PER CASE##
################################################################################

#get service data of outpatient cases only
outsvc<- semi_join(service, outpatient, by=c("Masked Case No", "Visit No")) 

#get cost per visit per case
totalapv <- outsvc %>%
  group_by(`Masked Case No`, `Visit No`) %>%
  summarize(`Total Amount` = sum(`Gross Amount`, `Tax Amount`, na.rm=T))

#get total subsidy per case
totalspv <- outsvc %>%
  group_by(`Masked Case No`, `Visit No`) %>%
  summarise(`Total Subsidy` = sum(`Subsidy Amount`, na.rm=T))

#add total amount as a variable
totalcpv <- outsvc %>%
  full_join(totalapv, by=c("Masked Case No", "Visit No")) %>%
  full_join(totalspv, by=c("Masked Case No", "Visit No"))

#omit other variables
outcpv <- totalcpv %>% 
  group_by(`Masked Case No`, `Visit No`) %>%
  slice(which.max(`Total Amount`)) %>%
  slice(which.max(`Total Subsidy`)) %>%
  select(c(`Subject Code`,`Masked Case No`,`Visit No`,`Total Amount`, `Total Subsidy`))

outcpv$percentsubsidy = (outcpv$`Total Subsidy`/outcpv$`Total Amount`) #new var for % subsidy

summary(outcpv)

################################################################################
## OUTPATIENT COST - DIAGNOSIS BASED ##
################################################################################
#diagnosis based inpatient data; combine inpatient cases with diagnosis (primary and secondary)
outdiag <- semi_join(diagnosis, outpatient, by=c("Masked Case No", "Visit No")) 

#get primary diagnoses for outpatient
outdiagpri <- outdiag %>%
  filter(`Diagnosis Type`=="P")

outdiagpri <- full_join(outdiagpri, outcpv, by=c("Masked Case No"))
outdiag <- full_join(outdiag, outcpv, by=c("Masked Case No")) 

topoutptdiag <- as.data.frame(summary(outdiagpri$`Diagnosis Description`)) #Get most common causes of outpatient visit
summary(outdiagpri$`Total Amount`) #get summary costs of outpatient visits


######OUTPATIENT DRUG ADMINISTRATION##### what is the best way to get outpatient drug admin?
##get only those with primary diagnosis
drugadmin <- outdiag %>%
  filter(`Diagnosis Code`=="C90" | `Diagnosis Code`=="C900" | `Diagnosis Code`=="C9000")
summary(drugadmin)

##get cost based on chemo administration procedure
chemo <- outsvc %>%
  filter(`Service Code`=="CM0013" | `Service Code`=="DS0016" | `Service Code`=="DS0015" | `Service Code`=="DS0014")
chemocost <- semi_join(outcpv, chemo, by=c("Masked Case No", "Visit No"))
summary(chemocost)


################################################################################
##COMBINE WITH IPHARM AND SERVICE DATA
################################################################################

###GET COSTS PER REGIMEN###
#get all dara patients
daraopt <- ipharmopt %>%
  filter(`Item Code`=="0010-10-180-W" | `Item Code`=="0010-10-179-E" | `Item Code`=="0310-10-017-T" | `Item Code`=="0310-10-018-T")

#get paid dara patients
darapaidopt <- ipharmopt %>%
  filter(`Item Code`=="0010-10-180-W" | `Item Code`=="0010-10-179-E")

#get len patients
lenopt <- ipharmopt %>%
  filter(`Item Code`=="0004-92-128-I" | `Item Code`=="0004-92-127-J"| `Item Code`=="0304-92-007-T")

#get thal patients
thalopt <- ipharmopt %>%
  filter(`Item Code`=="0004-92-062-I")

#get bort patients
bortopt <- ipharmopt %>%
  filter(`Item Code`=="0010-10-073-E"	| `Item Code`=="0375-16-053-R")

#get carf patients
carfopt <- ipharmopt %>%
  filter(`Item Code`=="0010-92-164-L")

#get ixa patients
ixaopt <- ipharmopt %>%
  filter(`Item Code`=="2125-54-96X-3" | `Item Code`=="4792-94-81X-0")

#get poma patients
pomaopt <- ipharmopt %>%
  filter(`Item Code`=="1843-35-43X-1" | `Item Code`=="0004-92-172-J")

#get elo patients
elopt <- ipharmopt %>%
  filter(`Item Code`=="0392-90-222-T")

#get cyclo patients
cyclopt <- ipharmopt %>%
  filter(`Item Code`=="0004-10-005-W" | `Item Code`=="0010-10-039-E")

#get dex patients
dexopt <- ipharmopt %>%
  filter(`Item Code`=="0004-68-083-K" | `Item Code`=="0004-68-007-E" | `Item Code`=="0010-68-006-E")


###get mean unit price of intervention and comparator
daracost <- service %>%
  filter(`Service Code` =="001010180W")
summary(daracost$`Service Unit Price`)

bortcost <- service %>%
  filter(`Service Code` =="001010073E")
summary(bortcost$`Total Standard Fixed Cost`)

lencost <- service %>%
  filter(`Service Code` =="000492128I")
summary(lencost$`Service Unit Price`) 

dexcost <- service %>%
  filter(`Service Code` =="000468083K")
summary(dexcost$`Service Unit Price`) # for 4mg, x 10 to get 40 mg as per MAIA regimen


###get mean costs of chemotherapy administration
chemotxmonitoring <- outsvc %>%
  filter(`Service Code` =="CM0034")
summary(chemotxmonitoring$`Service Unit Price`)

prepinfusionbag <- outsvc %>%
  filter(`Service Code` =="037516009R")
summary(prepinfusionbag$`Service Unit Price`)

prepinfusionadapter <- outsvc %>%
  filter(`Service Code` =="037516058R")
summary(prepinfusionadapter$`Service Unit Price`)

prepsyringes <- outsvc %>%
  filter(`Service Code` =="037516010R")
summary(prepsyringes$`Service Unit Price`)

ivinsertion <- outsvc %>%
  filter(`Service Code` =="CM0012")
summary(ivinsertion$`Service Unit Price`)

chemosetup <- outsvc %>%
  filter(`Service Code` =="CM0037")
summary(chemosetup$`Service Unit Price`)

ivfluid <- outsvc %>%
  filter(`Service Code` =="DS0030")
summary(ivfluid$`Service Unit Price`)


###get disease monitoring unit costs
crp <- outsvc %>%
  filter(`Service Code` =="BI1171")
summary(crp$`Service Unit Price`)

fbc <- outsvc %>%
  filter(`Service Code` =="HE1118")
summary(fbc$`Service Unit Price`)

renal <- outsvc %>%
  filter(`Service Code` =="BI1813")
summary(renal$`Service Unit Price`)

urea <- outsvc %>%
  filter(`Service Code` =="BI1939")
summary(urea$`Service Unit Price`)

sflc <- outsvc %>%
  filter(`Service Code` =="BI1568")
summary(sflc$`Service Unit Price`)

serummprotein <- outsvc %>%
  filter(`Service Code` =="BI1842")
summary(serummprotein$`Service Unit Price`)


################################################################################
## GET REGIMEN PER PATIENT
################################################################################
### FOR EASE OF SUMMARY ONLY --> CHANGE TO PAID DARA
D <- semi_join(totalcpv, daraopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "D") 
D <- D[, c(2,4,41)]
D <- unique(D)

PaidD <- semi_join(totalcpv, darapaidopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "PaidD")
PaidD <- PaidD[, c(2,4,41)]
PaidD <- unique(PaidD)

R <- semi_join(totalcpv, lenopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "R")
R <- R[, c(2,4,41)]
R <- unique(R)

Th <- semi_join(totalcpv, thalopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "Th")
Th <- Th[, c(2,4,41)]
Th <- unique(Th)

V <- semi_join(totalcpv, bortopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "V")
V <- V[, c(2,4,41)]
V <- unique(V)

K <- semi_join(totalcpv, carfopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "K")
K <- K[, c(2,4,41)]
K <- unique(K)

I <- semi_join(totalcpv, ixaopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "I")
I <- I[, c(2,4,41)]
I <- unique(I)

P <- semi_join(totalcpv, pomaopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "P")
P <- P[, c(2,4,41)]
P <- unique(P)

E <- semi_join(totalcpv, elopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "E")
E <- E[, c(2,4,41)]
E <- unique(E)

C <- semi_join(totalcpv, cyclopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "C")
C <- C[, c(2,4,41)]
C <- unique(C)

d <- semi_join(totalcpv, dexopt, by=c("Masked Case No", "Visit No")) %>%
  group_by(`Masked Case No`, `Visit No`) %>% mutate(Tx = "d")
d <- d[, c(2,4,41)]
d <- unique(d)


#set up regimen summary
regimen <- outpatient %>%
  full_join(D, by=c("Masked Case No", "Visit No")) %>%
  full_join(PaidD, by=c("Masked Case No", "Visit No")) %>%
  full_join(R, by=c("Masked Case No", "Visit No")) %>%
  full_join(Th, by=c("Masked Case No", "Visit No")) %>%
  full_join(V, by=c("Masked Case No", "Visit No")) %>%
  full_join(K, by=c("Masked Case No", "Visit No")) %>%
  full_join(I, by=c("Masked Case No", "Visit No")) %>%
  full_join(P, by=c("Masked Case No", "Visit No")) %>%
  full_join(E, by=c("Masked Case No", "Visit No")) %>%
  full_join(C, by=c("Masked Case No", "Visit No")) %>%
  full_join(d, by=c("Masked Case No", "Visit No"))

regimen <- regimen[apply(!is.na(regimen[,3:13]), 1, any),] #remove visits without specific drugs administered
regimen[is.na(regimen)] <- ""
regimen$regimen <- gsub(" ", "", paste(regimen$Tx.x, regimen$Tx.x.x, regimen$Tx.y.y, regimen$Tx.x.x.x, regimen$Tx.y.y.y, regimen$Tx.x.x.x.x, regimen$Tx.y.y.y.y, regimen$Tx.x.x.x.x.x, regimen$Tx.y.y.y.y.y, regimen$Tx))
regimen <- regimen %>%
  mutate(PaidDara = ifelse(Tx.y =="PaidD", 1,0))
regimen <- regimen[, -c(3:13)] 
regimen <- regimen %>% mutate_if(is.character, as.factor) #recode regimens

regimen <- regimen %>%
  full_join(outcpv, by=c("Masked Case No", "Visit No"))
regimen <- regimen[!is.na(regimen$regimen), ]


visitdate <- outsvc[, c(2,4,14)]
visitdate <- visitdate %>%
  group_by(`Masked Case No`, `Visit No`) %>%
  slice(which.min(`Service End Date`))

regimen <- regimen %>%
  full_join(visitdate, by=c("Masked Case No", "Visit No"))
regimen <- regimen[!is.na(regimen$regimen), ]


#summarize top regimens
topregimens <- as.data.frame(summary(regimen$regimen))

#get mean cost per regimen
meancostregimen <- regimen %>%
  group_by(`regimen`) %>%
  summarise(`Mean Cost per Regimen` = mean(`Total Amount`, na.rm=T))




################################################################################
## CREATE DATA FRAME FOR SUMMARY STATS ##
################################################################################

#GET SUMMARY STATS FOR ALL COMPLICATIONS
s.neut <- as.data.frame(describe(neut[, c(12:14, 17)], IQR=TRUE)) %>% mutate(ae = "neut")
s.pneu <- as.data.frame(describe(pneu[, c(12:14, 17)], IQR=TRUE)) %>% mutate(ae = "pneu")
s.anae <- as.data.frame(describe(anae[, c(12:14, 17)], IQR=TRUE)) %>% mutate(ae = "anae")
s.lymph <- as.data.frame(describe(lymph[, c(12:14, 17)], IQR=TRUE)) %>% mutate(ae = "lymph")
s.diar <- as.data.frame(describe(diar[, c(12:14, 17)], IQR=TRUE)) %>% mutate(ae = "diar")
s.neur <- as.data.frame(describe(neur[, c(12:14, 17)], IQR=TRUE)) %>% mutate(ae = "neur")

#GET SUMMARY STATS FOR COSTS AND LOS
s.tcpc<- as.data.frame(describe(totalcostpercase, IQR=TRUE)) %>% mutate(label = "tcpc")
s.tcpp<- as.data.frame(describe(totalcostperpx, IQR=TRUE)) %>% mutate(label = "tcpp")
s.ipcpc<- as.data.frame(describe(inptcpc, IQR=TRUE)) %>% mutate(label = "ipcpc")
s.ipharminptcpc<- as.data.frame(describe(ipharminptcpc, IQR=TRUE)) %>% mutate(label = "ipharminptcpc")
s.ipharmoptcpv<- as.data.frame(describe(ipharmoptcpv, IQR=TRUE)) %>% mutate(label = "ipharmoptcpv")
s.outcpv<- as.data.frame(describe(outcpv, IQR=TRUE)) %>% mutate(label = "outcpv")
