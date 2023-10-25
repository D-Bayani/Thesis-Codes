##Pathway analysis of succeeding lines## 

pathway <- service %>%
  filter(`Service Group Code`=="PRESCRIPT")

path.out <- drugadmin[, c(2)]
path.out <- unique(path.out)
path.ipt <- pmct[, c(2)]
path.ipt <- unique(path.ipt)

criteria <- rbind(path.out, path.ipt)
rm(path.out, path.ipt)

pathway <- semi_join(pathway, criteria, by = "Masked Case No")

path <- pathway[, c(2,4)]
path <- unique(path)

D.all <- full_join(daraopt, dara, by = "Masked Case No")
D.all <- D.all[, c(2,4)]
D.all <- unique(D.all) %>% mutate(Tx = "D") 

R.all <- full_join(lenopt, len, by = "Masked Case No")
R.all <- R.all[, c(2,4)]
R.all <- unique(R.all) %>% mutate(Tx = "R") 

Th.all <- full_join(thalopt, thal, by = "Masked Case No")
Th.all <- Th.all[, c(2,4)]
Th.all <- unique(Th.all) %>% mutate(Tx = "Th") 

V.all <- full_join(bortopt, bort, by = "Masked Case No")
V.all <- V.all[, c(2,4)]
V.all <- unique(V.all) %>% mutate(Tx = "V") 

K.all <- full_join(carfopt, carf, by = "Masked Case No")
K.all <- K.all[, c(2,4)]
K.all <- unique(K.all) %>% mutate(Tx = "K") 

I.all <- full_join(ixaopt, ixa, by = "Masked Case No")
I.all <- I.all[, c(2,4)]
I.all <- unique(I.all) %>% mutate(Tx = "I") 

P.all <- full_join(pomaopt, poma, by = "Masked Case No")
P.all <- P.all[, c(2,4)]
P.all <- unique(P.all) %>% mutate(Tx = "P") 

E.all <- elopt
E.all <- E.all[, c(2,4)]
E.all <- unique(E.all) %>% mutate(Tx = "E") 

C.all <- full_join(cyclopt, cyclo, by = "Masked Case No")
C.all <- C.all[, c(2,4)]
C.all <- unique(C.all) %>% mutate(Tx = "C") 

d.all <- full_join(dexopt, dexa, by = "Masked Case No")
d.all <- d.all[, c(2,4)]
d.all <- unique(d.all) %>% mutate(Tx = "d") 


#set up regimen summary
path <- path %>%
  full_join(D.all, by=c("Masked Case No", "Visit No")) %>%
  full_join(R.all, by=c("Masked Case No", "Visit No")) %>%
  full_join(Th.all, by=c("Masked Case No", "Visit No")) %>%
  full_join(V.all, by=c("Masked Case No", "Visit No")) %>%
  full_join(K.all, by=c("Masked Case No", "Visit No")) %>%
  full_join(I.all, by=c("Masked Case No", "Visit No")) %>%
  full_join(P.all, by=c("Masked Case No", "Visit No")) %>%
  full_join(E.all, by=c("Masked Case No", "Visit No")) %>%
  full_join(C.all, by=c("Masked Case No", "Visit No")) %>%
  full_join(d.all, by=c("Masked Case No", "Visit No"))

path <- path[apply(!is.na(path[,3:12]), 1, any),] #remove visits without specific drugs administered
path$`Visit No` <- as.character(path$`Visit No`)
path[is.na(path$`Visit No`), 2] <- "-"
path$`Visit No` <- as.factor(path$`Visit No`)
path[is.na(path)] <- ""
path$regimen <- gsub(" ", "", paste(path$Tx.x, path$Tx.y, path$Tx.x.x, path$Tx.y.y, path$Tx.x.x.x, path$Tx.y.y.y, path$Tx.x.x.x.x, path$Tx.y.y.y.y, path$Tx.x.x.x.x.x, path$Tx.y.y.y.y.y))
path <- path[, -c(3:12)] 
path <- path %>% mutate_if(is.character, as.factor) #recode regimens

path <- path %>%
  full_join(outcpv, by=c("Masked Case No", "Visit No")) %>%
  full_join(inptcpc, by=c("Masked Case No"))
path <- path[!is.na(path$regimen), ]

### ADD VISIT DATE FOR OUTPATIENT
visitdate <- service[, c(2,4,14)]
visitdate <- visitdate %>%
  group_by(`Masked Case No`, `Visit No`) %>%
  slice(which.min(`Service End Date`))

path <- path %>%
  full_join(visitdate, by=c("Masked Case No", "Visit No"))
path <- path[!is.na(path$regimen), ]

pathway <- pathway[, c(1,2,4)]

lot <- path[, c(1,2,3,15)]

txlines <- lot %>%
  full_join(pathway, by=c("Masked Case No", "Visit No"))
txlines <- unique(txlines)

txlines <- txlines[!is.na(txlines$regimen), ]
txlines <- txlines[!is.na(txlines$`Subject Code`), ]

txlines <- txlines[, c(5,1,2,3,4)]
txlines <- arrange(txlines, txlines$`Service End Date.y`, group_by = "Subject Code")

txlines <- txlines[order(txlines$`Subject Code`),]
txlines <- txlines %>% mutate(month = format(as.Date(txlines$`Service End Date.y`), "%Y-%m"))

subjcodecaseno <- service[, c(1,2,4)]
subjcodecaseno <- unique(subjcodecaseno)

##identify each patient exposed to drug
D.all <- D.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
D.all <- D.all %>% filter(Tx =="D")
D.all$Dara <- 1
D.all <- D.all[, c(4,6)]
D.all <- unique(D.all)

V.all <- V.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
V.all <- V.all %>% filter(Tx =="V")
V.all$V <- 1
V.all <- V.all[, c(4,6)]
V.all <- unique(V.all)

R.all <- R.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
R.all <- R.all %>% filter(Tx =="R")
R.all$R <- 1
R.all <- R.all[, c(4,6)]
R.all <- unique(R.all)

Th.all <- Th.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
Th.all <- Th.all %>% filter(Tx =="Th")
Th.all$Th <- 1
Th.all <- Th.all[, c(4,6)]
Th.all <- unique(Th.all)

K.all <- K.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
K.all <- K.all %>% filter(Tx =="K")
K.all$K <- 1
K.all <- K.all[, c(4,6)]
K.all <- unique(K.all)

I.all <- I.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
I.all <- I.all %>% filter(Tx =="I")
I.all$I <- 1
I.all <- I.all[, c(4,6)]
I.all <- unique(I.all)

P.all <- P.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
P.all <- P.all %>% filter(Tx =="P")
P.all$P <- 1
P.all <- P.all[, c(4,6)]
P.all <- unique(P.all)

E.all <- E.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
E.all <- E.all %>% filter(Tx =="E")
E.all$E <- 1
E.all <- E.all[, c(4,6)]
E.all <- unique(E.all)

C.all <- C.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
C.all <- C.all %>% filter(Tx =="C")
C.all$C <- 1
C.all <- C.all[, c(4,6)]
C.all <- unique(C.all)

d.all <- d.all %>% full_join(subjcodecaseno, by=c("Masked Case No"))
d.all <- d.all %>% filter(Tx =="d")
d.all$d <- 1
d.all <- d.all[, c(4,6)]
d.all <- unique(d.all)

drugs <- subjcodecaseno
drugs <- drugs[, c(1)]
drugs <- unique(drugs)

drugs <- drugs %>%
  full_join(D.all, by=c("Subject Code")) %>%
  full_join(V.all, by=c("Subject Code")) %>%
  full_join(R.all, by=c("Subject Code")) %>%
  full_join(Th.all, by=c("Subject Code")) %>%
  full_join(K.all, by=c("Subject Code")) %>%
  full_join(I.all, by=c("Subject Code")) %>%
  full_join(P.all, by=c("Subject Code")) %>%
  full_join(E.all, by=c("Subject Code")) %>%
  full_join(C.all, by=c("Subject Code")) %>%
  full_join(d.all, by=c("Subject Code"))
drugs[is.na(drugs)] <- 0


txlines <- txlines %>%
  full_join(drugs, by=c("Subject Code"))
txlines <- txlines[!is.na(txlines$regimen), ]

txlinescost <- txlines %>%
  full_join(outcpv, by=c("Masked Case No", "Visit No")) %>%
  full_join(inptcpc, by=c("Masked Case No"))
txlinescost <- txlinescost[!is.na(txlinescost$regimen), ]
txlinescost <- txlinescost[, -c(23,26:27,30:33)]
txlinescost[is.na(txlinescost)] <- 0

txlinescost.sum <- txlinescost
txlinescost.sum$TC = txlinescost.sum$`Total Amount.x` + txlinescost.sum$`Total Amount.y`
txlinescost.sum$TS = txlinescost.sum$`Total Subsidy.x` + txlinescost.sum$`Total Subsidy.y`
txlinescost.sum <- txlinescost.sum %>%
  mutate_if(is.character, as.factor)

##summary per drug/class
darabased <- txlinescost.sum %>%
  filter(reg=="D" | reg=="Dd"| reg=="DKd"| reg=="DPd"| reg=="DRd"| reg=="DVCd"| reg=="DVd") %>%
  group_by(`Subject Code.x`, month, reg) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)

vbased <- txlinescost.sum %>%
  filter(reg=="V" | reg=="VCd"| reg=="Vd"| reg=="VR"| reg=="VRd"| reg=="VthCd"| reg=="VThd") %>%
  group_by(`Subject Code.x`, month, reg) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)


DRd <- txlinescost.sum %>%
  filter(reg=="DRd") %>%
  group_by(`Subject Code.x`, month) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)

VRd <- txlinescost.sum %>%
  filter(reg=="VRd") %>%
  group_by(`Subject Code.x`, month) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)

LOT1 <- txlinescost.sum %>%
  filter(LOT1=="1", exclude!="1") %>%
  group_by(`Subject Code.x`, month, reg) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)

LOT2 <- txlinescost.sum %>%
  filter(LOT2=="1", exclude!="1") %>%
  group_by(`Subject Code.x`, month, reg) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)

LOT3 <- txlinescost.sum %>%
  filter(LOT3=="1", exclude!="1") %>%
  group_by(`Subject Code.x`, month, reg) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)

LOT4 <- txlinescost.sum %>%
  filter(LOT4=="1", exclude!="1") %>%
  group_by(`Subject Code.x`, month, reg) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)

#get succeeding lines
##DRd
postDRd <- txlinescost.sum %>%
  mutate(PostDRd = ifelse(LOT1=="1" & reg=="DRd", 1, 0))
postDRd <- postDRd[, c(1,29)]
postDRd <- unique(postDRd)
postDRd <- postDRd %>% filter(PostDRd ==1)

postDRdcost <- txlinescost.sum %>%
  semi_join(postDRd, by = "Subject Code.x")

postDRdcost <- postDRdcost %>%
  filter(LOT2=="1" | LOT3=="1" | LOT4=="1") %>%
  group_by(`Subject Code.x`, month, reg) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)

##VRd
postVRd <- txlinescost.sum %>%
  mutate(PostVRd = ifelse(LOT1=="1" & reg=="VRd", 1, 0))
postVRd <- postVRd[, c(1,29)]
postVRd <- unique(postVRd)
postVRd <- postVRd %>% filter(PostVRd ==1)

postVRdcost <- txlinescost.sum %>%
  semi_join(postVRd, by = "Subject Code.x")

postVRdcost <- postVRdcost %>%
  filter(LOT2=="1" | LOT3=="1" | LOT4=="1") %>%
  group_by(`Subject Code.x`, month, reg) %>%
  summarise(TotalCost = sum(TC, na.rm=T), TotalSubsidy = sum(TS, na.rm=T)) %>%
  mutate(PS = TotalSubsidy/TotalCost)


summary(darabased)
summary(vbased)
summary(DRd)
summary(VRd)
summary(LOT1)
summary(LOT2)
summary(LOT3)
summary(LOT4)
summary(postDRdcost)
summary(postVRdcost)

##get number of txlines
no.txlines <- txlinescost.sum[, c(1,7:11)]
no.txlines <- unique(no.txlines)
no.txlines <- no.txlines %>%
  filter(exclude!= "1") %>%
  group_by(`Subject Code.x`) %>%
  summarise(across(c(LOT1, LOT2, LOT3, LOT4), sum))
no.txlines <- no.txlines %>%
  mutate(no.txlines = LOT1 + LOT2 + LOT3 + LOT4)
summary(no.txlines)

#new summary per regimen
meancostreg <- txlinescost.sum[, c(1:6, 12, 23:28)]
meancostreg[meancostreg == 0] <- NA
meancostreg <- meancostreg %>%
  group_by(`Subject Code.x`, month, reg) %>%
  summarise(mctotal = mean(TC), mcoutpt = mean(`Total Amount.x`, na.rm=TRUE), mcinpt = mean(`Total Amount.y`, na.rm=TRUE))

meancostreg2 <- meancostreg %>%
  group_by(reg) %>%
  summarise(mctotal = mean(mctotal), mcoutpt = mean(mcoutpt, na.rm=TRUE), mcinpt = mean(mcinpt, na.rm=TRUE))

write.csv(meancostreg2, "meancost.reg.csv")

p.drd <- meancostreg %>%
  filter(reg== "VCd" | reg== "KCd" | reg== "VPd" | reg== "Kd" | reg== "PCd") 
p.drd <- as.data.frame(describe(p.drd[, c(3:6)], IQR=TRUE))

p.vrd <- meancostreg %>%
  filter(reg== "DKd" | reg== "Kd" | reg== "KCd" | reg== "DPd" | reg== "PCd" | reg== "Pd") 
p.vrd <- as.data.frame(describe(p.vrd[, c(3:6)], IQR=TRUE))

dara.based <- meancostreg %>%
  filter(reg == "Dd" | reg == "DKd" | reg == "DPd"| reg == "DRd" | reg == "DVCd" | reg == "DVd" | reg == "DVRd"| reg == "DVThd")
dara.based <- as.data.frame(describe(dara.based[, c(3:6)], IQR=TRUE))

bort.based <- meancostreg %>%
  filter(reg == "V" | reg == "VCd" | reg == "Vd"| reg == "VRd" | reg == "DVCd" | reg == "DVd" | reg == "DVRd"| reg == "DVThd"| reg == "VThCd")
bort.based <- as.data.frame(describe(bort.based[, c(3:6)], IQR=TRUE))

darabased.cost <- as.data.frame(describe(darabased[, c(3:6)], IQR=TRUE))
vbased.cost <- as.data.frame(describe(vbased[, c(3:6)], IQR=TRUE))
lot1.cost <- as.data.frame(describe(LOT1[, c(3:6)], IQR=TRUE))
lot2.cost <- as.data.frame(describe(LOT2[, c(3:6)], IQR=TRUE))
lot3.cost <- as.data.frame(describe(LOT3[, c(3:6)], IQR=TRUE))
lot4.cost <- as.data.frame(describe(LOT4[, c(3:6)], IQR=TRUE))