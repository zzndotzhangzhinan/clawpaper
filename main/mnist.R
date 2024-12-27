## ------------------------------------------------------------------------------------------
library(devtools)
library(randomForest)


MNIST <- readRDS("MNIST.rds")

RawTr = MNIST$train
RawTest = MNIST$test

All_0 = rbind(RawTest[which(RawTest$y==0),],  RawTr[which(RawTr$y==0),])
signal_6 = RawTest[which(RawTest$y==6),]
signal_6 = signal_6[1:120,]
signal_9 = RawTest[which(RawTest$y==9),]
signal_9 = signal_9[1:500,]
TestDat.1 = rbind(signal_6, All_0[1:980,])
TestDat.2 = rbind(signal_9, All_0[981:(980+1500),])

remain0 = All_0[-(1:980+1500),]
CalDat.1 = remain0[1:1100,]
CalDat.2 = remain0[1101:(1100+2000),]
Tr0 = rbind(remain0[-(1:(1100+2000)),],  RawTr[intersect(which(RawTr$y==1),1:20),], 
            RawTr[intersect(which(RawTr$y==2),1:20),])

bag1 = rbind(TestDat.1,CalDat.1)
bag2 = rbind(TestDat.2,CalDat.2)

set.seed(123321)

myrf1 = randomForest(x=Tr0[,1:784], y=as.factor(as.vector(Tr0[,785])), 
                     xtest = bag1[,1:784])
claw.te1 = myrf1[["test"]][["votes"]][1:1100,1]
claw.cal1 = myrf1[["test"]][["votes"]][1101:2200,1]

myrf2 = randomForest(x=Tr0[,1:784], y=as.factor(as.vector(Tr0[,785])), 
                     xtest = bag2[,1:784])
claw.te2 = myrf2[["test"]][["votes"]][1:2000,1]
claw.cal2 = myrf2[["test"]][["votes"]][2001:4000,1]


pool.rf = randomForest(x=Tr0[,1:784], y=as.factor(as.vector(Tr0[,785])), 
                       xtest = rbind(bag1,bag2)[,1:784])
pad.te1 = pool.rf[["test"]][["votes"]][1:1100,1]
pad.cal1 = pool.rf[["test"]][["votes"]][1101:2200,1]
pad.te2 = pool.rf[["test"]][["votes"]][2201:(2200+2000),1]
pad.cal2 = pool.rf[["test"]][["votes"]][4201:(4200+2000),1]
pad.te = c(pad.te1,pad.te2)
pad.cal = c(pad.cal1,pad.cal2)


######### PROCEDURES ############
source("group_func_storey.R")
al = 0.05
theta1 = c(rep(1,120),rep(0,980))
theta2 = c(rep(1,500),rep(0,1500))
theta = c(theta1,theta2)

#new mwthod
de.new = conf.q( c(claw.te1,claw.te2), c(claw.cal1,claw.cal2), al)
fdr.new = sum(de.new*(1-theta)) / max(1, sum(de.new) )
tp.new = sum(de.new*theta) #/ max(1, sum(theta) )


#pooled Adadetect
de.pad = Adadetect(pad.te,pad.cal,al)
fdr.pad = sum(de.pad*(1-theta)) / max(1, sum(de.pad) )
tp.pad = sum(de.pad*theta) 

#separate Adadetect
de.sad1 = Adadetect(claw.te1,claw.cal1,al)
de.sad2 = Adadetect(claw.te2,claw.cal2,al)
de.sad = c(de.sad1,de.sad2)
fdr.sad = sum(de.sad*(1-theta)) /max(1, sum(de.sad) )
tp.sad = sum(de.sad*theta) 




## ------------------------------------------------------------------------------------------
All_0 = rbind(RawTest[which(RawTest$y==0),],  RawTr[which(RawTr$y==0),])
signal_6 = RawTest[which(RawTest$y==8),]
signal_6 = signal_6[1:120,]
signal_9 = RawTest[which(RawTest$y==6),]
signal_9 = signal_9[1:500,]
TestDat.1 = rbind(signal_6, All_0[1:1080,])
TestDat.2 = rbind(signal_9, All_0[1081:(1080+1500),])

remain0 = All_0[-(1:1080+1500),]
CalDat.1 = remain0[1:1200,]
CalDat.2 = remain0[1201:(1200+2000),]
Tr0 = rbind(remain0[-(1:(1200+2000)),],  RawTr[intersect(which(RawTr$y==1),1:20),], 
            RawTr[intersect(which(RawTr$y==2),1:20),])

bag1 = rbind(TestDat.1,CalDat.1)
bag2 = rbind(TestDat.2,CalDat.2)

set.seed(123321)

myrf1 = randomForest(x=Tr0[,1:784], y=as.factor(as.vector(Tr0[,785])), 
                     xtest = bag1[,1:784])
claw.te1 = myrf1[["test"]][["votes"]][1:1200,1]
claw.cal1 = myrf1[["test"]][["votes"]][1201:2400,1]

myrf2 = randomForest(x=Tr0[,1:784], y=as.factor(as.vector(Tr0[,785])), 
                     xtest = bag2[,1:784])
claw.te2 = myrf2[["test"]][["votes"]][1:2000,1]
claw.cal2 = myrf2[["test"]][["votes"]][2001:4000,1]


pool.rf = randomForest(x=Tr0[,1:784], y=as.factor(as.vector(Tr0[,785])), 
                       xtest = rbind(bag1,bag2)[,1:784])
pad.te1 = pool.rf[["test"]][["votes"]][1:1200,1]
pad.cal1 = pool.rf[["test"]][["votes"]][1201:2400,1]
pad.te2 = pool.rf[["test"]][["votes"]][2401:(2400+2000),1]
pad.cal2 = pool.rf[["test"]][["votes"]][4401:(4400+2000),1]
pad.te = c(pad.te1,pad.te2)
pad.cal = c(pad.cal1,pad.cal2)


######### PROCEDURES ############
source("group_func_storey.R")
al = 0.05
theta1 = c(rep(1,120),rep(0,1080))
theta2 = c(rep(1,500),rep(0,1500))
theta = c(theta1,theta2)

#new mwthod
de.new = conf.q( c(claw.te1,claw.te2), c(claw.cal1,claw.cal2), al)
fdr.new = sum(de.new*(1-theta)) / max(1, sum(de.new) )
tp.new = sum(de.new*theta) 


#pooled Adadetect
de.pad = Adadetect(pad.te,pad.cal,al)
fdr.pad = sum(de.pad*(1-theta)) / max(1, sum(de.pad) )
tp.pad = sum(de.pad*theta) 

#separate Adadetect
de.sad1 = Adadetect(claw.te1,claw.cal1,al)
de.sad2 = Adadetect(claw.te2,claw.cal2,al)
de.sad = c(de.sad1,de.sad2)
fdr.sad = sum(de.sad*(1-theta)) /max(1, sum(de.sad) )
tp.sad = sum(de.sad*theta) 


