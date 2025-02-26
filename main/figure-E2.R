source("coordinate_simu2.R")


par(mfrow=c(3, 6), mgp=c(2, 0.5, 0), mar=c(1.5, 1.5, 1, 1))
image(thetaS, main="True States, mu=1.5")
image(thetaS.hat.new1, main="CLAW")
image(thetaS.hat.bh1, main="BH")
image(thetaS.hat.ad1, main="AdaDetect")
image(thetaS.hat.law1, main="LAWS")
image(thetaS.hat.sab1, main="SABHA")
image(thetaS, main="True States, mu=2")
image(thetaS.hat.new3, main="CLAW")
image(thetaS.hat.bh3, main="BH")
image(thetaS.hat.ad3, main="AdaDetect")
image(thetaS.hat.law3, main="LAWS")
image(thetaS.hat.sab3, main="SABHA")
image(thetaS, main="True States, mu=2.5")
image(thetaS.hat.new5, main="CLAW")
image(thetaS.hat.bh5, main="BH")
image(thetaS.hat.ad5, main="AdaDetect")
image(thetaS.hat.law5, main="LAWS")
image(thetaS.hat.sab5, main="SABHA")
