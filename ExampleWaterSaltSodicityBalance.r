# Test file to look at major rainfall events and their impact on the water balance
require(lattice)
source("functions/Start_up")
source("functions/Library_SaltFunctions.r")
#
set.seed(100)
# Poisson parameters derived by using Rodriguez-Iturbe [1984] for Oenpelli climate located in NT of Australia
alpha <- 1.5
lambda <- 0.4

# define a soil and vegetation type
soilpar <- Soil("SCL")      ## Checked the water balance for all other soils as coded in Librar_SaltFunctions.r. You can change soil type and can check wate balance
vegpar <- Veg("Grass",soilpar)
delta <- vegpar$delta
# Define groundwater depth
Z <- 125
# define the time 
time <- 3650 # 10 years

# use Teuling and Troch to calculate Emax
vegpar$E_max <- (1-exp(-vegpar$c_T*vegpar$LAI))*vegpar$Ep # To equate the Teuling and Troch function to the RI model
# root depth in cm
vegpar$Zr <- 25


## Generate Rain
R <- Precip(time,alpha,lambda,delta)

# --------------------------------------------------
# 1. first show Fig 2 in van der Zee et al. 2014
Conc <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3)
ESP_r <- seq(0,100)
Krel <- matrix(1,nrow=length(ESP_r),ncol=length(Conc))
for (i in 1:length(Conc)) {
  for (j in 2:length(ESP_r)) {
    Krel[j,i] <- Ksoil_Ezlit(Conc[i],ESP_r[j],Krel[j-1,i])[2]
  }
}
# plot the curves
plot(ESP_r,Krel[,1],xlab="ESP%",ylab="Relative Saturated Hydraulic conductivity",
     type="l",lwd=2)
for (i in 2:ncol(Krel)) {
  lines(ESP_r,Krel[,i],lty=i,col=i,lwd=2)  
}
lgd.txt <- paste(Conc,"molC/L")
legend("topright",lgd.txt,lty=1:7,lwd=2,col=1:7)
# this is the same as van der Zee et al. 2014 WRR
# ------------------------------------------------

# -------------------------------------------------
# 2. Salt and sodicity, no feedback
# push through salt/sodicity balance
foo <- watersaltsodic(R,time,Z-vegpar$Zr,soilpar,vegpar,FB="none")

# plot the different components
tiff("Output/sodicity_noFB_example.tif",width=720,height=720)
par(mfrow=c(2,1))
  with(foo,plot(1:time,s,type="l",xlab="time", ylab="saturation", ylim=c(0.5,1)))
  with(foo,lines(1:time,Svir,col="blue", lty=2, lwd=2))
  legend("topleft",c("s","Svir"),lty=c(1,2), col=c(1,"blue"))
  with(foo,plot(1:time,100*(1-caf),type="l",xlab="time", ylab="ESP%"))
dev.off()
# ---------------------------------------------------------------

#---------------------------------------------------------------
# 3. Salt and sodictity and feedbacks
foo_fb <-  watersaltsodic(R,time,Z-vegpar$Zr,soilpar,vegpar,FB="full")
foo_pfb <-  watersaltsodic(R,time,Z-vegpar$Zr,soilpar,vegpar,FB="partial")

# plot the different components
tiff("Output/sodicity_FB_example.tif",width=720,height=720)
par(mfrow=c(2,1))
with(foo,plot(1:time,Csalt,type="l",xlab="time", ylab="Concentration molc/L",
     ylim=c(0,0.2)))
with(foo_pfb,lines(1:time,Csalt,col="blue", lty=2, lwd=2))
with(foo_fb,lines(1:time,Csalt,col="red", lty=2, lwd=2))
legend("topleft",c("No FB","partial FB", "full FB"),lty=c(1,2,2), 
       col=c(1,"blue","red"))
with(foo,plot(1:time,100*(1-caf),type="l",xlab="time", ylab="ESP%",
     ylim=c(0,65)))
with(foo_pfb,lines(1:time,100*(1-caf),col="blue", lty=2, lwd=2))
with(foo_fb,lines(1:time,100*(1-caf),col="red", lty=2, lwd=2))
dev.off()

##########################################################################################################################################

