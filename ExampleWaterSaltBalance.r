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
    
# push through salt balance
foo <- watersalt_runoff(R,time,Z-vegpar$Zr,soilpar,vegpar)
      
# plot the different components
tiff("Output/water_and_saltexample.tif",width=720,height=720)
par(mfrow=c(2,1))
  with(foo,plot(time,s,type="l",xlab="time", ylab="saturation", ylim=c(0.5,1)))
  with(foo,lines(time,Svir,col="blue", lty=2, lwd=2))
  legend("topleft",c("s","Svir"),lty=c(1,2), col=c(1,"blue"))
  with(foo,plot(time,U,type="l",xlab="time", ylab="flux (cm/day)",ylim=c(-0.1,0.5)))
  with(foo,lines(time,E,col="blue", lty=2, lwd=2))
  legend("topleft",c("U","E"),lty=c(1,2), col=c(1,"blue"))
dev.off()


##########################################################################################################################################

