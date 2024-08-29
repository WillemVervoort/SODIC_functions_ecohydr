# SODIC: SODICITY, SALT AND WATER ECOHYDROLOGY MODEL
# RW Vervoort, HH Shah, SEATM van der Zee
# Version 0.9 September 2015

# Library salt functions
# This Script contains all the functions needed to run the salt balance calculations

# ------------------------------------------
# SOIL FUNCTION
# Soil data as a function
# To simplify soil input
# ########################################
Soil <- function(stype) {
  # wilting point
  psi_sh<--10
  soil <- stype
  if (soil == "L Med Clay") {
    # Medium Light Clay
    n<-0.418 # porosity
    # more soil variables for evaporation & losses
    K_s<-3.51 # cm/day
    b<-13.48 # neurotheta LMC
    nvg <- 1.089
    avg <- 0.0591
    s_fc<-0.364/n # Field capacity
    # This is the bubbling pressure
    psi_s_bar<--1.5E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    
  }
  
  if (soil == "S Clay Loam") {
    # Sandy Clay Loam
    n<-0.367 # porosity
    # more soil variables for evaporation & losses
    # Hydraulic conductivity
    K_s<-52.08 # cm/day
    # campbell's b
    b<-6.4069 # neurotheta sandy clay loam
    # van Genuchten parameters
    avg <- 0.0521
    nvg <- 1.237
    s_fc<-0.2677/n # Field capacity
    # This is the bubbling pressure
    psi_s_bar<--1.2E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
  }
  
  if (soil == "Loamy Sand") {
    # Loamy sand
    n<-0.37 # porosity
    # more soil variables for evaporation & losses
    # Hydraulic conductivity
    K_s<-175.3 # cm/day
    # Campbell's b
    b<-4.5206
    # Van Genuchten parameters
    avg <- 0.0641
    nvg <- 1.344
    s_fc<-0.2098/n # Field capacity
    
    psi_s_bar<--0.66E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    # This is the bubbling pressure
  }
  
  if (soil == "H Clay") {
    # Medium Heavy Clay
    n<-0.4473 # porosity
    # more soil variables for evaporation & losses
    # Hydraulic conductivity
    K_s<-2.82 # cm/day
    # Campbell's b
    b<-16.1501 # neurotheta Medium heavy clay
    # van Genuchten parameters
    avg <- 0.0613
    nvg <- 1.086
    s_fc<-0.3936/n # Field capacity
    # bubbling pressure
    psi_s_bar<--1.4E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
  }
  
  if (soil == "M Clay") {
    # Medium Clay
    n<-0.4391 # porosity
    # more soil variables for evaporation & losses
    # hydraulic conductivity
    K_s<-6.04 # cm/day
    # Campbell's b from neurotheta
    b<- 13.5127
    # Van Genuchten parameters
    avg <- 0.0507
    nvg <- 1.088
    s_fc<-0.3818/n # Field capacity
    # This is the bubbling pressure
    psi_s_bar<--1.75E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
  }
  if (soil == "C Sand") {
    # Coarse Sand
    n<-0.368 # porosity
    # more soil variables for evaporation & losses
    # Hydraulic conductivity
    K_s<-182.68 # cm/day
    # Campbell's b
    b<- 4.1152
    # van Genuchten parameters
    avg <- 0.0712   # Soil hydraulic parameters for van Genuchtan function
    nvg <- 1.392    # soil hydraulic papameters for van Genuchtan function
    s_fc<-0.1895/n # Field capacity
    
    psi_s_bar<--0.61E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
    # This is the bubbling pressure
  }
  if (soil == "HC") {   ## Soil types used for Sodicity paper
    # Coarse Sand
    n<-0.42 # porosity
    # more soil variables for evaporation & losses
    # Hydraulic conductivity
    K_s<-5 # cm/day
    # Campbell's b using neurotheta
    b<- 13.5
    # van Genuchten parameters
    nvg <- 1.089   # Soil hydraulic parameters for van Genuchtan function
    avg <- 0.0591  # Soil hydraulic parameters for van Genuchtan function
    s_fc<-0.364/n # Field capacity
    # This is the bubbling pressure
    psi_s_bar<--1.5E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
  }
  if (soil == "LC") {   ## Soil types used for Sodicity paper (Shah et al. 2013)
    # Coarse Sand
    n<-0.42 # porosity
    # more soil variables for evaporation & losses
    # Hydraulic conductivity
    K_s<-3.5 # cm/day
    # Campbell's b using neurotheta
    b<- 16
    # van Genuchten parameters
    nvg <- 1.089   # Soil hydraulic parameters for van Genuchtan function
    avg <- 0.0591  # Soil hydraulic parameters for van Genuchtan function
    s_fc<-0.364/n # Field capacity
    # This is the bubbling pressure
    psi_s_bar<--1.5E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
  }
  if (soil == "SCL") {    ## Soil types used for Sodicity paper (Shah et al. 2013)
    # Coarse Sand
    n<-0.42 # porosity
    # more soil variables for evaporation & losses
    # Hydraulic conductivity
    K_s<-50 # cm/day
    # Campbell's b using neurotheta
    b<- 13.5
    # van Genuchten parameters
    nvg <- 1.089   # Soil hydraulic parameters for van Genuchtan function
    avg <- 0.0591  # Soil hydraulic parameters for van Genuchtan function
    s_fc<-0.364/n # Field capacity
    # This is the bubbling pressure
    psi_s_bar<--1.5E-3 # IS IN HERE SO DON'T NEED TO DEFINE IN THE VEG FILES
  }
  
  # Other derived parameters
  s_h<-(psi_sh/psi_s_bar)^(-1/b)# soil moisture at hygroscopic point
  beta<-2*b+4 #page 714 Laio et al., 2001a
  
  # Define parameters for Eagleson function
  # Beta parameters
  beta1 <- 2+3/b
  # alpha parameter
  a1 <- 1+(3/2)/(beta1-1)
  # Create an output list
  soilpar <- list(n=n,K_s=K_s,b=b,psi_s_bar=psi_s_bar,s_fc=s_fc,s_h=s_h,beta=beta,beta1=beta1,a1=a1,nvg=nvg,avg=avg)
  return(soilpar)
}
######################################### end of soil function##################

# VEGETATION FUNCTION
######################################### start of vegetation function #########
# New vegetation input file with all vegs
Veg <- function(vtype,soilpar) {
  #ET at wilting
  E_w<-0.01 # cm/day
  k<-0.5  # index of plant resistance to water stress (page 739 second column, Porporato et al., 2001)
  if (vtype == "Grass") {
    #................................................
    # Vegetation 1 (Grass)
    # paspalum secateum F-I and R-I, 2004
    # Gives only parameters
    #.................................................
    
    Zr<-40 # soil depth (cm)   Also Table 2...Fernandez-Illescas and Rodriguez-Iturbe...2001
    delta <- 0.1  # cm interception depth (Fernandez-Illascez & Rodriguez-Iturbe, Table 2)
    
    # calculate s_w and s_star
    
    psi_ss<--0.09    ## soil matric potential at incipient stomatal closure...Also Table 2...Fernandez-Illescas and Rodriguez-Iturbe...2001
    s_star<-(psi_ss/soilpar$psi_s_bar)^(-1/soilpar$b)
    psi_sw<--4.5 # MPa   Also Table 2...Fernandez-Illescas and Rodriguez-Iturbe...2001
    
    #psi_sw<--1.5
    s_w<-(psi_sw/soilpar$psi_s_bar)^(-1/soilpar$b)
    E_max<-0.31606 # cm/day from Teuling and Troch equation with Ep=0.5cm/day,LAI=2.5,c=0.4 --->>(1-exp(-c*LAI)Ep
    
    LAI <- 2.5 ## Asner et al., 2003 Global Synthesis of LAI..implications for ecological and remote sensing studies.....  ##0.5
    c_T <- 0.4 # Tarrawarra dataset parameter in Teuling and Troch 2005
    fr <- 1  # root fraction in Teuling and Troch set to 1
    Ep <- 0.5 # cm/day BOM data average Moree  1. Also Table 2...Fernandez-Illescas and Rodriguez-Iturbe.2001..Also Teuling and Troch.2005,Vervoort and van der Zee.2008
    
    #new groundwater uptake function
    # fraction of roots close to the groundwater
    c1 <- 1.5
    
  }
  
  
  if (vtype == "Trees") {
    #......................
    # Vegetation 2 (trees)
    # proposis glandulosa F-I and R-I, 2004
    # Gives only parameters
    #............................
    
    Zr<-100 # root depth
    delta <- 0.2  # cm interception rate (Laio et al., 2001b page 748 for trees)
    
    psi_ss<--0.12  # MPa Table 1, F-I & R-I 2004
    s_star<-(psi_ss/soilpar$psi_s_bar)^(-1/soilpar$b)
    psi_sw<--2.5 # MPa (Eucalyptus Camaldulensis Whitehead and Beadle 2004)
    s_w<-(psi_sw/soilpar$psi_s_bar)^(-1/soilpar$b)
    
    E_max<-0.5 # 0.2 cm/day  Cooke et al 1998 0.3 cm/day Beadle and Whitehead 2004
    
    # Parameters for root uptake Teuling and Troch
    LAI <- 2.5 # 0.9 Huntley et al 2001 (table 1) 2.5 Whitehead and Beadle Table 4, native forest
    c_T <- 0.45 # Mid value from WAVES
    fr <- 1.0 # Set to 1 as this is unknown
    Ep <- 0.55 # cm/day BOM data average Moree
    #     Ep <- 0.9 # cm/day BOM data summer Moree
    
    #new groundwater uptake function
    # fraction of roots close to the groundwater
    c1 <- 0.8
    
    
    
  }
  return(list(E_max=E_max,E_w=E_w,k=k,Zr=Zr,delta=delta,s_w=s_w,s_star=s_star,LAI=LAI,c_T=c_T,fr=fr,Ep=Ep,c1=c1))
}
###################### end of vegetation function #############################
#############################################################################

#
# -------------------------------------
# DEFINITION OF GLOBAL FUNCTIONS
# --------------------------------------
# G function, capillary up flow
G <- function(b,hb,Z) {
  b1 <- 2+3/b
  a1 <- 1+(3/2)/(b1-1)
  H1 <- a1*(hb/Z)^b1
  return(H1)
}
# define heaviside function
H <- function(x) ifelse(x<0,0,1)


###############################################################################
# This is only the leakage part of the loss function (no EVAP)
# belonging to rho_new_1
# This is only the leackage part of the loss function (no EVAP)
# The moisture limit has been changed from s_cr to s_w because
# in the leakage function, there is no evaporation
# so we have seen that evaporation is much high as compared to the
# leakage/CR, it means there will be less loss as compared to the loss of
# the total loss function.

#******NOTE THIS WHOLE FUNCTION IS CHANGED
#******NOW SPLITS U AND L
#******************************************
L_n <- function(s,Z,soilpar,vegpar) {
  Zr <- vegpar$Zr
  hb <- soilpar$psi_s_bar*-10^4
  soilpar$s_fc <- (Z/hb)^(-1/soilpar$b)
  G1 <- G(soilpar$b,hb,Z) # using the G function 
  m <- soilpar$K_s/(soilpar$n*Zr*(exp(soilpar$beta*(1-soilpar$s_fc))-1))
  m2 <- (soilpar$K_s*G1-(soilpar$K_s*G1-vegpar$E_max)*H(
    soilpar$K_s*G1-vegpar$E_max))/(soilpar$n*Zr)
  m1 <- (soilpar$K_s*G1-(soilpar$K_s*G1-vegpar$E_max)*H(
    soilpar$K_s*G1-vegpar$E_max))/(soilpar$n*Zr*
                                     (1-exp(soilpar$beta*(vegpar$s_star-soilpar$s_fc))))
  
  eta_w<-vegpar$E_w/(soilpar$n*Zr)
  eta<-vegpar$E_max/(soilpar$n*Zr)
  # define s_cr: s-critical
  s_cr <- m2/eta*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
  # why we have taken s_w instead of s_cr(important question to ask)
  if (s <= vegpar$s_w) {
    u <-0 #*****
    l <- 0 #*****
  } else {
    if (s > vegpar$s_w && s <= vegpar$s_star)
    {
      u<--m2#*((s-vegpar$s_w)/(vegpar$s_star-vegpar$s_w))
      l <- 0 #*****
    } else {
      if (s > vegpar$s_star && s<=soilpar$s_fc) {
        u <--m1*(1-exp(soilpar$beta*(s-soilpar$s_fc)))#*****
        l <- 0 #*****
      } else {
        if (s>soilpar$s_fc && s<=1) {
          u <- 0 #*****
          l <-m*(exp(soilpar$beta*(s-soilpar$s_fc))-1)#*****
        }
      }
    }
  }
  return(list(U=u,L=l))#*****
}

################################################################################

# ---------------------------------------------------------
# Rainfall function
Precip <- function(time,alpha,lambda,delta) {
  # generate a vector of times between rainfall events > delta
  f_P<-floor(rexp(time,lambda*exp(-delta/alpha))) # vector of times between rainfall occurrences (equation 4 & 8)
  # generate a binary vector from this (impulse function)
  binary.vec <- unlist(lapply(1:time,function(i) c(rep(0,f_P[i]),1)))
  R <- rexp(length(binary.vec),1/alpha)*binary.vec 
  return(R[1:time])
}

###################################### start of E-Teuling function#############
# Evaporation function following Teuling and Troch 2005
E_Teuling <- function(s,vegpar) {
  if (s <= vegpar$s_w) {
    beta_T <- 0
  } else {
    if (s <= vegpar$s_star & s > vegpar$s_w) {
      beta_T <- (s-vegpar$s_w)/(vegpar$s_star-vegpar$s_w)
    } else {
      beta_T <- 1
    }
  }
  E <- vegpar$fr*beta_T*(1-exp(-vegpar$c_T*vegpar$LAI))*vegpar$Ep
  return(E)
}
############################# end of E-Teuling function ########################
# ----------------------------------------

#############################################
## combine water salt function for calculating numeric S due to osmotic effect
watersalt_runoff <- function(R,time,Z,soilpar,vegpar, C.in=0.02) {
  
  # R is rainfall input series
  # time is the maximum time to run the model
  # Z is the depth of the groundwater
  # soilpar is the list of soil data
  # vegpar is the list of vegetation data
  # C.in is the initial concentration of the groundwater
  deltat = 1/12
  # -----------------------------
  # Definitions
  ss<-rep(0,time*1/deltat)
  Phi<-vector(length=time*1/deltat)# Effective rainfall
  Roff<-vector(length=time*1/deltat) ## Calculates runoff 
  E<-vector(length=time*1/deltat)  # Calculates Et
  L<-vector(length=time*1/deltat) # calculates leakage loss without evaporation loss
  U<-vector(length=time*1/deltat) # ******** calculates upflow without evaporation loss
  
  Csalt<-vector(length=time*1/deltat) ## stores salt concentration values
  Sster<-vector(length=time*1/deltat) ## saturation without ET
  Svir<-vector(length=time*1/deltat) ## virtual saturation without ET
  ETvir<-vector(length=time*1/deltat) ## virtual ET
  Sosm<-vector(length=time*1/deltat) ## Saturation due to osmotic effect
  Leak<-vector(length=time*1/deltat) #** loss vector (L + U)
  # -----------------------
  # Initialisation
  Csalt[1]<-0
  newstartcon<-0
  R[length(R)+1] <- 0
  ss[1]<-soilpar$s_fc  # You need a begin value
  Svir<-soilpar$s_fc
  Sster[1] <- soilpar$s_fc #****
  Leak[1] <- 0 #*****
  b <- soilpar$b  #********* from soil the function
  h1bar <- -soilpar$psi_s_bar  #******* from the soil function
  # --------------------------------
  
  # Calculate the water balance
  for (i in 1:(length(ss)-1)) {
    # Rainfall input
    if ((round(i*deltat)+1)==(i*deltat+1)) {
      Phi[i]<-ifelse(R[round(i*deltat)+1]/(soilpar$n*vegpar$Zr)>(1-ss[i]),1-ss[i],R[round(i*deltat)+1]/(soilpar$n*vegpar$Zr))
    } else {
      Phi[i] <- 0
    }
    # Calculate runoff
    if ((round(i*deltat)+1)==(i*deltat+1)) {
      Roff[i]<-ifelse((R[round(i*deltat)+1])/(soilpar$n*vegpar$Zr)>(1-ss[i]),((R[round(i*deltat)+1]/(soilpar$n*vegpar$Zr))-(1-ss[i])),0)
    } else {
      Roff[i] <- 0
    }
    # Calculate leakage, Evaporation and upflow
    E[i]<-(do.call(E_Teuling,list(s=Svir[i],vegpar=vegpar)))
    L[i]<-(do.call(L_n,list(s=ss[i],Z=Z,soilpar=soilpar,vegpar=vegpar)))$L #**** use Z - Zr
    U[i]<-(do.call(L_n,list(s=ss[i],Z=Z,soilpar=soilpar,vegpar=vegpar)))$U #**** use Z - Zr for consistency
    # Calulate all the losses
    tot.loss <-(E[i]/(soilpar$n*vegpar$Zr)*deltat)+(L[i]*deltat)+(U[i]*deltat) #*********** NOTE U is always negative
    #     print(tot.loss)
    # Adjust the storage
    ss[i+1]<-ss[i]+Phi[i]- ifelse(tot.loss < 0,0,tot.loss) #*****************
    #     print(ss[i+1])
    # Correct the storage (moisture) for salt
    Sster[i+1]<-ss[i+1]+ (E[i]/(soilpar$n*vegpar$Zr)*deltat)
    # Caclulate teh virtual soil moisture
    Svir[i+1]<-((h1bar^(1/b))*((h1bar*(Sster[i+1])^(-b))+(3.6*Csalt[i]))^(-1/b))
    #ss[i+1]<-Svir[i+1]
    # Equation 4 in Shah et al
    Leak[i+1] <- (L[i]*deltat)+(U[i]*deltat) #**** for checking
    Csalt[i+1] <- Csalt[i] - U[i]*deltat*C.in/Sster[i]-L[i]*deltat*Csalt[i]/Sster[i]
    #*******************************************
    
  }
  
  # ----------------------
  # SUMMARISING
  days <- sort(rep(1:time,1/deltat))
  # Aggregate to daily values
  sat_out <- aggregate(ss,list(day=days),mean)
  # we are taking mean of L and ET, because we have not multiplied the E&L with deltat
  # if we multiply with deltat then we have to take sum.
  L_out <- aggregate(L,list(day=days),mean)
  U_out <- aggregate(U,list(day=days),mean)
  ET_out <- aggregate(E,list(day=days),mean)
  P_out <- aggregate(Phi*soilpar$n*vegpar$Zr,list(day=days),sum) # Note that this in only "inifltrated P" so P - Int - ROff
  # should be sum not  mean because Phi and Roff only exist in the first i
  Runoff_out <- aggregate(Roff,list(day=days),sum)
  #L_units<-L_out[,2]*(soilpar$n*vegpar$Zr)
  Csalt_out <- aggregate(Csalt,list(day=days),mean)
  Svir_out<- aggregate(Svir,list(day=days),mean)
  Sster_out<- aggregate(Sster,list(day=days),mean)
  Leak_out<- aggregate(Leak,list(day=days),sum) # error in original, should be sum not mean
  
  # return all information
  return(data.frame(time=seq(1,time),s=as.numeric(sat_out[,2]),L=as.numeric(L_out[,2]),
                    U=as.numeric(U_out[,2]),P=as.numeric(P_out[,2]),E=as.numeric(ET_out[,2]),Roff=as.numeric(Runoff_out[,2])
                    ,Csalt=as.numeric(Csalt_out[,2]),Svir=as.numeric(Svir_out[,2]),Sster=as.numeric(Sster_out[,2])
                    ,Leak=as.numeric(Leak_out[,2])))
}
##################################################
# END OSMOTIC EFFECT WB FUNCTION


####################################################
# Function that combines water salt and sodicity
# and effect on Ks
# Shah et al. WRR 2012 & van der Zee et al. 2014

watersaltsodic <- function(R,time,Z,soilpar,vegpar, 
                           C.in = 0.02, SoilCEC = 0.25,
                           fzGW = 0.05, f0 = 0.98,
                           C_soil_init = 0.00098, soilbd = 1560,
                           FB = c("none","full","partial")) {
  
  # R is rainfall input series
  # time is the maximum time to run the model
  # Z is the depth of the groundwater
  # soilpar is the list of soil data
  # vegpar is the list of vegetation data
  # C.in is the initial concentration of the groundwater
  # SoilCEC is gamma: Soil CEC(gama) = gama*(1-N)+gama*N
  # fzGW is the concentration of salt in the groundwater
  # f0 is the initial calcium fraction in the rootzone
  # C_soil_init is the initial salt concentration in the soil
  # soilbd is the soil mass in molc/kgsoil
  # FB is whether to include effects of sodicity feedback on the K_s
  
  deltat = 1/12
  # -----------------------------
  # Definitions of storage vectors
  KCESP<-vector(length=time*1/deltat)
  ss<-rep(0,time*1/deltat)
  Phi<-vector(length=time*1/deltat)# Effective rainfall
  Roff<-vector(length=time*1/deltat) ## Calculates runoff 
  E<-vector(length=time*1/deltat)  # Calculates Et
  L<-vector(length=time*1/deltat) # calculates leakage loss without evaporation loss
  U<-vector(length=time*1/deltat) # ******** calculates upflow without evaporation loss
  
  Csalt<-vector(length=time*1/deltat) ## stores salt concentration values
  Sster<-vector(length=time*1/deltat) ## saturation without ET
  Svir<-vector(length=time*1/deltat) ## virtual saturation without ET
  ETvir<-vector(length=time*1/deltat) ## virtual ET
  Sosm<-vector(length=time*1/deltat) ## Saturation due to osmotic effect
  Leak<-vector(length=time*1/deltat) #** loss vector (L + U)
  caf<-vector(length=time*1/deltat) ## Ca fraction in solution
  cax<-vector(length=time*1/deltat) ## Ca fraction in exchangeable complex
  dsdt<-vector(length=time*1/deltat) ## change in S over time
  dcdt<-vector(length=time*1/deltat) ## change in C over time
  dfdt<-vector(length=time*1/deltat) ## change in f over time
  # vector to store change in Ks and r1 factors
  Ks_corr <- vector(length=time*1/deltat)
  r1_st <- vector(length=time*1/deltat)
  # ---------------------------------------
  
  # ---------------------------------------------
  # Definition of variables
  #  variables that can be changed
  finna<-fzGW # Ca fraction coming into rootzone through capillary upflow
  gama <- SoilCEC ## Soil cation exchange capacity CEC(gama) = gama*(1-N)+gama*N
  k<-0.5 ## Gapon constant  (mol/l)^-1/2
  rho <- soilbd ## soil mass in molc/kgsoil
  b <- soilpar$b  #********* from soil the function
  h1bar <- -soilpar$psi_s_bar  #******* from the soil function
  # Store the original Ks for later use
  Ks_orig <- Ks_corr[1] <- soilpar$K_s
  
  # --- initialise vectors
  Csalt[1]<- C_soil_init # soil salt concentration
  caf[1]<- f0  ## initial calcium fraction in the rootzone
  cax[1] <- NA
  newstartcon<-0
  R[length(R)+1] <- 0
  ss[1]<-soilpar$s_fc  # You need a begin value
  Svir[1]<-soilpar$s_fc
  Sster[1] <- soilpar$s_fc #****
  Leak[1] <- 0 #*****
  # ------
  
  # Calculate the water balance
  for (i in 1:(length(ss)-1)) {
    # Rainfall input
    if ((round(i*deltat)+1)==(i*deltat+1)) {
      Phi[i]<-ifelse(R[round(i*deltat)+1]/(soilpar$n*vegpar$Zr)>(1-ss[i]),1-ss[i],R[round(i*deltat)+1]/(soilpar$n*vegpar$Zr))
    } else {
      Phi[i] <- 0
    }
    # Calculate runoff
    if ((round(i*deltat)+1)==(i*deltat+1)) {
      Roff[i]<-ifelse((R[round(i*deltat)+1])/(soilpar$n*vegpar$Zr)>(1-ss[i]),((R[round(i*deltat)+1]/(soilpar$n*vegpar$Zr))-(1-ss[i])),0)
    } else {
      Roff[i] <- 0
    }
    # Calculate leakage, Evaporation and upflow
    E[i]<-(do.call(E_Teuling,list(s=Svir[i],vegpar=vegpar)))
    
    # Include changes to Ks due to sodicity dependent on FB
    if (FB != "none") {
      # in here what happens with Ks
      soilpar$K_s <- Ks_corr[i]
      #print(Ks_corr[i])
    }
    L[i]<-(do.call(L_n,list(s=ss[i],Z=Z,soilpar=soilpar,vegpar=vegpar)))$L #**** use Z - Zr
    #print(L[i])
    # reset K_s to original if needed
    if (FB == "partial") soilpar$K_s <- Ks_orig
    U[i]<-(do.call(L_n,list(s=ss[i],Z=Z,soilpar=soilpar,vegpar=vegpar)))$U #**** use Z - Zr for consistency
    
    # Calulate all the losses
    tot.loss <-(E[i]/(soilpar$n*vegpar$Zr)*deltat)+(L[i]*deltat)+(U[i]*deltat) #*********** NOTE U is always negative
    #     print(tot.loss)
    # Adjust the storage
    ss[i+1]<-ss[i]+Phi[i]- ifelse(tot.loss < 0,0,tot.loss) #*****************
    #print(ss[i+1])
    # Correct the storage (moisture) for salt
    Sster[i+1]<-ss[i+1]+ (E[i]/(soilpar$n*vegpar$Zr)*deltat)
    # Caclulate teh virtual soil moisture
    Svir[i+1]<-((h1bar^(1/b))*((h1bar*(Sster[i+1])^(-b))+(3.6*Csalt[i]))^(-1/b))
    #print(Svir[i+1])
    #ss[i+1]<-Svir[i+1]
    # Equation 4 in Shah et al
    Leak[i+1] <- (L[i]*deltat)+(U[i]*deltat) #**** for checking
    Csalt[i+1] <- Csalt[i] - U[i]*deltat*C.in/Sster[i]-L[i]*deltat*Csalt[i]/Sster[i]
    
    ##*******Incorporate sodicity equation************************************
    # change in s
    dsdt[i]<-  ss[i+1] - ss[i]
    # change in concentration
    dcdt[i]<- Csalt[i+1] - Csalt[i]
    # calcium exchange
    cax[i+1]<-1/(1+k*(sqrt(2*Csalt[i]))*(1/sqrt(caf[i])-sqrt(caf[i])))
    ESP <- 100*(1 - caf[i])
    # change in calcium concentration
    dfdt[i]<-((-10*U[i]*finna*C.in*deltat)-(10*L[i]*caf[i]*Csalt[i]*deltat)-
                (10*soilpar$n*vegpar$Zr*caf[i]*Csalt[i]*dsdt[i])+
                (((1e-5*vegpar$Zr*rho*gama/(2*Csalt[i]))*
                (1-(1/(1+k*sqrt(2*Csalt[i])*(1/sqrt(caf[i])-sqrt(caf[i])))))*
                (1/(1+k*sqrt(2*Csalt[i])*(1/sqrt(caf[i])-sqrt(caf[i]))))-
                soilpar$n*vegpar$Zr*ss[i]*caf[i]/100)*
                (1000*dcdt[i])))/((10*soilpar$n*vegpar$Zr*ss[i]*Csalt[i])+
                (0.01*vegpar$Zr*rho*gama*k*sqrt(Csalt[i]/2)*(1/sqrt(caf[i])+
                1/(caf[i]*sqrt(caf[i])))*1/(1+k*sqrt(2*Csalt[i])*
                (1/sqrt(caf[i])-sqrt(caf[i])))^2))
    ##
    numf<-(caf[i]+dfdt[i])
    #caf[i+1]<-ifelse(numf<=2e-5 ,2e-5,numf) ## Controls Lower Limit
    if(numf<=2e-5){
      caf[i+1]<-2e-5
    } else {
      if(numf>1) {
        caf[i+1]<-1
      } else {
        caf[i+1]<-numf
      }
    }
    # Changes in K_s due to sodicity
    if (FB != "none") {
      Ks_corr[i+1] <- Ksoil_Ezlit(Csalt[i+1], ESP, Ks_corr[i])[1]
      r1_st[i+1] <- Ksoil_Ezlit(Csalt[i+1], ESP, Ks_corr[i])[2]
    }
  }
  # ----------------------
  # SUMMARISING
  days <- sort(rep(1:time,1/deltat))
  # Aggregate to daily values
  sat_out <- aggregate(ss,list(day=days),mean)
  # we are taking mean of L and ET, because we have not multiplied the E&L with deltat
  # if we multiply with deltat then we have to take sum.
  L_out <- aggregate(L,list(day=days),mean)
  U_out <- aggregate(U,list(day=days),mean)
  ET_out <- aggregate(E,list(day=days),mean)
  P_out <- aggregate(Phi*soilpar$n*vegpar$Zr,list(day=days),sum) # Note that this in only "inifltrated P" so P - Int - ROff
  # should be sum not  mean because Phi and Roff only exist in the first i
  Runoff_out <- aggregate(Roff,list(day=days),sum)
  #L_units<-L_out[,2]*(soilpar$n*vegpar$Zr)
  Csalt_out <- aggregate(Csalt,list(day=days),mean)
  Svir_out<- aggregate(Svir,list(day=days),mean)
  Sster_out<- aggregate(Sster,list(day=days),mean)
  Leak_out<- aggregate(Leak,list(day=days),sum) # error in original, should be sum not mean
  caf_out<- aggregate(caf,list(day=days),mean)
  cax_out<- aggregate(cax,list(day=days),mean)
  Ks_corr_out <- aggregate(Ks_corr,list(day=days),mean)
  r1_st_out <- aggregate(r1_st,list(day=days),mean)
  # -------------------------------
  
  # return all information
  return(list(s=as.numeric(sat_out[,2]),L=as.numeric(L_out[,2]),
              U=as.numeric(U_out[,2]),P=as.numeric(P_out[,2]),
              E=as.numeric(ET_out[,2]),Csalt=as.numeric(Csalt_out[,2]),
              Svir=as.numeric(Svir_out[,2]),Sster=as.numeric(Sster_out[,2]),
              Leak=as.numeric(Leak_out[,2]),caf=as.numeric(caf_out[,2]),
              cax=as.numeric(cax_out[,2]),Ks_corr=as.numeric(Ks_corr_out[,2]),
              r1=as.numeric(r1_st_out[,2])))
}
#############################################################
# End WB SALINITY AND SODICITY FUNCTION

###########################################################
# Ks function Ezlit
Ksoil_Ezlit <- function(Csalt_t,ESP_cal,K_prev) {
  #####________calculation of hydraulic conductivity start ....See also UNSATCHEM 16-21...also McNeal[1968] ....JULY 5 2011..
  Co <- Csalt_t*1000  ## Calculates salt concentration in mmolc/L...
  
  fmont <- 0.1 ## weight fraction montmorillonite in the soil...See also UNSATCHEM 16-21...also McNeal[1968]
  d_ster <- ifelse(Co > 300,0,356.4*Co ^(-0.5)+1.2) ##....See also UNSATCHEM 16-21...also McNeal[1968]
  ESP_ster <- max(0,ESP_cal-(1.24+11.63*log10(Co)))   ## calculates the ESP adjusted...See also UNSATCHEM 16-21...also McNeal[1968]
  
  swell <- fmont*3.6*10^(-4)*ESP_ster*d_ster  ## calculates the swelling factor....See also UNSATCHEM 16-21...also McNeal[1968]
  ## calculats the nn for the r1 calculation as an input...water retention parameters... .See also UNSATCHEM 16-21...also McNeal[1968]
  nn <- ((ESP_cal/100)^0.449+1.005)  ## Equation for n(ESP) developed by Ezlit_2010 at page 124..equation 5.24 
  cc <- 0.846*exp(10.967*ESP_cal/100) ## Equation for c(ESP) developed by Ezlit_2010 at page 124..equation 5.25
  r1 <- 1 - (cc*swell^nn)/(1+cc*swell^nn) ## calculates the relative hydraulic conductivity....See also UNSATCHEM 16-21...also McNeal[1968]
  ### This checks that if Hydraulic conductivity increases or decreases..If decreases that is Acceptable..
  ##If it increases..This is not acceptable..In this case keep the K same as in the previous time step..
  K_new <- ifelse(K_prev*r1 > K_prev, K_prev, K_prev*r1)
  return(c(K_new,r1))
}






