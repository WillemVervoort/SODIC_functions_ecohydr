#
# Script to load functions for RI model with groundwater
# Load at the beginning of a session
# 
# September 2008
# variation April 2012 (new DR functions)
# Willem Vervoort WIMEK project
# ################################################

# This now contains the following
# 1. calls to soilfunction.r and vegfunction.r  (line 16ff)
# 2. calls to Lattice and splancs (l 22ff)
# 2a. Dirac delta, heaviside and G function  (l 33 ff)
# 2b. Root weighting functions (l 56 ff)
# 2c. Rainfall function
# 3. definition of the old and new loss functions (rho) (l 76 ff)
# 3a. Definition of the root loss function  (l 202 ff)
# 4. definition of the function Mpar() (l 308 ff)
# 5. definition of the pdfs based on the loss functions (l 338 ff)
# 6. Integrals of the pdfs (l 399 ff)
# 7. Water stress functions (to be expanded)
# 8. Clusgen function
# 9. Some older functions (not really needed)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1. load the soil function
source("soilfunction.r")
# load veg function
source("vegfunction.r")

# 2. See Rowlingson and Diggle, 1993 Computers and Geosciences 19:627-655
#require(splancs)
## Original plotting was all using lattice
require(lattice)

# -------------------------------------
# DEFINITION OF GLOBAL FUNCTIONS
# --------------------------------------
# 2a
# The Dirac delta function
Dirac<-function(x){ifelse(x==0,1,0)}

# define heaviside function
H <- function(x) ifelse(x<0,0,1)


# H function
G <- function(b,hb,Z) {
  b1 <- 2+3/b
  a1 <- 1+(3/2)/(b1-1)
  H1 <- a1*(hb/Z)^b1
  return(H1)
}

# The root weighting function
# define function U (Feddes et al 2001)
# Deep root functions              #VvdZ 2009 eq.7
f <- function(x,c1) exp(-c1*x/100)
U <- function(z1,z2,c1) integrate(f,z2,z1,c1)$value/integrate(f,0,Inf,c1)$value   #eq.8
# root water uptake function
s.fun <- function(x,fs)exp(-fs*x)
#  Rc_B <- function(z1,c1,z2) U(z1,z2=z2,Z=z1,c1=c1)*z2/(z1-z2)                 #eq.7

# mean root density in top
#f.Zr <- function(Zr,c1) 1/Zr*integrate(f,0,Zr,c1)$value
# mean root density in bottom
#f.DR <- function(Zr,Z,c1) 1/(Z-Zr)*integrate(f,Zr,Z,c1)$value

# take this out, move to mpar
#Rc_B <- function(z1,c1,z2) U(z1,z2=z2,Z=z1,c1=c1)*z1/(z1-vegpar$Zr)
RWU <- function(z1,Zmean,fs) 2*dnorm(z1/100,mean=Zmean/100,s.fun(Zmean/100,fs))# test values
Rc_B <- function(z1,z2,c1,Zmean,fs) RWU(z1=z1,Zmean=Zmean,fs=fs)*U(z1=z1,z2=z2,c1=c1)

#---------------------------------------------------------------------------
# 2c Rainfall function (Poisson)
# --------------------------------------------------------
# Rainfall function
Rain.p <- function(time,alpha,lambda,delta) {
  # generate a vector of times between rainfall events > delta
  f_P<-floor(rexp(time,lambda*exp(-delta/alpha))) # vector of times between rainfall occurrences (equation 4 & 8)
  # generate a binary vector from this (impulse function)
  binary.vec <- unlist(lapply(1:time,function(i) c(rep(0,f_P[i]),1)))
  R <- rexp(length(binary.vec),1/alpha)*binary.vec 
  return(R[1:time])
}




# A temporary definition of s
s <- seq(0,1,length=200)
# ---------------------------------------------------------------------
# 3. Losses: Overall loss function (equation 16 Laio et al., 2001a)
# -----------------------------------------------------------------------
rho<-function(s,Z,soilpar,vegpar) {
    eta_w <-vegpar$E_w/(soilpar$n*vegpar$Zr)
    eta <- vegpar$E_max/(soilpar$n*vegpar$Zr)

     if (s <= soilpar$s_h) {
         r<-0
     } else {
         if (s > soilpar$s_h && s <= vegpar$s_w)
         {
            r<- eta_w*((s-soilpar$s_h)/(vegpar$s_w-soilpar$s_h))
         } else {
            if (s > vegpar$s_w && s <= vegpar$s_star)
            {
               r <- eta_w+(eta-eta_w)*((s - vegpar$s_w)/(vegpar$s_star-vegpar$s_w))
            } else {
               if (s > vegpar$s_star && s <= soilpar$s_fc) {
                  r <- eta
               } else {
                  if (s > soilpar$s_fc && s <= 1) {
                      r <- eta + soilpar$K_s/(soilpar$n*vegpar$Zr*(exp(soilpar$beta*(1-soilpar$s_fc))-1))*(exp(soilpar$beta*(s - soilpar$s_fc))-1)
                  }
               }
            }
         }
     }
  return(r)
}



# ...........................................
# New loss functions
#This version has 3 limits, VvdZ equation 9 & equation 11
rho_new_1<-function(s,Z,soilpar,vegpar) {
      Zr <- vegpar$Zr
      hb <- soilpar$psi_s_bar*-10^4
      soilpar$s_fc <- (Z/hb)^(-1/soilpar$b) # this is s_lim
      G1 <- G(soilpar$b,hb,Z) # using the G function
      m <- soilpar$K_s/(soilpar$n*Zr*(exp(soilpar$beta*(1-soilpar$s_fc))-1))
      m2 <- soilpar$K_s*G1/(soilpar$n*Zr)
      m1 <- soilpar$K_s*G1/(soilpar$n*Zr*
          (1-exp(soilpar$beta*(vegpar$s_star-soilpar$s_fc))))
          #(exp(soilpar$beta*(soilpar$s_fc-vegpar$s_star))-1))
          
      eta<-vegpar$E_max/(soilpar$n*Zr)

      r <- 0
 
      if (s > soilpar$s_fc) {
        r<-eta+m*(exp(soilpar$beta*(s-soilpar$s_fc))-1)
      }
      if (m1 < eta) {
      # define s_cr: s-critical
        s_cr <- m2/eta*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
          if (s > s_cr && s<=vegpar$s_star)
          {
            r<-(eta-m2)*((s-s_cr)/(vegpar$s_star-s_cr))
          } else {
            if (s>vegpar$s_star && s<=soilpar$s_fc) {
               r<-eta-m1*(1-exp(soilpar$beta*(s-soilpar$s_fc))) #(exp(soilpar$beta*(soilpar$s_fc-s))-1)
            }
          }
      } else {
          if (s > vegpar$s_star && s <= soilpar$s_fc) {
            r <- eta*exp(soilpar$beta*(s-soilpar$s_fc))
          }
      }
      return(r)
}


# This version has 4 limits, equation 15 VvdZ (dry end)
rho_new_2<-function(s,Z,soilpar,vegpar) {
      Zr <- vegpar$Zr
      hb <- soilpar$psi_s_bar*-10^4
      soilpar$s_fc <- (Z/hb)^(-1/soilpar$b)
      G1 <- G(soilpar$b,hb,Z) # using the G function
      m <- soilpar$K_s/(soilpar$n*Zr*(exp(soilpar$beta*(1-soilpar$s_fc))-1))
      m2 <- (soilpar$K_s*G1-vegpar$E_w-(soilpar$K_s*G1-vegpar$E_max)*H(
          soilpar$K_s*G1-vegpar$E_max))/(soilpar$n*Zr)
      m1 <- (soilpar$K_s*G1-vegpar$E_w-(soilpar$K_s*G1-vegpar$E_max)*H(
          soilpar$K_s*G1-vegpar$E_max))/(soilpar$n*Zr*
          (1-exp(soilpar$beta*(vegpar$s_star-soilpar$s_fc))))
          #(exp(soilpar$beta*(soilpar$s_fc-vegpar$s_star))-1))
          
      eta_w<-vegpar$E_w/(soilpar$n*Zr)
      eta<-vegpar$E_max/(soilpar$n*Zr)
      r <- 0
      if (s > soilpar$s_fc) {
       r<-eta+m*(exp(soilpar$beta*(s-soilpar$s_fc))-1)
      }
      if (m2 < eta_w) {
        s_cr <- m2/eta_w*(vegpar$s_w-soilpar$s_h)+soilpar$s_h
        if (s <= vegpar$s_w && s > s_cr) {
          s <- (eta_w - m2)*(s - s_cr)/(vegpar$s_w - s_cr)
        } else {
          if (s <= vegpar$s_star && s > vegpar$s_w) {
            s <- (eta_w - m2) + (eta - eta_w)*(s - vegpar$s_w)/(vegpar$s_star - vegpar$s_w)
          } else {
            if (s > vegpar$s_star && s <= soilpar$s_fc) {
               r <- eta-m1*(1-exp(soilpar$beta*(s-soilpar$s_fc))) #(exp(soilpar$beta*(soilpar$s_fc-s))-1)
            }
          }
        }
      } else {
        if (m1 < eta) {
          s_cr <- m2/eta*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
          if (s > s_cr && s<=vegpar$s_star)
          {
            r<-(eta-m2)*((s-s_cr)/(vegpar$s_star-s_cr))
          } else {
            if (s>vegpar$s_star && s<=soilpar$s_fc) {
               r<-eta-m1*(1-exp(soilpar$beta*(s-soilpar$s_fc))) #(exp(soilpar$beta*(soilpar$s_fc-s))-1)
            }
          }
        } else {
          if (s > vegpar$s_star && s <= soilpar$s_fc) {
            r <- eta*exp(soilpar$beta*(s-soilpar$s_fc))
          }
        }
      }
      return(r)
}

# This is the function with deep roots
# inserted 20080830
# loss function including deep root function
rho_root <- function(s,Z,soilpar,vegpar) {
      fr <- U(vegpar$Zr,0,Z,vegpar$c1)
      Zr <- vegpar$Zr
      hb <- soilpar$psi_s_bar*-10^4
      soilpar$s_fc <- (Z/hb)^(-1/soilpar$b)
      G1 <- G(soilpar$b,hb,Z) # using the G function
      m <- soilpar$K_s/(soilpar$n*Zr*(exp(soilpar$beta*(1-soilpar$s_fc))-1))
      m2 <- soilpar$K_s*G1/(soilpar$n*Zr)
     m1 <- (soilpar$K_s*G1-vegpar$E_w-(soilpar$K_s*G1-vegpar$E_max)*H(
          soilpar$K_s*G1-vegpar$E_max))/(soilpar$n*Zr*
          (1-exp(soilpar$beta*(vegpar$s_star-soilpar$s_fc))))
   
      eta<-vegpar$E_max/(soilpar$n*Zr)
      # define s_cr: s-critical
      s_cr <- m2/(eta)*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
      s_cr.r <- (1 - Rc_B(Z,vegpar$c1,Zr))/fr*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
      r <- 0

      if (s > soilpar$s_fc) {
        r<-(1 - Rc_B(Z,vegpar$c1,Zr))*(eta+m*(exp(soilpar$beta*(s-soilpar$s_fc))-1))
      }
      if (s > s_cr && s <= s_cr.r) r<- fr*(eta-m2)*(s - s_cr)/(s_cr.r - s_cr)
      if (s > s_cr.r && s <= vegpar$s_star) {
         r <- (1 - Rc_B(Z,vegpar$c1,Zr))*(eta-m2)
       }
      if (s > vegpar$s_star && s <= soilpar$s_fc) {
                r <- (1 - Rc_B(Z,vegpar$c1,Zr))*eta - m1*(1-exp(soilpar$beta*(s-soilpar$s_fc)))
      }
      return(r)
  }


# This is only the leakage part of the loss function (no EVAP)
# belonging to rho_new_1
L_n <- function(s,Z,soilpar,vegpar) {
      Zr <- vegpar$Zr
      hb <- soilpar$psi_s_bar*-10^4
      soilpar$s_fc <- (Z/hb)^(-1/soilpar$b)
      G1 <- G(soilpar$b,hb,Z) # using the G function
      m <- soilpar$K_s/(soilpar$n*Zr*(exp(soilpar$beta*(1-soilpar$s_fc))-1))
      m2 <- soilpar$K_s*G1/(soilpar$n*Zr)
      p <- (1-exp(soilpar$beta*(vegpar$s_star-soilpar$s_fc)))
      m1 <- m2/p
#          (exp(soilpar$beta*(soilpar$s_fc-vegpar$s_star))-1))
      eta_w<-vegpar$E_w/(soilpar$n*Zr)
      eta<-vegpar$E_max/(soilpar$n*Zr)
      # define s_cr: s-critical
      s_cr <- m2/eta*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
     
     if (s <= s_cr) {
         l<-0
     } else {
        if (s > s_cr && s <= vegpar$s_star)
        {
           l<--m2*((s-s_cr)/(vegpar$s_star-s_cr))
            } else {
               if (s > vegpar$s_star && s<=soilpar$s_fc) {
                  l <--m1*(1-exp(soilpar$beta*(s-soilpar$s_fc)))#(exp(soilpar$beta*(soilpar$s_fc-s))-1)
               } else {
                  if (s>soilpar$s_fc && s<=1) {
                      l<-m*(exp(soilpar$beta*(s-soilpar$s_fc))-1)
                  }
               }
         }
     }
  return(l)
}


# -------------------------------------------------------
# 4. define a function with all the model variables
# ----------------------------------------------------------
Mpar <- function(Z,Zmean,soilpar,vegpar,alpha,lambda) {
      l_prime <- lambda*exp(-vegpar$delta/alpha)
      gam_p<-soilpar$n*vegpar$Zr/alpha
      Zr <- vegpar$Zr
      hb <- soilpar$psi_s_bar*-10^4
      G1 <- do.call(G,list(b=soilpar$b,hb=hb,Z=Z-Zr)) # using the G function
      # define s_lim
      s_lim <-((Z-Zr)/hb)^(-1/soilpar$b)
#      s_lim <- soilpar$s_fc

      m <- soilpar$K_s/(soilpar$n*Zr*(exp(soilpar$beta*(1-soilpar$s_fc))-1))
      m_n <- soilpar$K_s/(soilpar$n*Zr*(exp(soilpar$beta*(1-s_lim))-1))
      m2 <- soilpar$K_s*G1/(soilpar$n*Zr)
      m1 <- m2/(1-exp(soilpar$beta*(vegpar$s_star-s_lim)))

      eta_w<-vegpar$E_w/(soilpar$n*Zr)
      eta<-vegpar$E_max/(soilpar$n*Zr)
      # define s_cr: s-critical
      # Note there are three versions of s_cr depending on what function you use
      # 3 limits only capil
      s_cr <- m2/eta*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
      # 4 limits only capil
      s_cr2 <- m2/eta_w*(vegpar$s_w-soilpar$s_h)+soilpar$s_h
      # 3 limits roots and capil.
      # insert new root functions WV 20120726
      fr <- U(vegpar$Zr,0,vegpar$c1)
      Rc <- Rc_B(Z,vegpar$c1,vegpar$Zr,Zmean,vegpar$fs)
   #   s_cr.r2 <-  m2/eta*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
   # variation to limit to lower than s_lim WV 26042012
      temp <- m2/eta*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
      s_cr.r2 <- ifelse(temp>s_lim,s_lim,temp)
      # test s_cr.r to vegpar$s_w
      #s_cr.r <- (1 - Rc)/fr*(vegpar$s_star-vegpar$s_w)+ vegpar$s_w
      s_cr.r <- (1 - Rc)/fr*(vegpar$s_star-vegpar$s_w)+vegpar$s_w
      mpar <- list(hb=hb,m=m,m_n=m_n,m1=m1,m2=m2,eta=eta,eta_w=eta_w,
                s_cr=s_cr,s_cr2=s_cr2,s_cr.r=s_cr.r,s_cr.r2=s_cr.r2,
                s_lim=s_lim,l_prime=l_prime,gam_p=gam_p,G1=G1,Rc=Rc,fr=fr)
return(mpar)
}

# -----------------------------------------------------------------
# 5. Probability density functions
# Based on Laio et al.
# and following Mathematica script from Amilcare Porporato
# --------------------------------------------------------------------
# s_h < s <= s_w
C<-1
pdf_s1<-function(s=s,soilpar,vegpar,mpar)
{
     ps<-exp(-mpar$gam_p*s)*C/mpar$eta_w*((s-soilpar$s_h)/(vegpar$s_w-soilpar$s_h))^
          (mpar$l_prime*(vegpar$s_w-soilpar$s_h)/mpar$eta_w-1)
     return(ps)
}
# s_w < s <= s*
pdf_s2<-function(s=s,soilpar,vegpar,mpar)
{
     ps<-exp(-mpar$gam_p*s)*C/mpar$eta_w*(1+(mpar$eta/mpar$eta_w-1)*(s-vegpar$s_w)/
          (vegpar$s_star-vegpar$s_w))^(mpar$l_prime*(vegpar$s_star-vegpar$s_w)/(mpar$eta-mpar$eta_w)-1)
     return(ps)
}

# s* < s <= s_fc
pdf_s3<-function(s=s,soilpar,vegpar,mpar)
{
     ps<-exp(-mpar$gam_p*s+(mpar$l_prime/mpar$eta)*(s-vegpar$s_star))*C/mpar$eta*
          (mpar$eta/mpar$eta_w)^(mpar$l_prime*(vegpar$s_star-vegpar$s_w)/(mpar$eta-mpar$eta_w))
     return(ps)
}

# s_fc < s <= 1
pdf_s4<-function(s=s,soilpar,vegpar,mpar)
{
      ps<-C/mpar$eta*exp(soilpar$beta*soilpar$s_fc)*exp(-(soilpar$beta+mpar$gam_p)*s
          )*exp(mpar$l_prime/mpar$eta*(soilpar$s_fc-vegpar$s_star))*(
          mpar$eta*exp(soilpar$beta*s)/((mpar$eta-mpar$m)*exp(soilpar$beta*soilpar$s_fc
          )+mpar$m*exp(soilpar$beta*s)))^(mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m)
          )+1)*(mpar$eta/mpar$eta_w)^(mpar$l_prime*(vegpar$s_star-vegpar$s_w)/(mpar$eta-mpar$eta_w))
      return(ps)
}

# .....................................................
# New functions with groundwater uptake: equations 10, 13 and 16 in VvdZ
# ....................................................
# s_h < s <= s_cr
# This part of the pdf does not exist in the simple model (rho_new_1) only for rho_new_2
pdf_s1_n2 <- function(s,soilpar,vegpar,mpar)
{
  if (mpar$eta_w < mpar$m2) {
    ps <- rep(0,length(s))
  } else {
    ps <- C/(mpar$eta_w-mpar$m2)*((s-mpar$s_cr2)/
    (vegpar$s_w-mpar$s_cr2))^(mpar$l_prime*(vegpar$s_w-mpar$s_cr2)/
    (mpar$eta_w-mpar$m2)-1)*exp(-mpar$gam_p*s)
  }
  return(ps)
}


# s_cr < s <= s*
# rho_new_1
pdf_s2_n<-function(s,soilpar,vegpar,mpar)
{
   if (mpar$eta-mpar$m1<=0) {
     ps <- rep(0,length(s))
   } else {
     B <- (mpar$l_prime*(vegpar$s_star-mpar$s_cr)/(mpar$eta-mpar$m2)-1)
     ps<-exp(-mpar$gam_p*s)*C/(mpar$eta-mpar$m2)*((s-mpar$s_cr)/(vegpar$s_star-mpar$s_cr))^B
   }
  return(ps)
}

# rho_new_2
pdf_s2_n2<-function(s,soilpar,vegpar,mpar)
{
   if (mpar$eta_w > mpar$m2) {
      ps <- C/(mpar$eta-mpar$m2)*exp(-mpar$gam_p*s)(1+(mpar$eta-mpar$eta_w)/(mpar$eta_w-mpar$m2)*
            (s-vegpar$s_w)/(vegpar$s_star-vegpar$s_w))^(l_prime*(vegpar$s_star-vegpar$s_w)/(mpar$eta-mpar$eta_w)-1)
   } else {
      if (mpar$eta > mpar$m1) {
       B <- (mpar$l_prime*(vegpar$s_star-mpar$s_cr)/(mpar$eta-mpar$m2)-1)
       ps<-exp(-mpar$gam_p*s)*C/(mpar$eta-mpar$m2)*((s-mpar$s_cr)/(vegpar$s_star-mpar$s_cr))^B
      } else {
       ps <- rep(0,length(s))
      }
   }
   return(ps)
}
 


# s* < s <= s_fc                         
#rho_new_1
pdf_s3_n<-function(s,soilpar,vegpar,mpar)
{
  if (mpar$eta-mpar$m1 <= 0 ) {
     ps <- C/mpar$eta*exp(-soilpar$beta*(s-mpar$s_lim))*exp(-mpar$gam_p*s)*
        exp(mpar$l_prime/(soilpar$beta*mpar$eta)*(1-exp(-soilpar$beta*(s-mpar$s_lim))))
  } else {   
     A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime/(mpar$eta-mpar$m1)*(s-vegpar$s_star))
     B <- (mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m1)))
     A <- (1-mpar$m2/mpar$eta)^B
     ps<-C/mpar$eta*A1*A*(1-mpar$m1/mpar$eta*(1-exp(soilpar$beta*(s-mpar$s_lim))))^(-B-1)
  }
  return(ps)          
}

#rho_new_2
pdf_s3_n2<-function(s,soilpar,vegpar,mpar)
{
     if (mpar$eta_w > mpar$m2) {
       A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime*(s-vegpar$s_star)/(mpar$eta-mpar$m1))
       B <- (mpar$l_prime*(vegpar$s_star-vegpar$s_w)/(mpar$eta-mpar$eta_w))
       A <- (1+(mpar$eta-mpar$eta_w)/(mpar$eta_w-mpar$m2))^B
       ps<-C/mpar$eta*(1-mpar$m1/mpar$eta*(1-exp(soilpar$beta*(s-mpar$s_lim)))
              )^(-mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m1))-1)*
              (1-mpar$m2/mpar$eta)^(mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m1)))*A1*A
     } else {
        if (mpar$eta > mpar$m1) {
         A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime/(mpar$eta-mpar$m1)*(s-vegpar$s_star))
         B <- (mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m1)))
         A <- (1-mpar$m2/mpar$eta)^B
         ps<-C/mpar$eta*A1*A*(1-mpar$m1/mpar$eta*(1-exp(soilpar$beta*(s-mpar$s_lim))))^(-B-1)
        } else {
         ps <- C/mpar$eta*exp(-soilpar$beta*(s-mpar$s_lim))*exp(-mpar$gam_p*s)*
            exp(mpar$l_prime/(soilpar$beta*mpar$eta)*(1-exp(-soilpar$beta*(s-mpar$s_lim))))
        }
     }
  return(ps)
}

# s_fc < s <= 1
# rho_new_1
pdf_s4_n<-function(s,soilpar,vegpar,mpar)
{
  if (mpar$eta-mpar$m1 <=0) {
      ps<-C/mpar$eta*exp(soilpar$beta*mpar$s_lim)*exp(-(soilpar$beta+mpar$gam_p)*s
      )*(mpar$eta*exp(soilpar$beta*s)/((mpar$eta-mpar$m_n)*exp(soilpar$beta*mpar$s_lim
      )+mpar$m_n*exp(soilpar$beta*s)))^(mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m_n)
      )+1)
  } else {
     B1 <- (mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m1)))
     A <-(1-mpar$m2/mpar$eta)^B1
     A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime/(mpar$eta-mpar$m_n)*(s-mpar$s_lim))*
           exp(-mpar$l_prime/(mpar$eta-mpar$m1)*(vegpar$s_star-mpar$s_lim))
     ps<- C/mpar$eta*(1+mpar$m_n/mpar$eta*(exp(soilpar$beta*(s-mpar$s_lim))-1))^
          (-mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m_n))-1)*A*A1
  }
     return(ps)
}

# rho_new_2
pdf_s4_n2<-function(s,soilpar,vegpar,mpar)
{
   if (mpar$eta_w > mpar$m2) {
       A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime*(mpar$s_lim-vegpar$s_star)/(mpar$eta-mpar$m1))*
            exp(-mpar$l_prime/(mpar$eta-mpar$m_n)*(s-mpar$s_lim))
       B <- (mpar$l_prime*(vegpar$s_star-vegpar$s_w)/(mpar$eta-mpar$eta_w))
       A <- (1+(mpar$eta-mpar$eta_w)/(mpar$eta_w-mpar$m2))^B
       ps<-C/mpar$eta*(1+mpar$m_n/mpar$eta*(exp(soilpar$beta*(s-mpar$s_lim))-1)
              )^(-mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m_n))-1)*
              (1-mpar$m2/mpar$eta)^(mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m1)))*A1*A
   } else {
     if (mpar$eta < mpar$m1) {
       B1 <- (mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m1)))
       A <-(1-mpar$m2/mpar$eta)^B1
       A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime/(mpar$eta-mpar$m_n)*(s-mpar$s_lim))*
             exp(-mpar$l_prime/(mpar$eta-mpar$m1)*(vegpar$s_star-mpar$s_lim))
       ps<- C/mpar$eta*(1+mpar$m_n/mpar$eta*(exp(soilpar$beta*(s-mpar$s_lim))-1))^
            (-mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m_n))-1)*A*A1
     } else  {
        ps<-C/mpar$eta*exp(soilpar$beta*mpar$s_lim)*exp(-(soilpar$beta+mpar$gam_p)*s
        )*(mpar$eta*exp(soilpar$beta*s)/((mpar$eta-mpar$m_n)*exp(soilpar$beta*mpar$s_lim
        )+mpar$m_n*exp(soilpar$beta*s)))^(mpar$l_prime/(soilpar$beta*(mpar$eta-mpar$m_n)
        )+1)
     }
   }
  return(ps)
}

# ...................................................................
# 5a Root pdfs  3 limits for model RF
# inserted 20080830 and adjusted 20090428
# ....................................................................
# between s_cr and s_cr.r
pdf_s1_r <- function(s,soilpar,vegpar,mpar,Z)
{
  if (mpar$fr*(mpar$eta-mpar$m2) <= 0 ) {
     ps <- rep(0,length(s))
  } else {
    ps<- ifelse(s - mpar$s_cr.r2>0,C/(mpar$fr*(mpar$eta-mpar$m2))*exp(-mpar$gam_p*s)*
          (vegpar$s_star - mpar$s_cr.r2)/(s - mpar$s_cr.r2)*
           ((s - mpar$s_cr.r2)/(mpar$s_cr.r - mpar$s_cr.r2))^(mpar$l_prime*(vegpar$s_star - 
             mpar$s_cr.r2)/(mpar$fr*(mpar$eta-mpar$m2))),0)
  }
  return(ps)
}


# between s_cr.r and s_star
pdf_s2_r <- function(s,soilpar,vegpar,mpar,Z)
{
  if ((1-mpar$Rc)*(mpar$eta-mpar$m2) <= 0 ) {
     ps <- rep(0,length(s))
  } else {
    ps<-C/((1-mpar$Rc)*(mpar$eta-mpar$m2))*exp(-mpar$gam_p*s)*
        exp(mpar$l_prime*(s-mpar$s_cr.r)/((1-mpar$Rc)*(mpar$eta-mpar$m2)))
  }
  return(ps)
}

# between s_star and s_lim
pdf_s3_r <- function(s,soilpar,vegpar,mpar,Z)
{
  if ((1-mpar$Rc)*(mpar$eta-mpar$m1) <= 0 ) {
     ps <- C/(mpar$eta*(1-mpar$Rc))*exp(-soilpar$beta*(s-mpar$s_lim))*exp(-mpar$gam_p*s)*
        exp(mpar$l_prime/(soilpar$beta*mpar$eta*(1-mpar$Rc))*(1-exp(-soilpar$beta*(s-mpar$s_lim))))
  } else {   
     A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime/((1-mpar$Rc)*(mpar$eta-mpar$m1))*(s-vegpar$s_star))
     B <- mpar$l_prime/(soilpar$beta*((1-mpar$Rc)*(mpar$eta-mpar$m1)))
     A <- (1-mpar$m2/mpar$eta)^B
     B1 <- exp(mpar$l_prime*(vegpar$s_star-mpar$s_cr.r)/((1-mpar$Rc)*(mpar$eta-mpar$m2)))
     ps<-C/(mpar$eta*(1-mpar$Rc))*A1*A*B1*(1-mpar$m1/mpar$eta*(1-exp(soilpar$beta*(s-mpar$s_lim))))^(-B-1)
  }
  return(ps)          
}

# between s_lim and 1
pdf_s4_r <- function(s,soilpar,vegpar,mpar,Z)
{
  if ((1-mpar$Rc)*(mpar$eta-mpar$m1) <=0) {
      B <- -mpar$l_prime/(soilpar$beta*((1-mpar$Rc)*mpar$eta-mpar$m_n))
      A <- exp(mpar$l_prime/((1-mpar$Rc)*(mpar$eta-mpar$m_n))*(s-mpar$s_lim))*exp(-mpar$gam_p*s)
      ps<-C/(mpar$eta*(1-mpar$Rc))*(1+mpar$m_n/(mpar$eta)*A*(exp(soilpar$beta*(s-mpar$s_lim))-1))^(B-1)
  } else {
     B <- mpar$l_prime/(soilpar$beta*(1-mpar$Rc)*(mpar$eta-mpar$m1))
     A <-(1-mpar$m2/mpar$eta)^B
     A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime/((1-mpar$Rc)*mpar$eta-mpar$m_n)*(s-mpar$s_lim))*
           exp(-mpar$l_prime/((1-mpar$Rc)*(mpar$eta-mpar$m1))*(vegpar$s_star-mpar$s_lim))
     B1 <- exp(mpar$l_prime*(vegpar$s_star-mpar$s_cr.r)/((1-mpar$Rc)*(mpar$eta-mpar$m2)))
     ps<- C/((1-mpar$Rc)*mpar$eta)*(1+mpar$m_n/((1-mpar$Rc)*mpar$eta)*(exp(soilpar$beta*(s-mpar$s_lim))-1))^
          (-mpar$l_prime/(soilpar$beta*((1-mpar$Rc)*mpar$eta-mpar$m_n))-1)*A*A1*B1
  }
return(ps)
}
# ...................................................................
# 5c Root pdfs  3 limits for model FB and model E
# inserted 20080830 and adjusted 20090428
# ....................................................................
# between s_cr and s_star
pdf_s2_fr <- function(s,soilpar,vegpar,mpar)
{
   if (mpar$fr*(mpar$eta-mpar$m2)<=0) {
     ps <- rep(0,length(s))
   } else {
     B <- (mpar$l_prime*(vegpar$s_star-mpar$s_cr.r2)/(mpar$fr*(mpar$eta-mpar$m2))-1)
     ps<-exp(-mpar$gam_p*s)*C/(mpar$fr*(mpar$eta-mpar$m2))*((s-mpar$s_cr.r2)/(vegpar$s_star-mpar$s_cr.r2))^B
   }
  return(ps)
}

# between s_star and s_lim
pdf_s3_fr<-function(s,soilpar,vegpar,mpar)
{
  if (mpar$fr*(mpar$eta-mpar$m1) <= 0 ) {
     ps <- C/(mpar$fr*mpar$eta)*exp(-soilpar$beta*(s-mpar$s_lim))*exp(-mpar$gam_p*s)*
        exp(mpar$l_prime/(soilpar$beta*mpar$fr*mpar$eta)*(1-exp(-soilpar$beta*(s-mpar$s_lim))))
  } else {   
     A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime/(mpar$fr*(mpar$eta-mpar$m1))*(s-vegpar$s_star))
     B <- (mpar$l_prime/(soilpar$beta*mpar$fr*(mpar$eta-mpar$m1)))
     A <- (1-mpar$m2/mpar$eta)^B
     ps<-C/(mpar$fr*mpar$eta)*A1*A*(1-mpar$m1/mpar$eta*(1-exp(soilpar$beta*(s-mpar$s_lim))))^(-B-1)
  }
  return(ps)          
}

# between s_lim and 1
pdf_s4_fr <- function(s,soilpar,vegpar,mpar)
{
  if (mpar$fr*(mpar$eta-mpar$m1) <=0) {
      B <- -mpar$l_prime/(soilpar$beta*(mpar$fr*mpar$eta-mpar$m_n))
      A <- exp(mpar$l_prime/(mpar$fr*mpar$eta-mpar$m_n)*(s-mpar$s_lim))*exp(-mpar$gam_p*s)
      ps<-C/(mpar$fr*mpar$eta)*(1+mpar$m_n/(mpar$fr*mpar$eta))*A*(exp(soilpar$beta*(s-mpar$s_lim))-1)^(B-1)
  } else {
     B1 <- (mpar$l_prime/(soilpar$beta*(mpar$fr*(mpar$eta-mpar$m1))))
     A <-(1-mpar$m2/mpar$eta)^B1
     A1 <- exp(-mpar$gam_p*s)*exp(mpar$l_prime/(mpar$fr*mpar$eta-mpar$m_n)*(s-mpar$s_lim))*
           exp(-mpar$l_prime/(mpar$fr*(mpar$eta-mpar$m1))*(vegpar$s_star-mpar$s_lim))
     ps<- C/(mpar$fr*mpar$eta)*(1+mpar$m_n/(mpar$fr*mpar$eta)*(exp(soilpar$beta*(s-mpar$s_lim))-1))^
          (-mpar$l_prime/(soilpar$beta*(mpar$fr*mpar$eta-mpar$m_n))-1)*A*A1
  }
     return(ps)
}


# .........................................................
# 6. Integration functions for new and old functions
# .........................................................

# Integrate old functions to find cc
I_old <- function(soilpar,vegpar,mpar) {
    Pd_s1<-integrate(pdf_s1,soilpar$s_h+0.001,vegpar$s_w,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
    Pd_s2<-integrate(pdf_s2,vegpar$s_w,vegpar$s_star,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
    Pd_s3<-integrate(pdf_s3,vegpar$s_star,soilpar$s_fc,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
    Pd_s4<-integrate(pdf_s4,soilpar$s_fc,1,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
    cc<-1/(Pd_s1$value+Pd_s2$value+Pd_s3$value+Pd_s4$value)
    return(cc)
}

# rho_new_1
I_new <- function(soilpar,vegpar,mpar) {
    # Integrate new functions to find cc_n
    Pd_s2_n<-integrate(pdf_s2_n,mpar$s_cr,vegpar$s_star,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
    Pd_s3_n<-integrate(function(s,soilpar,vegpar,mpar) sapply(s,pdf_s3_n,soilpar,vegpar,mpar),vegpar$s_star,mpar$s_lim,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
    Pd_s4_n<-integrate(pdf_s4_n,mpar$s_lim,1,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
    cc_n<-1/(Pd_s2_n$value+Pd_s3_n$value+Pd_s4_n$value)
    return(cc_n)
}

# rho_new_2
I_new2 <- function(soilpar,vegpar,mpar) {
      Pd_s1_n<-integrate(pdf_s1_n2,soilpar$s_h+0.001,mpar$s_cr2,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
      Pd_s2_n<-integrate(pdf_s2_n2,mpar$s_cr2,vegpar$s_star,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
      Pd_s3_n<-integrate(pdf_s3_n2,vegpar$s_star,mpar$s_lim,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
      Pd_s4_n<-integrate(pdf_s4_n2,mpar$s_lim,1,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
      cc_n<-1/(Pd_s1_n$value+Pd_s2_n$value+Pd_s3_n$value+Pd_s4_n$value)
      return(cc_n)
}

# roots
# RF model
I_root <- function(soilpar,vegpar,mpar,Z) {
      Pd_s1_r<-integrate(pdf_s1_r,mpar$s_cr.r2+0.0001,mpar$s_cr.r,soilpar=soilpar,vegpar=vegpar,mpar=mpar,Z=Z)
      Pd_s2_r<-integrate(pdf_s2_r,mpar$s_cr.r,vegpar$s_star,soilpar=soilpar,vegpar=vegpar,mpar=mpar,Z=Z)
      Pd_s3_r<-integrate(pdf_s3_r,vegpar$s_star,mpar$s_lim,soilpar=soilpar,vegpar=vegpar,mpar=mpar,Z=Z)
      Pd_s4_r<-integrate(pdf_s4_r,mpar$s_lim,1,soilpar=soilpar,vegpar=vegpar,mpar=mpar,Z=Z)
      cc_n<-1/(Pd_s1_r$value+Pd_s2_r$value+Pd_s3_r$value+Pd_s4_r$value)
      return(cc_n)
}
# feedback model
I_root2 <- function(soilpar,vegpar,mpar) {
      Pd_s2_fr<-integrate(pdf_s2_fr,mpar$s_cr.r2+0.0001,vegpar$s_star,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
      Pd_s3_fr<-integrate(pdf_s3_fr,vegpar$s_star,mpar$s_lim,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
      Pd_s4_fr<-integrate(pdf_s4_fr,mpar$s_lim,1,soilpar=soilpar,vegpar=vegpar,mpar=mpar)
      cc_n<-1/(Pd_s2_fr$value+Pd_s3_fr$value+Pd_s4_fr$value)
      return(cc_n)
}


# Inserted 20080830 from RIGWup_plotpdfs.r
# Integrate and define pdfs
# Integrate and define pdfs
pdfs <- function(soilpar,vegpar,mpar)
{
      cc<-I_old(soilpar,vegpar,mpar)
      s1<-seq(soilpar$s_h,vegpar$s_w,length=50)
      s2<-seq(vegpar$s_w,vegpar$s_star,length=50)
      s3<-seq(vegpar$s_star,soilpar$s_fc,length=50)
      s4<-seq(soilpar$s_fc,1,length=50)

      p1 <- sapply(s1,pdf_s1,soilpar,vegpar,mpar)
      p2 <- sapply(s2,pdf_s2,soilpar,vegpar,mpar)
      p3 <- sapply(s3,pdf_s3,soilpar,vegpar,mpar)
      p4 <- sapply(s4,pdf_s4,soilpar,vegpar,mpar)

      return(cbind(c(s1,s2,s3,s4),cc*c(p1,p2,p3,p4)))
}

     # new functions
pdfs_n <- function(soilpar,vegpar,mpar) {
      cc_n <- do.call(I_new,list(soilpar=soilpar,vegpar=vegpar,mpar=mpar))

      s1<-seq(soilpar$s_h,mpar$s_cr,length=50)
      s2<-seq(mpar$s_cr,vegpar$s_star,length=50)
      s3<-seq(vegpar$s_star,mpar$s_lim,length=50)
      s4<-seq(mpar$s_lim,1,length=50)

      #p1 <- pdf_s1_n2(s1,soilpar,vegpar,mpar)
      p1 <- rep(0,length(s1))
      p2 <- sapply(s2,pdf_s2_n,soilpar,vegpar,mpar)
      p3 <- sapply(s3,pdf_s3_n,soilpar,vegpar,mpar)
      p4 <- sapply(s4,pdf_s4_n,soilpar,vegpar,mpar)

      s_n <- c(s1,s2,s3,s4)

      return(cbind(s_n,cc_n*c(p1,p2,p3,p4)))
}
# roots
# This is the rootfraction model
pdfs_r <- function(soilpar,vegpar,mpar,Z) {
      cc_n <- do.call(I_root,list(soilpar=soilpar,vegpar=vegpar,mpar=mpar,Z=Z))

      s0 <- seq(soilpar$s_h,mpar$s_cr.r2,length=25)
      s1<-seq(mpar$s_cr.r2,mpar$s_cr.r,length=25)
      s2<-seq(mpar$s_cr.r,vegpar$s_star,length=50)
      s3<-seq(vegpar$s_star,mpar$s_lim,length=50)
      s4<-seq(mpar$s_lim,1,length=50)

      #p1 <- pdf_s1_n2(s1,soilpar,vegpar,mpar)
      p0 <- rep(0,25)
      p1 <- sapply(s1,pdf_s1_r,soilpar,vegpar,mpar,Z)
      p2 <- sapply(s2,pdf_s2_r,soilpar,vegpar,mpar,Z)
      p3 <- sapply(s3,pdf_s3_r,soilpar,vegpar,mpar,Z)
      p4 <- sapply(s4,pdf_s4_r,soilpar,vegpar,mpar,Z)

      s_n <- c(s0,s1,s2,s3,s4)

      return(cbind(s_n,cc_n*c(p0,p1,p2,p3,p4)))
}

# this is the feedback model
pdfs_r2 <- function(soilpar,vegpar,mpar) {
      cc_n <- do.call(I_root2,list(soilpar=soilpar,vegpar=vegpar,mpar=mpar))

      s1<-seq(soilpar$s_h,mpar$s_cr.r2,length=50)
      s2<-seq(mpar$s_cr.r2,vegpar$s_star,length=50)
      s3<-seq(vegpar$s_star,mpar$s_lim,length=50)
      s4<-seq(mpar$s_lim,1,length=50)

      #p1 <- pdf_s1_n2(s1,soilpar,vegpar,mpar)
      p1 <- rep(0,50)
      p2 <- sapply(s2,pdf_s2_fr,soilpar,vegpar,mpar)
      p3 <- sapply(s3,pdf_s3_fr,soilpar,vegpar,mpar)
      p4 <- sapply(s4,pdf_s4_fr,soilpar,vegpar,mpar)

      s_n <- c(s1,s2,s3,s4)

      return(cbind(s_n,cc_n*c(p1,p2,p3,p4)))
}

# ........................................................
# 7. calculate zeta (equation 9 & 10, Porporato et al.)
# water stress functions need to be redefined
# ..........................................................
# calculate zeta (equation 9 & 10, Porporato et al.)
zeta<-function(s,vegpar) {((vegpar$s_star-s)/(vegpar$s_star-vegpar$s_w))^q}
zeta_n<-function(s,vegpar,mpar) {((vegpar$s_star-s)/(vegpar$s_star-mpar$s_cr))^q}

#........................................
# Water stress probability density function
# equation 13 Porporato et al., 2001
#........................................

pdf_z<-function(s,vegpar,mpar) {
   ps<-1/mpar$eta_w*exp(mpar$gam_p*((vegpar$s_star-vegpar$s_w)*zeta(s,vegpar)^(1/q)-
   vegpar$s_star))*(((1-mpar$eta/mpar$eta_w)*zeta(s,vegpar)^(1/q)+(mpar$eta/mpar$eta_w))^
   ((mpar$l_prime*(vegpar$s_star-vegpar$s_w)/mpar$eta-mpar$eta_w)-1))
   return(ps)
}

pdf_z_n<-function(s,vegpar,mpar)
{
   if (mpar$eta-mpar$m2<=0) {
     ps <- rep(0,length(s))
   } else {
     B <- (mpar$l_prime*(vegpar$s_star-mpar$s_cr)/(mpar$eta-mpar$m2)-1)
#     ps<-exp(-gam_p*s)*C/eta_w*(1+(eta-m2)/eta_w*(s-s_w)/(s_star-s_w))^(B-1)
     ps<-exp(-mpar$gam_p*(vegpar$s_star-(zeta_n(s,vegpar,mpar))^(1/q)*(vegpar$s_star-mpar$s_cr)))*
     1/(mpar$eta-mpar$m2)*(((vegpar$s_star-(zeta_n(s,vegpar,mpar))^(1/q)*
     (vegpar$s_star-mpar$s_cr))-mpar$s_cr)/(vegpar$s_star-mpar$s_cr))^B
   }
     return(ps)
}

zeta_bar<-function(s,vegpar) {zeta(s,vegpar)*Cc_z*do.call(pdf_z,list(s=s,vegpar=vegpar,mpar=mpar))}
zeta_bar_n<-function(s,vegpar,mpar) {zeta(s,vegpar)*Cc_z_n*do.call(pdf_z_n,list(s=s,vegpar=vegpar,mpar=mpar))}


# Mean water stress (equation 16 Porporato et al., 2001)
temp_f<-function(s) {
     temp1<-vector(length=length(s))
     for (i in 1:length(temp1)){
         temp1[i]<-zeta(s[i])*Cc_z*pdf_z(s[i])
     }
     return(temp1)
}

#..........................................................................
# 8. Define the function "clusgen" following Rowlingson and Diggle (1993)
# ..........................................................................
clusgen<-function(cent,npoi,sd,poly) {
# cent = collection of centers
# npoi = number of points (integer!)
# sd = standard deviation of normal distribution
# poly = polygon
  centers<-cent
  np<-floor(1+npts(cent)*runif(npoi))
  xydis<-spoints(rnorm(npoi*2,0,sd))
  points<-centers[np,] + ifelse(xydis>=0,xydis,xydis)
  list(centers=centers,points=points)
# output is a list with centers and points
}

# -----------------------------------------------------------------------
# 9. Definition of some older obsolete functions more for completeness than necessity
# -----------------------------------------------------------------------
#...............................................
# Losses: Evaporation (equation 12 Laio et al., 2001b)
#................................................
E<-function(s)
{
  et<-vector(length=length(s))
  for (i in 1:length(s)) {
     if (s[i]>s_h && s[i]<=s_w)
     {
        et[i]<-E_w*((s[i]-s_h)/(s_w-s_h))
     } else {
         if (s[i]>s_w && s[i]<=s_star)
         {
            et[i]<-E_w+(E_max-E_w)*((s[i]-s_w)/(s_star-s_w))
         } else {
            if (s[i]>s_star && s[i]<=1) {
              et[i]<-E_max
            } else {
              et[i] <- 0
            }
         }
     }
  }
  return(et)
}

#..................................................
# Losses: Leakage (equation 14 Laio et al., 2001b)
#..................................................
K<-function(s,beta=NULL,s_fc=0.1)
{
   k<-vector(length=length(s))
   for (i in 1:length(s)) {
         k[i]<-ifelse(s[i]>s_fc && s[i]<=1,K_s/(exp(beta*(1-s_fc))-1)*(exp(beta*(s[i]-s_fc))-1),0)
   }
   return(k)
}

# Also define soil water characteristic (equation 13, Laio et al., 2001a)
# not really needed though just for completeness

psi<-function(s,b,psi_bar)
{
    psi_bar*s^b
}
