setwd()
rm(list=ls())
load("dat_deident.Rda")
bob<-dat_deident
#Packagues
library(psych)
library(lavaan)
library(semTools)

#Descriptives ######################

#Subject 96 with wrong puntuation for BMI at pre-treatment
bob$bmi_0[96] 
#Subsitute for correct value. Does not influence reliability tests for EDEQ
bob$bmi_0[96] <- 16.85989

#Age and sex
summary(bob[c("age","sex")])
describe(bob$age)

#EDEQ and BMI without standarization
describe(bob[c(paste0("edeq_", c(0:2)),paste0("bmi_", c(0, 5, 10)))]) [c("n", "mean", "sd","min", "max", "skew","kurtosis")]

#Correlations
round(cor(bob[c(paste0("edeq_", c(0:2)),paste0("bmi_", c(0, 5, 10)))],
          use="pairwise.complete.obs"),3)

#Correlation tests: BMI and EDEQ in pre-ttm, 5weeks, 10weeks
cor.test(bob$bmi_0, bob$edeq_0, method="pearson", conf.level = 0.95)
cor.test(bob$bmi_5, bob$edeq_1, method="pearson", conf.level = 0.95)
cor.test(bob$bmi_10, bob$edeq_2, method="pearson", conf.level = 0.95)

#Standardize variables ######################
#BMI
bmiMn <- mean(bob$bmi_0, na.rm = TRUE) 
bmiSd <- sd(bob$bmi_0, na.rm = TRUE)

bmivars <- names(bob)[grepl("bmi", names(dat_deident))] # save old names
newbmivars <- paste0(bmivars, "_std") # new names

bob[newbmivars] <- (bob[bmivars] - bmiMn)/bmiSd # Standardization

#EDEQ
edeqMn <- mean(bob$edeq_0, na.rm = TRUE) 
edeqSd <- sd(bob$edeq_0, na.rm = TRUE)

edeqvars <- paste0("edeq_", c(0:2)) # save old names
newedeqvars <- paste0(edeqvars, "_std") # new names

bob[newedeqvars] <- (bob[edeqvars] - edeqMn)/edeqSd # Standardization

rm(bmiMn,bmiSd,bmivars,edeqMn,edeqSd,edeqvars)

#Descriptive for standardized variables
describe(bob[newbmivars])[c("n", "mean", "sd","skew")]
describe(bob[newedeqvars])[c("n", "mean", "sd","skew")]

# Re-name for the code ##################

#BMI as X
names(bob)[names(bob) == 'bmi_0_std'] <- 'ox1'
names(bob)[names(bob) == 'bmi_5_std'] <- 'ox2'
names(bob)[names(bob) == 'bmi_10_std'] <- 'ox3'
#EDEQ as Y
names(bob)[names(bob) == 'edeq_0_std'] <- 'oy1'
names(bob)[names(bob) == 'edeq_1_std'] <- 'oy2'
names(bob)[names(bob) == 'edeq_2_std'] <- 'oy3'

#Model 0/Most restrictive one ####################################
LCSM0 = "
# Declaring latent level
y1	=~ 1* oy1
y2	=~ 1* oy2
y3	=~ 1* oy3

x1	=~ 1* ox1
x2	=~ 1* ox2
x3	=~ 1* ox3

# Auto-regression
y2	~ 1* y1
y3	~ 1* y2

x2	~ 1* x1
x3	~ 1* x2

# Define latent change
dy2	 =~ 1* y2
dy3	 =~ 1* y3

dx2	 =~ 1* x2
dx3	 =~ 1* x3

# Auto-proportions
dy2	 ~ 0 * y1
dy3	 ~ 0 * y2

dx2	 ~ 0 * x1
dx3	 ~ 0 * x2

# Couplings
dy2	 ~ 0 * x1
dy3	 ~ 0 * x2

dx2	 ~ 0 * y1
dx3	 ~ 0 * y2

# Latent intercepts and slopes
yInt =~ 1 * y1
ySlp =~ 1*dy2  + 1*dy3 

xInt =~ 1 * x1
xSlp =~ 1*dx2  + 1*dx3 

# Latent means
yInt ~ yInMn * 1
ySlp ~ ySlMn * 1

xInt ~ xInMn * 1
xSlp ~ xSlMn * 1

# Observed means fixed at 0
oy1 ~ 0 * 1
oy2 ~ 0 * 1
oy3 ~ 0 * 1

ox1 ~ 0 * 1
ox2 ~ 0 * 1
ox3 ~ 0 * 1

# Latent variances and covariances
yInt ~~ yInV * yInt   #varianza interseccion Y
ySlp ~~ 0 * ySlp      #varianza pendiente Y
yInt ~~ 0 * ySlp      #covarianza int slp y

xInt ~~ xInV * xInt   #varianza intersecci0n X
xSlp ~~ 0 * xSlp      #varianza pendiente X
xInt ~~ 0 * xSlp      #covarianza int slp x

yInt ~~ xyIntCv * xInt
ySlp ~~ 0 * xSlp
yInt ~~ 0 * xSlp
xInt ~~ 0 * ySlp

#Residual variance (Errores (varianzas) de latentes)
y1  ~~ 0 * y1
y2  ~~ 0 * y2
y3  ~~ 0 * y3

x1  ~~ 0 * x1
x2  ~~ 0 * x2
x3  ~~ 0 * x3

# Dynamic errors
dy2  ~~ 0 * dy2
dy3  ~~ 0 * dy3
dx2  ~~ 0 * dx2
dx3  ~~ 0 * dx3
# Covariance between dynamic errors
dy2 ~~ 0 * dx2
dy3 ~~ 0 * dx3

#Error de medida observada (MerX, MerY)
oy1  ~~ MerY * oy1
oy2  ~~ MerY * oy2
oy3  ~~ MerY * oy3

ox1  ~~ MerX * ox1
ox2  ~~ MerX * ox2
ox3  ~~ MerX * ox3

#Covarianza de errores de medida observada (MerXY)
oy1  ~~ MerXY * ox1
oy2  ~~ MerXY * ox2
oy3  ~~ MerXY * ox3

"
# Fit model
LCSM0.fit = lavaan(LCSM0, data=bob, missing="ML")
LCSM0.sum <- summary(LCSM0.fit, fit.measures=T) # Obtain summary
LCSM0.pars <- LCSM0.sum$PE
LCSM0.pars <- LCSM0.pars[!is.na(LCSM0.pars$z),]
LCSM0.pars <- LCSM0.pars[!duplicated(round(LCSM0.pars$est,5)),]
LCSM0.pars

#Fit Measures
fitmeasures(LCSM0.fit)
#Parameters estimation
cbind(LCSM0.pars[c("lhs", "op", "rhs", "label")],
      round(LCSM0.pars[c("est", "se", "z", "pvalue")], 3))
#Model 1 #########################################
LCSM1 = "
# Declaring latent level
y1	=~ 1* oy1
y2	=~ 1* oy2
y3	=~ 1* oy3

x1	=~ 1* ox1
x2	=~ 1* ox2
x3	=~ 1* ox3

# Auto-regression
y2	~ 1* y1
y3	~ 1* y2

x2	~ 1* x1
x3	~ 1* x2

# Define latent change
dy2	 =~ 1* y2
dy3	 =~ 1* y3

dx2	 =~ 1* x2
dx3	 =~ 1* x3

#M1 Estimaci?n de par?metro de autorregresion nivel-cambio
# Auto-proportions
dy2	 ~ b_y * y1
dy3	 ~ b_y * y2

dx2	 ~ b_x * x1
dx3	 ~ b_x * x2

# Couplings
dy2	 ~ 0 * x1
dy3	 ~ 0 * x2

dx2	 ~ 0 * y1
dx3	 ~ 0 * y2

# Latent intercepts and slopes
yInt =~ 1 * y1
ySlp =~ 1*dy2  + 1*dy3 

xInt =~ 1 * x1
xSlp =~ 1*dx2  + 1*dx3 

# Latent means
yInt ~ yInMn * 1
ySlp ~ ySlMn * 1

xInt ~ xInMn * 1
xSlp ~ xSlMn * 1

# Observed means fixed at 0
oy1 ~ 0 * 1
oy2 ~ 0 * 1
oy3 ~ 0 * 1

ox1 ~ 0 * 1
ox2 ~ 0 * 1
ox3 ~ 0 * 1

# Latent variances and covariances
yInt ~~ yInV * yInt #varianza intersecci?n Y
ySlp ~~ 0 * ySlp #varianza pendiente Y
yInt ~~ 0 * ySlp #covarianza int slp y

xInt ~~ xInV * xInt #varianza intersecci?n X
xSlp ~~ 0 * xSlp #varianza pendiente X
xInt ~~ 0 * xSlp #covarianza int slp x

yInt ~~ xyIntCv * xInt
ySlp ~~ 0 * xSlp
yInt ~~ 0 * xSlp
xInt ~~ 0 * ySlp

#Errores (varianzas) de latentes
y1  ~~ 0 * y1
y2  ~~ 0 * y2
y3  ~~ 0 * y3

x1  ~~ 0 * x1
x2  ~~ 0 * x2
x3  ~~ 0 * x3

# Dynamic errors
dy2  ~~ 0 * dy2
dy3  ~~ 0 * dy3
dx2  ~~ 0 * dx2
dx3  ~~ 0 * dx3
# Covariance between dynamic errors
dy2 ~~ 0 * dx2
dy3 ~~ 0 * dx3

#Error de medida observada (MerX, MerY)
oy1  ~~ MerY * oy1
oy2  ~~ MerY * oy2
oy3  ~~ MerY * oy3

ox1  ~~ MerX * ox1
ox2  ~~ MerX * ox2
ox3  ~~ MerX * ox3

#Covarianza de errores de medida observada (MerXY)
oy1  ~~ MerXY * ox1
oy2  ~~ MerXY * ox2
oy3  ~~ MerXY * ox3

"
# Fit model
LCSM1.fit = lavaan(LCSM1, data=bob, missing="ML")
LCSM1.sum <- summary(LCSM1.fit, fit.measures=T) # Obtain summary
LCSM1.pars <- LCSM1.sum$PE
LCSM1.pars <- LCSM1.pars[!is.na(LCSM1.pars$z),]
LCSM1.pars <- LCSM1.pars[!duplicated(round(LCSM1.pars$est,5)),]
LCSM1.pars

#Fit Measures
fitmeasures(LCSM1.fit)
#Parameters estimation
cbind(LCSM1.pars[c("lhs", "op", "rhs", "label")],
      round(LCSM1.pars[c("est", "se", "z", "pvalue")], 3))

#Model 2 #########################################
LCSM2 = "
# Declaring latent level
y1	=~ 1* oy1
y2	=~ 1* oy2
y3	=~ 1* oy3

x1	=~ 1* ox1
x2	=~ 1* ox2
x3	=~ 1* ox3

# Auto-regression
y2	~ 1* y1
y3	~ 1* y2

x2	~ 1* x1
x3	~ 1* x2

# Define latent change
dy2	 =~ 1* y2
dy3	 =~ 1* y3

dx2	 =~ 1* x2
dx3	 =~ 1* x3

# Auto-proportions
dy2	 ~ b_y * y1
dy3	 ~ b_y * y2

dx2	 ~ b_x * x1
dx3	 ~ b_x * x2

#M2 Estimaci?n de par?metro de emparejamiento nivel-cambio
# Couplings
dy2	 ~ g_y * x1
dy3	 ~ g_y * x2

dx2	 ~ g_x * y1
dx3	 ~ g_x * y2

# Latent intercepts and slopes
yInt =~ 1 * y1
ySlp =~ 1*dy2  + 1*dy3 

xInt =~ 1 * x1
xSlp =~ 1*dx2  + 1*dx3 

# Latent means
yInt ~ yInMn * 1
ySlp ~ ySlMn * 1

xInt ~ xInMn * 1
xSlp ~ xSlMn * 1

# Observed means fixed at 0
oy1 ~ 0 * 1
oy2 ~ 0 * 1
oy3 ~ 0 * 1

ox1 ~ 0 * 1
ox2 ~ 0 * 1
ox3 ~ 0 * 1

# Latent variances and covariances
yInt ~~ yInV * yInt #varianza intersecci?n Y
ySlp ~~ 0 * ySlp #varianza pendiente Y
yInt ~~ 0 * ySlp #covarianza int slp y

xInt ~~ xInV * xInt #varianza intersecci?n X
xSlp ~~ 0 * xSlp #varianza pendiente X
xInt ~~ 0 * xSlp #covarianza int slp x

yInt ~~ xyIntCv * xInt
ySlp ~~ 0 * xSlp
yInt ~~ 0 * xSlp
xInt ~~ 0 * ySlp

#Errores (varianzas) de latentes
y1  ~~ 0 * y1
y2  ~~ 0 * y2
y3  ~~ 0 * y3

x1  ~~ 0 * x1
x2  ~~ 0 * x2
x3  ~~ 0 * x3

# Dynamic errors
dy2  ~~ 0 * dy2
dy3  ~~ 0 * dy3
dx2  ~~ 0 * dx2
dx3  ~~ 0 * dx3
# Covariance between dynamic errors
dy2 ~~ 0 * dx2
dy3 ~~ 0 * dx3

#Error de medida observada (MerX, MerY)
oy1  ~~ MerY * oy1
oy2  ~~ MerY * oy2
oy3  ~~ MerY * oy3

ox1  ~~ MerX * ox1
ox2  ~~ MerX * ox2
ox3  ~~ MerX * ox3

#Covarianza de errores de medida observada (MerXY)
oy1  ~~ MerXY * ox1
oy2  ~~ MerXY * ox2
oy3  ~~ MerXY * ox3

"
# Fit model
LCSM2.fit = lavaan(LCSM2, data=bob, missing="ML")
LCSM2.sum <- summary(LCSM2.fit, fit.measures=T) # Obtain summary
LCSM2.pars <- LCSM2.sum$PE
LCSM2.pars <- LCSM2.pars[!is.na(LCSM2.pars$z),]
LCSM2.pars <- LCSM2.pars[!duplicated(round(LCSM2.pars$est,5)),]
LCSM2.pars

#Fit Measures
fitmeasures(LCSM2.fit)
#Parameters estimation
cbind(LCSM2.pars[c("lhs", "op", "rhs", "label")],
      round(LCSM2.pars[c("est", "se", "z", "pvalue")], 3))

#Model 3 #########################################
LCSM3 = "
# Declaring latent level
y1	=~ 1* oy1
y2	=~ 1* oy2
y3	=~ 1* oy3

x1	=~ 1* ox1
x2	=~ 1* ox2
x3	=~ 1* ox3

# Auto-regression
y2	~ 1* y1
y3	~ 1* y2

x2	~ 1* x1
x3	~ 1* x2

# Define latent change
dy2	 =~ 1* y2
dy3	 =~ 1* y3

dx2	 =~ 1* x2
dx3	 =~ 1* x3

# Auto-proportions
dy2	 ~ b_y * y1
dy3	 ~ b_y * y2

dx2	 ~ b_x * x1
dx3	 ~ b_x * x2

# Couplings
dy2	 ~ g_y * x1
dy3	 ~ g_y * x2

dx2	 ~ g_x * y1
dx3	 ~ g_x * y2

#M3 estimaci?n de par?metro de autorregresion del cambio
#Phi (autorregresion del cambio) 
dy3	 ~ p_y * dy2

dx3	 ~ p_x * dx2

# Latent intercepts and slopes
yInt =~ 1 * y1
ySlp =~ 1*dy2  + 1*dy3 

xInt =~ 1 * x1
xSlp =~ 1*dx2  + 1*dx3 

# Latent means
yInt ~ yInMn * 1
ySlp ~ ySlMn * 1

xInt ~ xInMn * 1
xSlp ~ xSlMn * 1

# Observed means fixed at 0
oy1 ~ 0 * 1
oy2 ~ 0 * 1
oy3 ~ 0 * 1

ox1 ~ 0 * 1
ox2 ~ 0 * 1
ox3 ~ 0 * 1

# Latent variances and covariances
yInt ~~ yInV * yInt #varianza intersecci?n Y
ySlp ~~ 0 * ySlp #varianza pendiente Y
yInt ~~ 0 * ySlp #covarianza int slp y

xInt ~~ xInV * xInt #varianza intersecci?n X
xSlp ~~ 0 * xSlp #varianza pendiente X
xInt ~~ 0 * xSlp #covarianza int slp x

yInt ~~ xyIntCv * xInt
ySlp ~~ 0 * xSlp
yInt ~~ 0 * xSlp
xInt ~~ 0 * ySlp

#Errores (varianzas) de latentes
y1  ~~ 0 * y1
y2  ~~ 0 * y2
y3  ~~ 0 * y3

x1  ~~ 0 * x1
x2  ~~ 0 * x2
x3  ~~ 0 * x3

# Dynamic errors
dy2  ~~ 0 * dy2
dy3  ~~ 0 * dy3
dx2  ~~ 0 * dx2
dx3  ~~ 0 * dx3
# Covariance between dynamic errors
dy2 ~~ 0 * dx2
dy3 ~~ 0 * dx3

#Error de medida observada (MerX, MerY)
oy1  ~~ MerY * oy1
oy2  ~~ MerY * oy2
oy3  ~~ MerY * oy3

ox1  ~~ MerX * ox1
ox2  ~~ MerX * ox2
ox3  ~~ MerX * ox3

#Covarianza de errores de medida observada (MerXY)
oy1  ~~ MerXY * ox1
oy2  ~~ MerXY * ox2
oy3  ~~ MerXY * ox3

"
# Fit model
LCSM3.fit = lavaan(LCSM3, data=bob, missing="ML")
LCSM3.sum <- summary(LCSM3.fit, fit.measures=T) # Obtain summary
LCSM3.pars <- LCSM3.sum$PE
LCSM3.pars <- LCSM3.pars[!is.na(LCSM3.pars$z),]
LCSM3.pars <- LCSM3.pars[!duplicated(round(LCSM3.pars$est,5)),]
LCSM3.pars

#Fit Measures
fitmeasures(LCSM3.fit)
#Parameters estimation
cbind(LCSM3.pars[c("lhs", "op", "rhs", "label")],
      round(LCSM3.pars[c("est", "se", "z", "pvalue")], 3))

#Model 4 #########################################
LCSM4 = "
# Declaring latent level
y1	=~ 1* oy1
y2	=~ 1* oy2
y3	=~ 1* oy3

x1	=~ 1* ox1
x2	=~ 1* ox2
x3	=~ 1* ox3

# Auto-regression
y2	~ 1* y1
y3	~ 1* y2

x2	~ 1* x1
x3	~ 1* x2

# Define latent change
dy2	 =~ 1* y2
dy3	 =~ 1* y3

dx2	 =~ 1* x2
dx3	 =~ 1* x3

# Auto-proportions
dy2	 ~ b_y * y1
dy3	 ~ b_y * y2

dx2	 ~ b_x * x1
dx3	 ~ b_x * x2

# Couplings
dy2	 ~ g_y * x1
dy3	 ~ g_y * x2

dx2	 ~ g_x * y1
dx3	 ~ g_x * y2

#Phi (autorregresion del cambio) 
dy3	 ~ p_y * dy2

dx3	 ~ p_x * dx2

#M4 estimaci?n de par?metro de emperejamiento cambio-cambio
#Ksi (coupling del cambio)
dy3	 ~ k_y * dx2

dx3	 ~ k_x * dy2

# Latent intercepts and slopes
yInt =~ 1 * y1
ySlp =~ 1*dy2  + 1*dy3 

xInt =~ 1 * x1
xSlp =~ 1*dx2  + 1*dx3 

# Latent means
yInt ~ yInMn * 1
ySlp ~ ySlMn * 1

xInt ~ xInMn * 1
xSlp ~ xSlMn * 1

# Observed means fixed at 0
oy1 ~ 0 * 1
oy2 ~ 0 * 1
oy3 ~ 0 * 1

ox1 ~ 0 * 1
ox2 ~ 0 * 1
ox3 ~ 0 * 1

# Latent variances and covariances
yInt ~~ yInV * yInt #varianza intersecci?n Y
ySlp ~~ 0 * ySlp #varianza pendiente Y
yInt ~~ 0 * ySlp #covarianza int slp y

xInt ~~ xInV * xInt #varianza intersecci?n X
xSlp ~~ 0 * xSlp #varianza pendiente X
xInt ~~ 0 * xSlp #covarianza int slp x

yInt ~~ xyIntCv * xInt
ySlp ~~ 0 * xSlp
yInt ~~ 0 * xSlp
xInt ~~ 0 * ySlp

#Errores (varianzas) de latentes
y1  ~~ 0 * y1
y2  ~~ 0 * y2
y3  ~~ 0 * y3

x1  ~~ 0 * x1
x2  ~~ 0 * x2
x3  ~~ 0 * x3

# Dynamic errors
dy2  ~~ 0 * dy2
dy3  ~~ 0 * dy3
dx2  ~~ 0 * dx2
dx3  ~~ 0 * dx3
# Covariance between dynamic errors
dy2 ~~ 0 * dx2
dy3 ~~ 0 * dx3

#Error de medida observada (MerX, MerY)
oy1  ~~ MerY * oy1
oy2  ~~ MerY * oy2
oy3  ~~ MerY * oy3

ox1  ~~ MerX * ox1
ox2  ~~ MerX * ox2
ox3  ~~ MerX * ox3

#Covarianza de errores de medida observada (MerXY)
oy1  ~~ MerXY * ox1
oy2  ~~ MerXY * ox2
oy3  ~~ MerXY * ox3

"
# Fit model
LCSM4.fit = lavaan(LCSM4, data=bob, missing="ML")
LCSM4.sum <- summary(LCSM4.fit, fit.measures=T) # Obtain summary
LCSM4.pars <- LCSM4.sum$PE
LCSM4.pars <- LCSM4.pars[!is.na(LCSM4.pars$z),]
LCSM4.pars <- LCSM4.pars[!duplicated(round(LCSM4.pars$est,5)),]
LCSM4.pars

#Fit Measures
fitmeasures(LCSM4.fit)
#Parameters estimation
cbind(LCSM4.pars[c("lhs", "op", "rhs", "label")],
      round(LCSM4.pars[c("est", "se", "z", "pvalue")], 3))

#Model 5 #########################################
LCSM5 = "
# Declaring latent level
y1	=~ 1* oy1
y2	=~ 1* oy2
y3	=~ 1* oy3

x1	=~ 1* ox1
x2	=~ 1* ox2
x3	=~ 1* ox3

# Auto-regression
y2	~ 1* y1
y3	~ 1* y2

x2	~ 1* x1
x3	~ 1* x2

# Define latent change
dy2	 =~ 1* y2
dy3	 =~ 1* y3

dx2	 =~ 1* x2
dx3	 =~ 1* x3

# Auto-proportions
dy2	 ~ b_y * y1
dy3	 ~ b_y * y2

dx2	 ~ b_x * x1
dx3	 ~ b_x * x2

# Couplings
dy2	 ~ g_y * x1
dy3	 ~ g_y * x2

dx2	 ~ g_x * y1
dx3	 ~ g_x * y2

#Phi (autorregresion del cambio) 
dy3	 ~ p_y * dy2

dx3	 ~ p_x * dx2

#Ksi (coupling del cambio)
dy3	 ~ k_y * dx2

dx3	 ~ k_x * dy2

# Latent intercepts and slopes
yInt =~ 1 * y1
ySlp =~ 1*dy2  + 1*dy3 

xInt =~ 1 * x1
xSlp =~ 1*dx2  + 1*dx3 

# Latent means
yInt ~ yInMn * 1
ySlp ~ ySlMn * 1

xInt ~ xInMn * 1
xSlp ~ xSlMn * 1

# Observed means fixed at 0
oy1 ~ 0 * 1
oy2 ~ 0 * 1
oy3 ~ 0 * 1

ox1 ~ 0 * 1
ox2 ~ 0 * 1
ox3 ~ 0 * 1

#M5 liberaci?n de varianzas y covarianzas xy del componente aditivo
# Latent variances and covariances
yInt ~~ yInV   * yInt #varianza intersecci?n Y
ySlp ~~ ySlV   * ySlp #varianza pendiente Y
yInt ~~ 0* ySlp #covarianza int slp y

xInt ~~ xInV   * xInt #varianza intersecci?n X
xSlp ~~ xSlpV * xSlp #varianza pendiente X
xInt ~~ 0 * xSlp #covarianza int slp x

yInt ~~ xyIntCv * xInt
ySlp ~~ xySlpCv * xSlp
yInt ~~ 0 * xSlp
xInt ~~ 0 * ySlp

#Errores (varianzas) de latentes
y1  ~~ 0 * y1
y2  ~~ 0 * y2
y3  ~~ 0 * y3

x1  ~~ 0 * x1
x2  ~~ 0 * x2
x3  ~~ 0 * x3

# Dynamic errors
dy2  ~~ 0 * dy2
dy3  ~~ 0 * dy3
dx2  ~~ 0 * dx2
dx3  ~~ 0 * dx3
# Covariance between dynamic errors
dy2 ~~ 0 * dx2
dy3 ~~ 0 * dx3

#Error de medida observada (MerX, MerY)
oy1  ~~ MerY * oy1
oy2  ~~ MerY * oy2
oy3  ~~ MerY * oy3

ox1  ~~ MerX * ox1
ox2  ~~ MerX * ox2
ox3  ~~ MerX * ox3

#Covarianza de errores de medida observada (MerXY)
oy1  ~~ MerXY * ox1
oy2  ~~ MerXY * ox2
oy3  ~~ MerXY * ox3

"
# Fit model
LCSM5.fit = lavaan(LCSM5, data=bob, missing="ML")
LCSM5.sum <- summary(LCSM5.fit, fit.measures=T) # Obtain summary
LCSM5.pars <- LCSM5.sum$PE
LCSM5.pars <- LCSM5.pars[!is.na(LCSM5.pars$z),]
LCSM5.pars <- LCSM5.pars[!duplicated(round(LCSM5.pars$est,5)),]
LCSM5.pars

#Fit Measures
fitmeasures(LCSM5.fit)
#Parameters estimation
cbind(LCSM5.pars[c("lhs", "op", "rhs", "label")],
      round(LCSM5.pars[c("est", "se", "z", "pvalue")], 3))

#Model 6 #########################################
LCSM6 = "
# Declaring latent level
y1	=~ 1* oy1
y2	=~ 1* oy2
y3	=~ 1* oy3

x1	=~ 1* ox1
x2	=~ 1* ox2
x3	=~ 1* ox3

# Auto-regression
y2	~ 1* y1
y3	~ 1* y2

x2	~ 1* x1
x3	~ 1* x2

# Define latent change
dy2	 =~ 1* y2
dy3	 =~ 1* y3

dx2	 =~ 1* x2
dx3	 =~ 1* x3

# Auto-proportions
dy2	 ~ b_y * y1
dy3	 ~ b_y * y2

dx2	 ~ b_x * x1
dx3	 ~ b_x * x2

# Couplings
dy2	 ~ g_y * x1
dy3	 ~ g_y * x2

dx2	 ~ g_x * y1
dx3	 ~ g_x * y2

#Phi (autorregresion del cambio) 
dy3	 ~ p_y * dy2

dx3	 ~ p_x * dx2

#Ksi (coupling del cambio)
dy3	 ~ k_y * dx2

dx3	 ~ k_x * dy2

# Latent intercepts and slopes
yInt =~ 1 * y1
ySlp =~ 1*dy2  + 1*dy3 

xInt =~ 1 * x1
xSlp =~ 1*dx2  + 1*dx3 

# Latent means
yInt ~ yInMn * 1
ySlp ~ ySlMn * 1

xInt ~ xInMn * 1
xSlp ~ xSlMn * 1

# Observed means fixed at 0
oy1 ~ 0 * 1
oy2 ~ 0 * 1
oy3 ~ 0 * 1

ox1 ~ 0 * 1
ox2 ~ 0 * 1
ox3 ~ 0 * 1

# Latent variances and covariances
yInt ~~ yInV   * yInt #varianza intersecci?n Y
ySlp ~~ ySlV   * ySlp #varianza pendiente Y
yInt ~~ yInSlCv* ySlp #covarianza int slp y

xInt ~~ xInV   * xInt #varianza intersecci?n X
xSlp ~~ xSlpV * xSlp #varianza pendiente X
xInt ~~ xInSlCv * xSlp #covarianza int slp x

#M6 Liberaci?n de covarianzas:
yInt ~~ yxInCv   * xInt
ySlp ~~ xySlCv * xSlp
yInt ~~ yIntxSlCv * xSlp
xInt ~~ xInySlCv   * ySlp

#Errores (varianzas) de latentes
y1  ~~ 0 * y1
y2  ~~ 0 * y2
y3  ~~ 0 * y3

x1  ~~ 0 * x1
x2  ~~ 0 * x2
x3  ~~ 0 * x3

# Dynamic errors
dy2  ~~ 0 * dy2
dy3  ~~ 0 * dy3
dx2  ~~ 0 * dx2
dx3  ~~ 0 * dx3
# Covariance between dynamic errors
dy2 ~~ 0 * dx2
dy3 ~~ 0 * dx3

#Error de medida observada (MerX, MerY)
oy1  ~~ MerY * oy1
oy2  ~~ MerY * oy2
oy3  ~~ MerY * oy3

ox1  ~~ MerX * ox1
ox2  ~~ MerX * ox2
ox3  ~~ MerX * ox3

#Covarianza de errores de medida observada (MerXY)
oy1  ~~ MerXY * ox1
oy2  ~~ MerXY * ox2
oy3  ~~ MerXY * ox3

"
# Fit model
LCSM6.fit = lavaan(LCSM6, data=bob, missing="ML")
#  lavaan WARNING:
#  The variance-covariance matrix of the estimated parameters (vcov)
#  does not appear to be positive definite! The smallest eigenvalue
#  (= 3.485592e-18) is close to zero. This may be a symptom that the
#  model is not identified.
LCSM6.sum <- summary(LCSM6.fit, fit.measures=T) # Obtain summary
LCSM6.pars <- LCSM6.sum$PE
LCSM6.pars <- LCSM6.pars[!is.na(LCSM6.pars$z),]
LCSM6.pars <- LCSM6.pars[!duplicated(round(LCSM6.pars$est,5)),]
LCSM6.pars

#Fit Measures
fitmeasures(LCSM6.fit)
#Parameters estimation
cbind(LCSM6.pars[c("lhs", "op", "rhs", "label")],
      round(LCSM6.pars[c("est", "se", "z", "pvalue")], 3))

#Grouping of results. Model comparisons #########################

modelL <- list()
modelL$m0 <- LCSM0
modelL$m1 <- LCSM1
modelL$m2 <- LCSM2
modelL$m3 <- LCSM3
modelL$m4 <- LCSM4
modelL$m5 <- LCSM5
modelL$m6 <- LCSM6

fitL <- list()
fitL$m0 <- lavaan(data=bob, modelL$m0, missing="ML")
fitL$m1 <- lavaan(data=bob, modelL$m1, missing="ML")
fitL$m2 <- lavaan(data=bob, modelL$m2, missing="ML")
fitL$m3 <- lavaan(data=bob, modelL$m3, missing="ML")
fitL$m4 <- lavaan(data=bob, modelL$m4, missing="ML")
fitL$m5 <- lavaan(data=bob, modelL$m5, missing="ML")
fitL$m6 <- lavaan(data=bob, modelL$m6, missing="ML")


fitms <- c("npar","logl","aic",
           "chisq", "df", "pvalue",
           "rmsea", "cfi", "srmr")
fitmeasL <- list()
fitmeasL$m0 <- fitmeasures(fitL$m0, fit.measures = fitms)
fitmeasL$m1 <- fitmeasures(fitL$m1, fit.measures = fitms)
fitmeasL$m2 <- fitmeasures(fitL$m2, fit.measures = fitms)
fitmeasL$m3 <- fitmeasures(fitL$m3, fit.measures = fitms)
fitmeasL$m4 <- fitmeasures(fitL$m4, fit.measures = fitms)
fitmeasL$m5 <- fitmeasures(fitL$m5, fit.measures = fitms)
fitmeasL$m6 <- fitmeasures(fitL$m6, fit.measures = fitms)

summodels <- sapply(fitmeasL, function(f) {round(f,3)})
summodels

compmodels <- compareFit(fitL$m0, fitL$m1, fitL$m2, fitL$m3,
                         fitL$m4, fitL$m5,
                         nested=TRUE) 
summary(compmodels)

#Best model is 5. Coupling paramenters do not fit better
#Try new model 5 B. Coupling parameters (g_ and k_) constrained to zero

#Model 5 Optimized (M5B) #########################
LCSM5B = "
# Declaring latent level
y1	=~ 1* oy1
y2	=~ 1* oy2
y3	=~ 1* oy3

x1	=~ 1* ox1
x2	=~ 1* ox2
x3	=~ 1* ox3

# Auto-regression
y2	~ 1* y1
y3	~ 1* y2

x2	~ 1* x1
x3	~ 1* x2

# Define latent change
dy2	 =~ 1* y2
dy3	 =~ 1* y3

dx2	 =~ 1* x2
dx3	 =~ 1* x3

# Auto-proportions
dy2	 ~ b_y * y1
dy3	 ~ b_y * y2

dx2	 ~ b_x * x1
dx3	 ~ b_x * x2

# Couplings RESTRINGIDOS
dy2	 ~ 0 * x1
dy3	 ~ 0 * x2

dx2	 ~ 0 * y1
dx3	 ~ 0 * y2

#Phi (autorregresion del cambio) 
dy3	 ~ p_y * dy2

dx3	 ~ p_x * dx2

#RESTRICCION
#Ksi (coupling del cambio)
dy3	 ~ 0 * dx2

dx3	 ~ 0 * dy2

# Latent intercepts and slopes
yInt =~ 1 * y1
ySlp =~ 1*dy2  + 1*dy3 

xInt =~ 1 * x1
xSlp =~ 1*dx2  + 1*dx3 

# Latent means
yInt ~ yInMn * 1
ySlp ~ ySlMn * 1

xInt ~ xInMn * 1
xSlp ~ xSlMn * 1

# Observed means fixed at 0
oy1 ~ 0 * 1
oy2 ~ 0 * 1
oy3 ~ 0 * 1

ox1 ~ 0 * 1
ox2 ~ 0 * 1
ox3 ~ 0 * 1

#M5 liberaci?n de varianzas y covarianzas xy del componente aditivo
# Latent variances and covariances
yInt ~~ yInV   * yInt #varianza intersecci?n Y
ySlp ~~ ySlV   * ySlp #varianza pendiente Y
yInt ~~ 0* ySlp #covarianza int slp y

xInt ~~ xInV   * xInt #varianza intersecci?n X
xSlp ~~ xSlpV * xSlp #varianza pendiente X
xInt ~~ 0 * xSlp #covarianza int slp x

yInt ~~ xyIntCv * xInt
ySlp ~~ xySlpCv * xSlp
yInt ~~ 0 * xSlp
xInt ~~ 0 * ySlp

#Errores (varianzas) de latentes
y1  ~~ 0 * y1
y2  ~~ 0 * y2
y3  ~~ 0 * y3

x1  ~~ 0 * x1
x2  ~~ 0 * x2
x3  ~~ 0 * x3

# Dynamic errors
dy2  ~~ 0 * dy2
dy3  ~~ 0 * dy3
dx2  ~~ 0 * dx2
dx3  ~~ 0 * dx3
# Covariance between dynamic errors
dy2 ~~ 0 * dx2
dy3 ~~ 0 * dx3

#Error de medida observada (MerX, MerY)
oy1  ~~ MerY * oy1
oy2  ~~ MerY * oy2
oy3  ~~ MerY * oy3

ox1  ~~ MerX * ox1
ox2  ~~ MerX * ox2
ox3  ~~ MerX * ox3

#Covarianza de errores de medida observada (MerXY)
oy1  ~~ MerXY * ox1
oy2  ~~ MerXY * ox2
oy3  ~~ MerXY * ox3

"
# Fit model
LCSM5B.fit = lavaan(LCSM5B, data=bob, missing="ML")
LCSM5B.sum <- summary(LCSM5B.fit, fit.measures=T) # Obtain summary
LCSM5B.pars <- LCSM5B.sum$PE
LCSM5B.pars <- LCSM5B.pars[!is.na(LCSM5B.pars$z),]
LCSM5B.pars <- LCSM5B.pars[!duplicated(round(LCSM5B.pars$est,5)),]
LCSM5B.pars

#Fit Measures
fitmeasures(LCSM5B.fit)
#Parameters estimation
cbind(LCSM5B.pars[c("lhs", "op", "rhs", "label")],
      round(LCSM5B.pars[c("est", "se", "z", "pvalue")], 3))

estM5B <- cbind(LCSM5B.pars[c("lhs", "op", "rhs", "label")],
                       round(LCSM5B.pars[c("est", "se", "z", "pvalue")], 3))
#Comparison M5 y M5B ################
modelL5 <- list()
modelL5$m5 <- LCSM5
modelL5$m5B <- LCSM5B
fitL5 <- list()
fitL5$m5 <- lavaan(data=bob, modelL5$m5, missing="ML")
fitL5$m5B <- lavaan(data=bob, modelL5$m5B, missing="ML")
summ_list <- list()
fitms <- c("npar","logl","aic",
           "chisq", "df", "pvalue",
           "rmsea", "cfi", "srmr")
fitmeasL5 <- list()
fitmeasL5$m5 <- fitmeasures(fitL5$m5, fit.measures = fitms)
fitmeasL5$m5B <- fitmeasures(fitL5$m5B, fit.measures = fitms)

sum55b <- cbind(round(fitmeasL5$m5,3), round(fitmeasL5$m5B,3))

#Summary
sum55b 
comp55b <- compareFit(fitL5$m5, fitL5$m5B, nested=TRUE) 
#Comparison
summary(comp55b)

# The increment of maladjustment is not significant
# M5B fits as well as M5 with less parameters

#-------- Results from all analysis ---------
#____________________________________________

# Fit index from nested models (M0:M6)
summodels
# Comparison of nested models except M6 (it does not converge)
summary(compmodels)
# Fit index from M5 and M5B 
sum55b
# M5 and M5B comparison
summary(comp55b)
# M5B parameters to interpret
estM5B
