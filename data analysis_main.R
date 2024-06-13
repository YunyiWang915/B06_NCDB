###############################
######  import packages  ######
###############################
rm(list=ls())
packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos='https://cran.rstudio.com/')
  sapply(pkg, require, character.only = TRUE)
}
packages(c("utils", "arsenal", "dplyr", "gtsummary"))
######################################################################

###############################
########  import data  ########
###############################
path <- "/Users/mojamoja0218max/Desktop/Project2/application/"

# 1) B06 dataset
data_B06 <- read.csv(paste0(path, "nsabpb06.csv"))
colnames(data_B06)
data_B06 <- data_B06 %>% select("RACE", "AGE", "TRT", "PNOD", "MCSIZE", "ER4", "PR4", "DEATH", "DSS", "TLFUP")
dim(data_B06) # m = 1851

data_B06_new <- data_B06 %>% select("AGE", "DEATH", "TLFUP")

data_B06_new$TLFUP <- data_B06_new$TLFUP / 12 # convert months to years
data_B06_new$RACE <- NA
data_B06_new$RACE[data_B06$RACE == 1] <- "1"
data_B06_new$RACE[data_B06$RACE == 2] <- "2"
data_B06_new$RACE[data_B06$RACE == 3] <- "3"
data_B06_new$RACE <- factor(data_B06$RACE, levels = c("1", "2", "3"), labels = c("White", "Black", "Other"))

data_B06_new$TRT <- factor(data_B06$TRT, levels = c("1", "3", "2"), labels = c("TM", "BCT", "Other"))

data_B06_new$PNOD <- NA
data_B06_new$PNOD[data_B06$PNOD > 0] <- "Positive"
data_B06_new$PNOD[data_B06$PNOD == 0] <- "Negative"

data_B06_new$TUMOR_SIZE <- NA
data_B06_new$TUMOR_SIZE[data_B06$MCSIZE < 2] <- "1"
data_B06_new$TUMOR_SIZE[data_B06$MCSIZE >= 2] <- "2"
data_B06_new$TUMOR_SIZE <- factor(data_B06_new$TUMOR_SIZE, levels = c("1", "2"), labels = c("<2cm", ">=2cm"))

data_B06_new$HR <- "Negative"
data_B06_new$HR[data_B06$ER4 >= 1 | data_B06$PR4 >= 1] <- "Positive"

data_B06_new$CD <- NA
data_B06_new$CD[data_B06$DSS == 3] <- "0"
data_B06_new$CD[data_B06$DSS == 4] <- "1"

# 2) NCDB dataset
data_NCDB <- read.csv(paste0(path, "forarctsept2019.csv"))
data_NCDB <- data_NCDB %>% select("RACE", "AGE", "Treatment", "REGIONAL_NODES_POSITIVE", "TUMOR_SIZE", "HR", "death", "DX_LASTCONTACT_DEATH_MONTHS")

data_NCDB_new <- data_NCDB %>% select("AGE", "HR")

data_NCDB_new$TLFUP <- data_NCDB$DX_LASTCONTACT_DEATH_MONTHS / 12 # convert months to years

data_NCDB_new$RACE <- "3"
data_NCDB_new$RACE[data_NCDB$RACE == 1] <- "1"
data_NCDB_new$RACE[data_NCDB$RACE == 2] <- "2"
data_NCDB_new$RACE <- factor(data_NCDB_new$RACE, levels = c("1", "2", "3"), labels = c("White", "Black", "Other"))

data_NCDB_new$TRT <- factor(data_NCDB$Treatment, levels = c("Mastectomy", "BCT", "BCS"), labels = c("TM", "BCT", "Other"))

data_NCDB_new$PNOD <- NA
data_NCDB_new$PNOD[data_NCDB$REGIONAL_NODES_POSITIVE == 0] <- "Negative"
data_NCDB_new$PNOD[data_NCDB$REGIONAL_NODES_POSITIVE %in% c(90, 95, 97)] <- "Positive"

data_NCDB_new$DEATH <- data_NCDB$death

data_NCDB_new$TUMOR_SIZE <- NA
data_NCDB_new$TUMOR_SIZE[data_NCDB$TUMOR_SIZE < 2 | data_NCDB$TUMOR_SIZE %in% c(989, 990, 991, 992)] <- "1"
data_NCDB_new$TUMOR_SIZE[data_NCDB$TUMOR_SIZE >= 2 | data_NCDB$TUMOR_SIZE %in% c(993, 994, 995, 998)] <- "2"
data_NCDB_new$TUMOR_SIZE <- factor(data_NCDB_new$TUMOR_SIZE, levels = c("1", "2"), labels = c("<2cm", ">=2cm"))


###############################
######## data summary  ########
###############################
formula <- as.formula(factor(`DEATH`, levels = c(0, 1), labels = c("Alive", "Dead")) ~ .)
mycontrol <- tableby.control(numeric.test = "kwt", cat.test = "chisq", simulate.p.value = TRUE,
                             numeric.stats = c("Nmiss", "median", "iqr", "range"),
                             cat.stats = c("Nmiss", "countpct"),
                             stats.labels = list(N = "Count", Nmiss = "Missing Values (N)",
                                                 median = "Median", iqr = "IQR", range = "Range", countpct = "Count(Pct)"),
                             digits = 1)
tab1 <- tableby(formula, data = data_B06_new, control = mycontrol)
labels(tab1)[c(3,6,9)] <- c("Time to Last Follow-up", "Positive Axillary Nodes", "Cause of Death")
summary(tab1, text = TRUE)

tab2 <- tableby(formula, data = data_NCDB_new, control = mycontrol)
labels(tab2)[c(4,7)] <- c("Time to Last Follow-up", "Positive Axillary Nodes")
summary(tab2, text = TRUE)


###############################
########## functions ##########
###############################
################ Estimating function for gamma ################
gamma.fun = function(data, ind.wt){
  mydata = data[data$DEATH == 1,]
  TLFUP = mydata$TLFUP
  TLFUP_0 = log(TLFUP)
  TLFUP_2 = TLFUP^2
  
  if (ind.wt == 0){
    gamma.mat = data.frame(cbind(TLFUP_0,TLFUP,TLFUP_2, mydata[c("CD","AGE","PNOD","TUMOR_SIZE","RACE","TRT","HR")]))
    mod = glm((1 - as.numeric(CD)) ~ .,
              data = gamma.mat, family = "binomial")
  }
  if (ind.wt == 1){
    gamma.mat = data.frame(cbind(TLFUP_0,TLFUP,TLFUP_2, mydata[c("CD","AGE","PNOD","TUMOR_SIZE","RACE","TRT","HR","rwt")]))
    mod = glm((1 - as.numeric(CD)) ~ .,
              data = gamma.mat, family = "binomial", weights = rwt)
  }
  
  coeff = coef(summary(mod))[,1]
  coeff_se = coef(summary(mod))[,2]
  
  vcov = vcov(mod)
  
  return(list(coeff = coeff, coeff_se = coeff_se, vcov = vcov))
}


################ Estimating function for beta ################
# negative log partial likelihood function
neglog_pl_beta = function(odat, beta, gammas.hat, events_id, ind.wt){
  Y = odat$Y
  Y.mat = cbind(rep(1,length(Y)), log(Y), Y, Y^2)
  
  gamma1 = gammas.hat[1:4]; gamma2 = gammas.hat[-c(1:4)]
  
  # Estimate rho0
  rho0.hat = exp(as.vector(Y.mat%*%gamma1))
  odat = cbind(odat,rho0.hat)
  
  log_parlike = 0
  for (l in 1:length(events_id)){
    id = events_id[l]
    Yi = odat$Y[id] # Yi
    Xi = c(odat$X1[id], odat$X2[id]) # Xi
    Wi = odat$rwt[id]
    risk.id = which(odat$Y >= Yi) # risk set R(Yi)
    Xj = cbind(odat$X1[risk.id], odat$X2[risk.id])
    
    rho0.hat = odat$rho0.hat[id]
    
    if (ind.wt == 0){
      log_parlike = log_parlike - Xi%*%beta - log(1+rho0.hat*exp(Xi%*%gamma2)) + log(sum(exp(Xj%*%beta)*(1+rho0.hat*exp(Xj%*%gamma2))))
    }
    if (ind.wt == 1){
      log_parlike = log_parlike - Wi*(Xi%*%beta + log(1+rho0.hat*exp(Xi%*%gamma2)) - log(sum(exp(Xj%*%beta)*(1+rho0.hat*exp(Xj%*%gamma2)))))
    }
  }
  return(log_parlike)
}


################ Logit link function ################
g.logit <- function(xx){exp(xx) / (1 + exp(xx))}
# gd.logit <- function(xx){exp(xx) / (1+exp(xx))^2} # derivative of logit link function


###############################
######## data analysis ########
###############################
#1) Estimate gamma's and their se
mod_gamma = gamma.fun(data = data_B06_new, ind.wt = 0)
gammas.hat = mod_gamma$coeff
gammas.se = mod_gamma$coeff_se
gammas.cov = mod_gamma$vcov
  
#2) Generate observational study dataset and order it by observed event time
odat = odat[order(odat$Y),] # order by event time
events_id = which(odat$delta == 1) # observed events index
  
#4) Estimate beta's
par1 = c(0.1,0.1)
beta.hat = nlm(f = neglog_pl_beta, par1, odat = odat, ind.rho_est, gammas.hat = gammas.hat, 
               events_id = events_id, ind.wt = 0)$estimate # Use nlm() to minimize negative log partial likelihood
  
#5) Compute beta's se using perturbation
beta.hat.pr = data.frame(matrix(nrow = 200, ncol = 2)) 
  for (k in 1:200) {
    cdat.pr = cdat
    cdat.pr$rwt = rexp(Nc,1) # generating random weights using standard exponential distribution
    mod_gamma.pr = gamma.fun(data = cdat.pr, ind.wt = 1)
    gammas.hat.pr = mod_gamma.pr$coeff
    
    odat.pr = odat
    odat.pr$rwt = rexp(No,1)
    odat.pr = odat.pr[order(odat.pr$Y),]
    events_id.pr = which(odat.pr$delta == 1)
    
    beta.hat.pr[k,] = nlm(f = neglog_pl_beta, par1, odat = odat.pr, ind.rho_est, gammas.hat = gammas.hat.pr, 
                          events_id = events_id.pr, ind.wt = 1)$estimate
}
  
beta.se = apply(beta.hat.pr,2,sd)
coeff.est = c(gammas.hat, beta.hat, gammas.se, beta.se)
return(list(coeff_est = coeff.est,
              gamma1_cov = gammas.cov[1:4,1:4]))



