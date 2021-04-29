##@@ Myriam BROSSARD  2020-10
##@@ Example of R script to generate the DCCT-based simulated data under a scenario of complex genetic architecture involving:
##@@ M causal SNPs with direct/indirect associations 
##@@ K=2 time-to-event traits (retinopathy, nephropathy)
##@@ L=3 longitudinal risk factors: 2 observed (HbA1c, SBP) & one simulated (U)

# test2 

set.seed(13857)
dir_wk=c("C:/Users/bross/Documents/GitHub/Joint_Model_for_compleU_genetic_architecture_characterization")
setwd(dir_wk)
long_data=read.delim("/home/bulllab/mbrossard/Joint_models/DCCT_data/2.Data_Processed/DCCT_longdata_for_simu_nomiss.tsv",sep=" ")
#long_data=long_data[,-c(12,13)]
long_data$SEU=long_data$SEx-1
n=length(unique(long_data$PATIENT))
pkg_install <- c("iterators","foreach","doSNOW","snow","survival","eha","mvtnorm","JM","nlme","psych","MASS")
#install.packages("iterators", repos = "http://cran.r-project.org", type="source")

## Parameters for simulations of genotype data for SNP1, SNP5 and SNP3
R<-1000 # no. of replicates
maf_snp1<-0.3
maf_SNP2<-0.1
maf_snp3<-0.4
maf_snp4<-0.3
maf_snp5<-0.2

beta_snp1<-0.7
beta_snp3<-0.8
beta_snp5<-0.4
gamma_snp2<-0.8
gamma_snp4<-0.7
gamma_snp5<-7

alpha_l1k1<-alpha_l1k2<-alpha_l2k2<-0.2
alpha_Uk1<-alpha_Uk2<-0.4

# test 3

#############################################################
### Part 1 : Generation of the genotypes for M=5 SNPs & simulated shared longitudinal risk factors
############################################################# 

## The L longitudinal risk factor values are centered
simdata_QT$hbac <- simdata_QT$hba-mean(long_data$hba) 
simdata_QT$SBPc <- simdata_QT$SBP-mean(long_data$SBP) 

## Simulation of U (unobserved longitudinal trait)
parms_U=list(beta0=12, beta1=0.09, tau0=1.52, tau1=0.3, tau01=-0.5, sigma2=0.5^2)
simdata_QT<-do.call(rbind,by(long_data,INDICES=list(long_data$PATIENT),
  FUN=function(U){
  	U<-U[order(U$obstime),]
  	n<-nrow(U)
  	U<-cbind(rep(1,n),sort(U$obstime))
  	Z<-cbind(rep(1,n),sort(U$obstime))
  	parms_U$Sigma=matriU(c(parms_U$tau0^2,parms_U$tau01*parms_U$tau0*parms_U$tau1,parms_U$tau01*parms_U$tau0*parms_U$tau1, parms_U$tau1^2),nrow=2)
  	vcov<-Z%*%parms_U$Sigma%*%t(Z)+parms_U$sigma2*diag(n)
  	parms_U$beta <-c(parms_U$beta0,parms_U$beta1)
  	U$sim_U <- mvrnorm(1, mu= U%*%parms_U$beta,  Sigma = vcov)
  	return(U)
  }
  ))
simdata_QT$sim_Uc <- simdata_QT$sim_U-mean(simdata_QT$sim_U) 


# Fit the LMM to determine realistic parameters values for the simulation study
lmefitl1 <- lme( hbac ~ obstime , random = ~ 1+obstime|PATIENT, 
                  data = simdata_QT[,c("hbac","obstime","PATIENT")],
                  control = lmeControl(opt = "optim"))
lmefitl2 <- lme( SBPc ~ obstime + SEX , random = ~ 1+obstime|PATIENT,
                  data = simdata_QT[,c("SBPc","obstime","PATIENT","SEX")],
                  control = lmeControl(opt = "optim"))
lmefitU <- lme( sim_Uc ~ obstime, random = ~ 1+obstime|PATIENT, 
                data = simdata_QT[,c("sim_Uc","obstime","PATIENT")],
                control = lmeControl(opt = "optim"))

# Set up parameters values in LMM for SNP1
beta_0<-summary(lmefitl1)$coefficients$fixed[1] 
beta_t<- summary(lmefitl1)$coefficients$fixed[2]  
parms_SNP1A <-list(
  beta=c(beta_0,beta_t,beta_g=beta_snp1),
  Sigma=matriU(c(getVarCov( lmefil1)[1,1],getVarCov(lmefitl1)[1,2],getVarCov( lmefitl1)[1,2],getVarCov( lmefitl1)[2,2]),nrow=2), 
  sigma2=as.numeric(VarCorr(lmefitl1)[3,1]))

# Set up parameters values in LMM for SNP3
parms_SNP3A <-list(
  beta=c(beta_0,beta_t,beta_g=beta_snp3), 
  Sigma=matriU(c(getVarCov( lmefitU)[1,1],getVarCov( lmefitU)[1,2],getVarCov( lmefitU)[1,2],getVarCov( lmefitU)[2,2]),nrow=2), 
  sigma2=as.numeric(VarCorr(lmefitU)[3,1]))

# Set up parameters values in LMM for SNP5
beta_0 <-summary( lmefitl2 )$coefficients$fixed[1] 
beta_t <- summary( lmefitl2 )$coefficients$fixed[2]
beta_sex <- summary( lmefitl2 )$coefficients$fixed[3] 
parms_SNP5A <-list( 
  beta=c(beta_0,beta_t,beta_sex, beta_g=beta_snp5), 
  Sigma=matriU(c(getVarCov( lmefitl2)[1,1],getVarCov(lmefitl2)[1,2],getVarCov(lmefitl2)[1,2],getVarCov(lmefitl2)[2,2]),nrow=2), 
  sigma2=as.numeric(VarCorr(lmefitl2)[3,1]))


pmf_g<-function(g,data,parms,maf,log=F, colnum.QT){
  if(g %in% 0:2){
    pG<-ifelse(g==0,(1-maf)^2,ifelse(g==1,2*maf*(1-maf),maf^2)) 
    n_i<-nrow(data)  
    U<-cbind(rep(1,n_i), data$obstime, rep(g,n_i)) 
	  Z<-cbind(rep(1,n_i),data$obstime) 
    vcov<-Z%*%parms$Sigma%*%t(Z)+parms$sigma2*diag(n_i) 
     f_G_Y<-(log==F)*(dmvnorm(data[,colnum.QT], mean = U%*%parms$beta, sigma = vcov)*pG) + 
       (log==T)*(dmvnorm(data[,colnum.QT], mean = U%*%parms$beta, sigma = vcov, log=T)+ log(pG))
  }else stop("G not in the right range")
  return(f_G_Y)
} # P(Y,G) = P(Y|G) U P(G)

pmf_l2<-function(g,data,parms,maf,log=F, colnum.QT){
  if(g %in% 0:2){
    pG<-ifelse(g==0,(1-maf)^2,ifelse(g==1,2*maf*(1-maf),maf^2)) 
    n_i<-nrow(data)  
    U<-cbind(rep(1,n_i), data$obstime, rep(g,n_i), rep(unique(data$SEX),n_i)) 
	  Z<-cbind(rep(1,n_i),data$obstime) 
    vcov<-Z%*%parms$Sigma%*%t(Z)+parms$sigma2*diag(n_i) 
	  f_G_Y<-(log==F)*(dmvnorm(data[,colnum.QT], mean = U%*%parms$beta, sigma = vcov)*pG) + 
	    (log==T)*(dmvnorm(data[,colnum.QT], mean = U%*%parms$beta, sigma = vcov, log=T)+ log(pG))
  }else stop("G not in the right range")
  return(f_G_Y)
} # P(Y,G) = P(Y|G) U P(G)

Vpmf_g<-Vectorize(pmf_g, vectorize.args = c("g")) # Computation of probability P(Y|G=g) for each possible g 
Vpmf_sbp<-Vectorize(pmf_l2, vectorize.args = c("g")) # Computation of probability P(Y|G=g) for each possible g 

Generate_SNP1 <-do.call(rbind, by(simdata_QT,INDICES = list(simdata_QT$PATIENT),
  function(data,reps,parmsA,parmsN,maf){
    pA<-Vpmf_g(0:2,data,parmsA,maf, colnum.QT=14)
    if(any(pA==0)){ 
      p0A<-Vpmf_g(0:2,data,parmsA,maf,log=T) 
      p0bA<-p0A-maU(p0A) 
      p1A<-eUp(p0bA)/sum(eUp(p0bA))
    }else p1A<-pA/sum(pA)
     pN<-Vpmf_g(0:2,data,parmsN,maf, colnum.QT=14)
    if(any(pN==0)){ 
      p0N<-Vpmf_g(0:2,data,parmsN,maf,log=T) 
      p0bN<-p0N-maU(p0N) 
      p1N<-eUp(p0bN)/sum(eUp(p0bN))
    }else p1N<-pN/sum(pN)
  data.frame(PATIENT=data$PATIENT[1],
             rep=1:reps,
             SNP1A=sample(0:2,reps,prob=p1A,rep=T)) 
  },
reps=R, parmsA=parms_SNP1A,maf=maf_snp1)) # Generation of  P(G|Y)=P(Y|G)/P(Y)

Generate_SNP3<-do.call(rbind,by(simdata_QT,INDICES = list(simdata_QT$PATIENT),
  function(data,reps,parmsA,parmsN,maf){
    pA<-Vpmf_g(0:2,data,parmsA,maf, colnum.QT=13)
    if(any(pA==0)){ 
      p0A<-Vpmf_g(0:2,data,parmsA,maf,log=T) 
      p0bA<-p0A-maU(p0A) 
      p1A<-eUp(p0bA)/sum(eUp(p0bA))
    }else p1A<-pA/sum(pA)
    
    pN<-Vpmf_g(0:2,data,parmsN,maf, colnum.QT=13)
    if(any(pN==0)){ 
      p0N<-Vpmf_g(0:2,data,parmsN,maf,log=T) 
      p0bN<-p0N-maU(p0N) 
      p1N<-eUp(p0bN)/sum(eUp(p0bN))
    }else p1N<-pN/sum(pN)
   
   data.frame(PATIENT=data$PATIENT[1],rep=1:reps,
              SNP3A=sample(0:2,reps,prob=p1A,rep=T)) 
  },reps=R,parmsA=parms_SNP3A, maf=maf_snp3))
 
Generate_SNP5<-do.call(rbind,by(simdata_QT,INDICES = list(simdata_QT$PATIENT),
                                   function(data,reps,parmsA,parmsN,maf){
                                     pA<-Vpmf_sbp(0:2,data,parmsA,maf, colnum.QT=15)
                                     if(any(pA==0)){ 
                                       p0A<-Vpmf_sbp(0:2,data,parmsA,maf,log=T) 
                                       p0bA<-p0A-maU(p0A) 
                                       p1A<-eUp(p0bA)/sum(eUp(p0bA))
                                     }else p1A<-pA/sum(pA)
                                     
                                     pN<-Vpmf_sbp(0:2,data,parmsN,maf, colnum.QT=15)
                                     if(any(pN==0)){ 
                                       p0N<-Vpmf_sbp(0:2,data,parmsN,maf,log=T) 
                                       p0bN<-p0N-maU(p0N) 
                                       p1N<-eUp(p0bN)/sum(eUp(p0bN))
                                     }else p1N<-pN/sum(pN)
                                     data.frame(PATIENT=data$PATIENT[1],
                                                rep=1:reps,
                                                SNP5A=sample(0:2,reps,prob=p1A,rep=T))
                                   },reps=R,parmsA=parms_SNP5A,maf=maf_SNP5))

# Simulation of the M=5 SNPs under the null (assuming HWE)
N=length(unique(simdata_QT$PATIENT))
dat_SNP<-do.call(rbind,lapply(1:R,function(i){data.frame(rep=i,PATIENT=unique(simdata_QT$PATIENT),
	SNP1R=t(c(0,1,2)%*%rmultinom(N, 1, c((1-maf_snp1)^2,2*maf_snp1*(1-maf_snp1),maf_snp1^2))),
	SNP5R=t(c(0,1,2)%*%rmultinom(N, 1, c((1-maf_SNP5)^2,2*maf_SNP5*(1-maf_SNP5),maf_SNP5^2))),
	SNP3R=t(c(0,1,2)%*%rmultinom(N, 1, c((1-maf_snp3)^2,2*maf_snp3*(1-maf_snp3),maf_snp3^2))), 
	SNP4A=t(c(0,1,2)%*%rmultinom(N, 1, c((1-maf_snp4)^2,2*maf_snp4*(1-maf_snp4),maf_snp4^2))),
	SNP5A=t(c(0,1,2)%*%rmultinom(N, 1, c((1-maf_snp5)^2,2*maf_snp5*(1-maf_snp5),maf_snp5^2))),
	SNP4R=t(c(0,1,2)%*%rmultinom(N, 1, c((1-maf_snp4)^2,2*maf_snp4*(1-maf_snp4),maf_snp4^2))),
	SNP5R=t(c(0,1,2)%*%rmultinom(N, 1, c((1-maf_snp5)^2,2*maf_snp5*(1-maf_snp5),maf_snp5^2))))})
	)

# Gather the data
TMP1=merge(Generate_SNP1,Generate_SNP5,by=c("PATIENT", "rep"))
TMP2=merge(TMP1,Generate_SNP3,by=c("PATIENT", "rep"))
TMP3=merge(TMP2,dat_SNP,by=c("PATIENT", "rep"))
TMP4=merge(simdata_QT,TMP3, by=c("PATIENT"))
longQT.SNP=dfOrder(TMP4,c(1,13,10))
#save.image("/home/bulllab/mbrossard/Joint_models/Scenario1_March2018/Simulation_longpart_SNP_longtcentered.RData") # step1
head()

######################################################################
## Part 2:  Simulation of the K=2 non-independent time-to-T1DC traits
# using Weibull survival sub-models
# association structures for HbA1c & SBP: Contemporaneous measures 
######################################################################

# time-to-retinopathy ()
parms_K <-list(
  gamma_SNP2=gamma_SNP2,
  gamma_SNP4=gamma_SNP4, 
  gamma_SNP5=gamma_SNP5,
  nu_k1=1.01,
  lambda_k1=0.1,
  alpha_l1k1=alpha_l1k1,
  alpha_Uk1=alpha_Uk1,
  alpha_t1dur_k1=0.2,
  nu_k2=1.01, 
  lambda_k2=0.01, 
  alpha_l1k2=alpha_l1k2, 
  alpha_l2k2=alpha_l1k2,
  alpha_Uk2=alpha_Uk2, 
  alpha_t1dur.k2=0.2, 
)

survival_K <- NULL

####	
## Simulation of K=2 non-independent time-to-event traits
foreach ( r in 1:R ) { # for each patient i & for each replicate r
  
  longQT.SNP.byrep=longQT.SNP[which(longQT.SNP$rep==r),]
  
  lmefitl1 <- lme( hbac  ~ obstime + SNP1A , random = ~ 1+obstime|PATIENT, 
                 control = lmeControl(opt = "optim"), data = longQT.SNP.byrep) 
  lmefitl1.RE <- summary(lmefitl1)$coefficients$random$PATIENT
  lmefitl1.FE <- summary(lmefitl1)$coefficients$fiUed
  
  lmefitU <- lme( sim_Uc ~ obstime + SNP3A , random = ~ 1+obstime|PATIENT, 
                  control = lmeControl(opt = "optim"), data = longQT.SNP.byrep) 
  lmefitU.RE <- summary(lmefitU)$coefficients$random$PATIENT
  lmefitU.FE <- summary(lmefitU)$coefficients$fiUed

  QT.SNP.byrep<-unique(longQT.SNP.byrep[,c(1,5,7,10,16:29)])
  
  survival_sim <- do.call(rbind,by(QT.SNP.byrep, INDICES = list(QT.SNP.byrep$PATIENT),
                                   function(data, parms, tmin, tmax){
                                     all_iter<-NULL
                                     
                                     lmefitl1.REPAT<-lmefitl1.RE[which(rownames(lmefitl1.RE)==data$PATIENT),]	
                                     trajl1<-function(t) { 
                                       lmefitl1.FE[1] + lmefitl1.REPAT[1] + (lmefitl1.FE[2] + lmefitl1.REPAT[2])*t +
                                       lmefitl1.FE[3]*data$SNP1A }
                                     
                                     lmefitl2.REPAT<-lmefitl2.RE[which(rownames(lmefitl2.RE)==data$PATIENT),]
                                     trajl2<-function(t) {
                                       lmefitl2.FE[1]+ lmefitl2.REPAT[1] + (lmefitl2.FE[2] + lmefitl2.REPAT[2])*t +
                                       lmefitl2.FE[3]*data$SNP5A + lmefitl2.FE[4]*data$SEX	}
                                     
                                     lmefitU.REPAT<-lmefitU.RE[which(rownames(lmefitU.RE)==data$PATIENT),]
                                     trajU<-function(t) { 
                                       lmefitU.FE[1] + lmefitU.REPAT[1] + (lmefitU.FE[2] + lmefitU.REPAT[2])*t +
                                       lmefitU.FE[3]*data$SNP3A }
                                     
                                     ## Simulation of the Time-to-retinopathy outcome (k=1)
                                     baseline_k1<-function(t) {	lambda_k1*nu_k1*t^(nu_k1-1) }
                                     hazard_k1 <- function(t) { 
                                       baseline_k1(t)*exp( alpha_l1k1*traj.l1(t) +  alpha_t1dur_k1*data$DURATION_Years+ 
                                       gamma_SNP2_k1*data$SNP2A + alpha_Uk1*trajU(t) )} 
                                     HazCum_k1 <- function(t) { integrate(hazard_k1, lower=0, upper=t, subdivisions=10000, rel.tol=1.e-05)$value  }
                                     InvHazard_k1 <- function(HazCum_k1, hazard_k1, tmin, tmax) { 
                                       U <- runif(1, 0, 1)
                                       rootfct <- function(t) {   U - exp(-HazCum_k1(t)) }
                                       root <- try(uniroot(rootfct, interval = c(tmin, tmax))$root, silent = TRUE)
                                       root <- if (inherits(root, "try-error")) { tmax + 0.01} else { root }
                                       return(root)								
                                     }
                                     
                                     time.k1 <- InvHazard(HazCum.k1, hazard.k1, tmin, tmaU) 
                                     cumhaz.k1 <- HazCum.k1(time.k1)
                                     survprob.k1 <- exp(-cumhaz.k1 )
                                     
                                     ##Simulation of the Time-to-nephropathy outcome (k=2) 
                                     baseline_k2 <- function(t) { lambda_k2*nu_k2*t^(nu_k2-1) }
                                     hazard_k2 <- function(t) { baseline(t)*exp(alpha_l1k1*traj_l1(t)+
                                                                alpha_l2k2*traj_l2(t) + alpha_Uk2*traj.U(t) + 
                                                                alpha_t1dur_k2*data$DURATION_Years +
                                                                gamma_SNP2*data$SNP2A + gamma_SNP5*data$SNP5A)} 
                                     HazCum_k2 <- function(t) { integrate(hazard_k2, lower=0, upper=t, subdivisions=10000, rel.tol=1.e-05)$value  }
                                     InvHazard_k2 <- function(HazCum.k2, hazard_k2, tmin, tmax) { 
                                       U <- runif(1, 0, 1)
                                       rootfct <- function(t) {   U - exp(-HazCum_k2(t)) }
                                       root <- try(uniroot(rootfct, interval = c(0, tmax))$root, silent = TRUE) 
                                       root <- if (inherits(root, "try-error")) { tmax + 0.01} else { root }
                                       return(root)								}
                                     
                                     time_k2<- InvHazard.k2(HazCum_k2, hazard_k2, tmin, tmax)
                                     cumhaz_k2 <- HazCum.k2(time_k2)
                                     survprob_k2 <- exp(-cumhaz_k2 )
                                     
                                     ## censoring
                                     maxobst <- max(longQT.SNP.byrep[which(longQT.SNP.byrep$PATIENT==data$PATIENT),]$obstime)
                                     cens <- runif(1,min=0,max=maUobst)
                                     event.k1 <- as.numeric(time_k1 <= cens )
                                     t_obs.k1 <- min(time_k1,cens)
                                     event.k2 <- as.numeric(time_k2 <= cens )
                                     t_obs.k2 <- min(time_k2,cens)
                                     
                                     subj_iter<-data.frame(data$PATIENT, r, event_k1, t_obs_k1, event_k2, t_obs_k2 )
                                     all_iter<-rbind(all_iter, subj_iter)
                                     
                                   },parms=parms_K, tmin=1e-5, tmax=100))
  
  survival_K<-rbind(survival_K, survival_sim)
}
colnames(survival_K)<-c("PATIENT", "rep", "event_k1", "tobs_k1","event_k2", "tobs_k2")

####	
# 


survivalK <- merge(survivalk1,survivalk2,by=c("rep","PATIENT"))
merge with the SNPs, SNPR




