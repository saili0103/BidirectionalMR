source('~/Dropbox/Ting-Sai/Bi-directional MR/Focused-main.R')
library(TwoSampleMR)
library(MendelianRandomization)
set.seed(123)
#test the causal effect of BMI on T2D
dat<-readRDS('bmi.t2d.rds')
p<-length(dat$beta.out)
re.fmed.dy<-Focused.med.test(beta.out=dat$beta.out, se.out=dat$se.out, beta.exp=dat$beta.exp, se.exp=dat$se.exp,
                          tau.f=1.5, tau.s=qnorm(1-1/p)) #focused median
re.fivw.dy<-Focused.ivw.test(beta.out=dat$beta.out, se.out=dat$se.out, beta.exp=dat$beta.exp, se.exp=dat$se.exp,
                     tau.f=1.5, tau.s=qnorm(1-1/p)) #focused IVW
Sd.hat<-which(abs(dat$beta.exp)/dat$se.exp>=qnorm(1-1/p))
re.med.dy<-MendelianRandomization::mr_median(mr_input(bx = dat$beta.exp[Sd.hat], bxse = dat$se.exp[Sd.hat],
                                                     by = dat$beta.out[Sd.hat], byse = dat$se.out[Sd.hat]), 
                                            iterations=1000) #MR-Median
re.egger.dy<-MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exp[Sd.hat], bxse = dat$se.exp[Sd.hat],
                                          by = dat$beta.out[Sd.hat], byse = dat$se.out[Sd.hat]))#MR-Egger
re.ivw.dy<-MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exp[Sd.hat], bxse = dat$se.exp[Sd.hat],
                                        by = dat$beta.out[Sd.hat], byse = dat$se.out[Sd.hat])) #MR-IVW


#test the causal effect of T2D on BMI
dat<-readRDS('t2d.bmi.rds')
re.fmed.yd<-Focused.med.test(beta.out=dat$beta.out, se.out=dat$se.out, beta.exp=dat$beta.exp, se.exp=dat$se.exp,
                             tau.f=1.5, tau.s=qnorm(1-1/p))#focused median
re.fivw.yd<-Focused.ivw.test(beta.out=dat$beta.out, se.out=dat$se.out, beta.exp=dat$beta.exp, se.exp=dat$se.exp,
                             tau.f=1.5, tau.s=qnorm(1-1/p))#focused IVW
Sd.hat<-which(abs(dat$beta.exp)/dat$se.exp>=qnorm(1-1/p))
re.med.yd<-MendelianRandomization::mr_median(mr_input(bx = dat$beta.exp[Sd.hat], bxse = dat$se.exp[Sd.hat],
                                                      by = dat$beta.out[Sd.hat], byse = dat$se.out[Sd.hat]), 
                                             iterations=1000) #MR-Median
re.egger.yd<-MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exp[Sd.hat], bxse = dat$se.exp[Sd.hat],
                                                       by = dat$beta.out[Sd.hat], byse = dat$se.out[Sd.hat]))#MR-Egger
re.ivw.yd<-MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exp[Sd.hat], bxse = dat$se.exp[Sd.hat],
                                                   by = dat$beta.out[Sd.hat], byse = dat$se.out[Sd.hat])) #MR-IVW

