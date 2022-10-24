dat<-readRDS("sim1.ldl.cad.rds")
intersect(dat$SNP[dat$dataset=='d2y'],dat$SNP[dat$dataset=='y2d'])
source("Focused-main.r")
library(MendelianRandomization)
pi.D <- dat$beta.d
pi.Y <- dat$beta.y
pi.D.sd<-dat$se.d
pi.Y.sd<- dat$se.y
p<-length(pi.D)
Vdy<-which(pi.Y==0 & pi.D!=0)
Vyd<-which(pi.D==0 & pi.Y!=0)
S.D<-which(pi.D!=0)
S.Y<-which(pi.Y!=0)

set.seed(123)
#data generation
betady = -0.2 #or 0 or 0.2
p<-length(pi.D)
Niter= 3000
for(betayd in c(-0.2, 0, 0.2)){
  B<-matrix(c(1,-betady, -betayd,1),ncol=2)
  gam.Y<-cbind(pi.Y,pi.D)%*%solve(B)[,1]
  gam.D<-cbind(pi.Y,pi.D)%*%solve(B)[,2]
  sig.gamY<-(pi.Y.sd^2 + betady^2*pi.D.sd^2)/(1-betady*betayd)^2
  sig.gamD<-(pi.D.sd^2 + betayd^2*pi.Y.sd^2)/(1-betady*betayd)^2
  prop.mat1<-matrix(NA, ncol=6,nrow=Niter) #focused ivw with tau_f=1.5
  prop.mat2<-matrix(NA, ncol=6,nrow=Niter)#focused ivw with tau_f=1.2
  prop.mat3<-matrix(NA, ncol=6,nrow=Niter)#focused median with tau_f=1.5
  other.mat<-matrix(NA, ncol=6,nrow=Niter)
  for(iter in 1:Niter){
    sig.gamY<-(pi.Y.sd^2 + betady^2*pi.D.sd^2)/(1-betady*betayd)^2
    sig.gamD<-(pi.D.sd^2 + betayd^2*pi.Y.sd^2)/(1-betady*betayd)^2
    gam.D.hat<-sapply(1:p, function(k) rnorm(1,gam.D[k], sqrt(sig.gamD[k])))
    gam.Y.hat<-sapply(1:p, function(k) rnorm(1,gam.Y[k], sqrt(sig.gamY[k])))
    ks.dy<-Focused.ivw.test(beta.out=gam.Y.hat, se.out=sqrt(sig.gamY), beta.exp=gam.D.hat, se.exp=sqrt(sig.gamD),
                  tau.f=1.5, tau.s=qnorm(1-1/p))
    ks.yd<-Focused.ivw.test(beta.out=gam.D.hat, se.out=sqrt(sig.gamD), beta.exp=gam.Y.hat, se.exp=sqrt(sig.gamY),
                  tau.f=1.5, tau.s=qnorm(1-1/p))
    prop.mat1[iter,]<-c(ks.dy$rej, ks.yd$rej, length(intersect(Vdy,ks.dy$Vtau))/length(ks.dy$Vtau), 
                     length(intersect(Vyd,ks.yd$Vtau))/length(ks.yd$Vtau), ks.dy$IVW.focus, ks.yd$IVW.focus)
    ks.dy2<-Focused.ivw.test(beta.out=gam.Y.hat, se.out=sqrt(sig.gamY), beta.exp=gam.D.hat, se.exp=sqrt(sig.gamD),
                  tau.f=1.2, tau.s=qnorm(1-1/p))
    ks.yd2<-Focused.ivw.test(beta.out=gam.D.hat, se.out=sqrt(sig.gamD), beta.exp=gam.Y.hat, se.exp=sqrt(sig.gamY),
                  tau.f=1.2, tau.s=qnorm(1-1/p))
    prop.mat2[iter,]<-c(ks.dy2$rej, ks.yd2$rej, length(intersect(Vdy,ks.dy2$Vtau))/length(ks.dy2$Vtau),
                     length(intersect(Vyd,ks.yd2$Vtau))/length(ks.yd2$Vtau), ks.dy2$IVW.focus, ks.dy2$IVW.focus)
    ks.dy3<-Focused.med.test(beta.out=gam.Y.hat, se.out=sqrt(sig.gamY), beta.exp=gam.D.hat, se.exp=sqrt(sig.gamD),
                   tau.f=1.5, tau.s=qnorm(1-1/p))
    ks.yd3<-Focused.med.test(beta.out=gam.D.hat, se.out=sqrt(sig.gamD), beta.exp=gam.Y.hat, se.exp=sqrt(sig.gamY),
                   tau.f=1.5, tau.s=qnorm(1-1/p))
    prop.mat3[iter,]<-c(ks.dy3$rej, ks.yd3$rej, length(intersect(Vdy,ks.dy3$Vtau))/length(ks.dy3$Vtau),
                        length(intersect(Vyd,ks.yd3$Vtau))/length(ks.yd3$Vtau),ks.dy3$Med.focus, ks.yd3$Med.focus)
    
     Sd.hat<-which(abs(gam.D.hat)/sqrt(sig.gamD) >=  qnorm(1-1/p))
     Sy.hat<-which(abs(gam.Y.hat)/sqrt(sig.gamY) >=  qnorm(1-1/p))
     egger.re<-MendelianRandomization::mr_egger(mr_input(bx = gam.D.hat[Sd.hat], bxse = sqrt(sig.gamD[Sd.hat]),
                                  by = gam.Y.hat[Sd.hat], byse = sqrt(sig.gamY[Sd.hat])))
     pval.egg.dy<-egger.re$Pvalue.Est
     median.re.dy<-MendelianRandomization::mr_median(mr_input(bx = gam.D.hat[Sd.hat], bxse = sqrt(sig.gamD[Sd.hat]),
                         by = gam.Y.hat[Sd.hat], byse = sqrt(sig.gamY[Sd.hat])), iterations=500)
     pval.med.dy<-median.re.dy$Pvalue
     
     egger.re<-MendelianRandomization::mr_egger(mr_input(bx = gam.Y.hat[Sy.hat], bxse = sqrt(sig.gamY[Sy.hat]),
                                  by = gam.D.hat[Sy.hat], byse = sqrt(sig.gamD[Sy.hat])))
     pval.egg.yd<-egger.re$Pvalue.Est
     median.re.yd<-MendelianRandomization::mr_median(mr_input(bx = gam.Y.hat[Sy.hat], bxse = sqrt(sig.gamY[Sy.hat]),
                           by = gam.D.hat[Sy.hat], byse = sqrt(sig.gamD[Sy.hat])), iterations=500)
     pval.med.yd<-median.re.yd$Pvalue
     other.mat[iter,]<-c(pval.med.dy<=0.05, pval.med.yd<=0.05,pval.egg.dy<=0.05,pval.egg.yd<=0.05, 
                       length(intersect(Vdy,Sd.hat))/length(Sd.hat), length(intersect(Vyd,Sy.hat))/length(Sy.hat))
   }
  cat('iter=',iter, colMeans(prop.mat1[1:iter,]), '\n')
  cat('iter=',iter, colMeans(prop.mat2[1:iter,]), '\n')
  cat('iter=',iter, colMeans(prop.mat3[1:iter,]), '\n')
}


