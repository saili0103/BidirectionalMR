setwd("/Users/saili/Dropbox/Ting-Sai/Bi-directional MR/Code")
source("Focused-main.R")
library(MendelianRandomization)
 #data generation
set.seed(123)
n.D=10^4; n.Y=10^4
betady = 0 #or 0 or 0.2
p=1000
 Niter= 500
 V0<-sample(1:p,p)
 Vdy<-V0[1:(0.3*p)]
 Vyd<-V0[(0.3*p+1):(0.6*p)]
 Vnull<-V0[(0.6*p+1):(0.8*p)]
 Vpl<-V0[(0.8*p+1):p]
 pi.D<-rep(0,p)
 pi.Y<-rep(0,p)
 pi.D[c(Vdy,Vpl)]<-rnorm(0.5*p)*sqrt(2*log(p)/n.D)
 pi.Y[c(Vyd,Vpl)]<-rnorm(0.5*p)*sqrt(2*log(p)/n.Y)
 re1<-NULL; re2<-NULL; re3<-NULL; re.other<-NULL
for(betayd in c(-0.2, 0, 0.2)){
  B<-matrix(c(1,-betady, -betayd,1),ncol=2)
  prop.mat1<-matrix(NA, ncol=6,nrow=Niter) #focused ivw with tau_f=1.5
  prop.mat2<-matrix(NA, ncol=6,nrow=Niter)#focused ivw with tau_f=1.2
  prop.mat3<-matrix(NA, ncol=6,nrow=Niter)#focused median with tau_f=1.5
  other.mat<-matrix(NA, ncol=6,nrow=Niter)
  iter=1
  for(iter in 1:Niter){
    freq <- sample(c(0.1,0.2,0.3), p, replace = TRUE)
    Z<-sapply(1:p, function(x) rbinom(n.D+n.Y, 2, freq[x]))
    Z.D<-Z[1:n.D,]
    Z.Y<- Z[-(1:n.D),]
    U<-cbind(rnorm(n.D,0.1),rnorm(n.Y,0.1))
    U<-U%*%matrix(c(sqrt(1.1)/2+sqrt(0.9)/2, sqrt(1.1)/2-sqrt(0.9)/2,
                    sqrt(1.1)/2-sqrt(0.9)/2, sqrt(1.1)/2+sqrt(0.9)/2),nrow=2)
    D<-Z.D%*%cbind(pi.Y,pi.D)%*%solve(B)[,2]+U%*%solve(B)[,2]
    Y<-Z.Y%*%cbind(pi.Y,pi.D)%*%solve(B)[,1]+U%*%solve(B)[,1]
    #mean((Z.Y%*%cbind(pi.Y,pi.D)%*%solve(B)[,1])^2); mean((U%*%solve(B)[,1])^2)
    #mean((Z.D%*%cbind(pi.Y,pi.D)%*%solve(B)[,2])^2); mean((U%*%solve(B)[,2])^2)

    summ.D<-sapply(1:p, function(k) summary(lm(D~Z.D[,k]))$coef[2,1:2])
    summ.Y<-sapply(1:p, function(k) summary(lm(Y~Z.Y[,k]))$coef[2,1:2])
    gam.D.hat<-summ.D[1,]; gam.Y.hat<-summ.Y[1,]
    se.D<-summ.D[2,]; se.Y<- summ.Y[2,]
    ks.dy<-Focused.ivw.test(beta.out=gam.Y.hat, se.out=se.Y, beta.exp=gam.D.hat, se.exp=se.D,
                  tau.f=1.5, tau.s=qnorm(1-1/p))
    ks.yd<-Focused.ivw.test(beta.out=gam.D.hat, se.out=se.D, beta.exp=gam.Y.hat, se.exp=se.Y,
                  tau.f=1.5, tau.s=qnorm(1-1/p))
    prop.mat1[iter,]<-c(ks.dy$rej, ks.yd$rej, length(intersect(Vdy,ks.dy$Vtau))/length(ks.dy$Vtau),
                     length(intersect(Vyd,ks.yd$Vtau))/length(ks.yd$Vtau), ks.dy$IVW.focus, ks.yd$IVW.focus)
    ks.dy2<-Focused.ivw.test(beta.out=gam.Y.hat, se.out=se.Y, beta.exp=gam.D.hat, se.exp=se.D,
                  tau.f=1.2, tau.s=qnorm(1-1/p))
    ks.yd2<-Focused.ivw.test(beta.out=gam.D.hat, se.out=se.D, beta.exp=gam.Y.hat, se.exp=se.Y,
                  tau.f=1.2, tau.s=qnorm(1-1/p))
    prop.mat2[iter,]<-c(ks.dy2$rej, ks.yd2$rej, length(intersect(Vdy,ks.dy2$Vtau))/length(ks.dy2$Vtau),
                     length(intersect(Vyd,ks.yd2$Vtau))/length(ks.yd2$Vtau), ks.dy2$IVW.focus, ks.dy2$IVW.focus)
    ks.dy3<-Focused.med.test(beta.out=gam.Y.hat, se.out=se.Y, beta.exp=gam.D.hat, se.exp=se.D,
                   tau.f=1.5, tau.s=qnorm(1-1/p))
    ks.yd3<-Focused.med.test(beta.out=gam.D.hat, se.out=se.D, beta.exp=gam.Y.hat, se.exp=se.Y,
                   tau.f=1.5, tau.s=qnorm(1-1/p))
    prop.mat3[iter,]<-c(ks.dy3$rej, ks.yd3$rej, length(intersect(Vdy,ks.dy3$Vtau))/length(ks.dy3$Vtau),
                        length(intersect(Vyd,ks.yd3$Vtau))/length(ks.yd3$Vtau),ks.dy3$Med.focus, ks.yd3$Med.focus)

     Sd.hat<-which(abs(gam.D.hat)/se.D>=  qnorm(1-1/p))
     Sy.hat<-which(abs(gam.Y.hat)/se.Y >=  qnorm(1-1/p))
     egger.re<-MendelianRandomization::mr_egger(mr_input(bx = gam.D.hat[Sd.hat], bxse = se.D[Sd.hat],
                                  by = gam.Y.hat[Sd.hat], byse = se.Y[Sd.hat]))
     pval.egg.dy<-egger.re$Pvalue.Est
     median.re.dy<-MendelianRandomization::mr_median(mr_input(bx = gam.D.hat[Sd.hat], bxse = se.D[Sd.hat],
                         by = gam.Y.hat[Sd.hat], byse = se.Y[Sd.hat]), iterations=500)
     pval.med.dy<-median.re.dy$Pvalue

     egger.re<-MendelianRandomization::mr_egger(mr_input(bx = gam.Y.hat[Sy.hat], bxse = se.Y[Sy.hat],
                                  by = gam.D.hat[Sy.hat], byse = se.D[Sy.hat]))
     pval.egg.yd<-egger.re$Pvalue.Est
     median.re.yd<-MendelianRandomization::mr_median(mr_input(bx = gam.Y.hat[Sy.hat], bxse = se.Y[Sy.hat],
                           by = gam.D.hat[Sy.hat], byse = se.D[Sy.hat]), iterations=500)
     pval.med.yd<-median.re.yd$Pvalue
     other.mat[iter,]<-c(pval.med.dy<=0.05, pval.med.yd<=0.05,pval.egg.dy<=0.05,pval.egg.yd<=0.05,
                       length(intersect(Vdy,Sd.hat))/length(Sd.hat), length(intersect(Vyd,Sy.hat))/length(Sy.hat))
    if(iter>=2){
      cat('iter=',iter, colMeans(prop.mat1[1:iter,]),'\n')
    }

  }
  re1<-c(re1, colMeans(prop.mat1))
  re2<-c(re2, colMeans(prop.mat2))
  re3<-c(re3, colMeans(prop.mat3))
  re.other<-c(re.other, colMeans(other.mat))
}

write.table(re1,'Raw-re1-betady0.txt')
write.table(re2,'Raw-re2-betady0.txt')
write.table(re3,'Raw-re3-betady0.txt')
write.table(re.other,'Raw-re-other-betady0.txt')

