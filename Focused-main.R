#focused IVW method
#input: effect size on the outcome "beta.out" and its standard error "se.out"; 
#        effect size on the exposure "beta.exp" and its standard error "se.exp"; 
#        threshold for focusing "tau.f", threshold for strong IV "tau.s", significance level "alpha"
#output: "Vtau": the focused set; "rej.region": the rejection region; "pval": the p-value; 
#         "rej": reject(true) or not; IVW.focus: the value of test statistic.
Focused.ivw.test<-function(beta.out, beta.exp, se.out, se.exp, tau.f=NULL,tau.s=NULL, alpha=0.95){
  p<-length(beta.out)
  if(is.null(tau.f)){
    tau.f= qnorm(0.9)
  }
  Vtau<-which(abs(beta.out)/se.out <= tau.f)
  if(is.null(tau.s)){
    tau.s= qnorm(1-1/p)
  }
  Vtau<-intersect(Vtau, which(abs(beta.exp)/se.exp>=tau.s))
  if(length(Vtau)<=2){
    return(list(Vtau=Vtau, ctr.bs=NA, 
                rej= TRUE, Q.focus=NA, pval=0))
  }
  wts0<-(beta.exp[Vtau]/se.out[Vtau])^2
  wts<-1/beta.exp[Vtau]*wts0/sum(wts0)
  Q.focus=sum(beta.out[Vtau]*wts)
  V.trunc<-1-2*tau.f*dnorm(tau.f)/(pnorm(tau.f)-pnorm(-tau.f))
  sd.psi<-sqrt(V.trunc)/sqrt(sum(wts0))*1.2 #adjusted for remaining bias
  pval.psi=2*(1-pnorm(abs(Q.focus)/sd.psi))
  list(Vtau=Vtau, rej.region=qnorm(0.5+alpha/2)*sd.psi*c(-1,1), pval=pval.psi,
       rej=(pval.psi<=1-alpha), IVW.focus=Q.focus)

}

#focused Median method
#input: effect size on the outcome "beta.out" and its standard error "se.out"; 
#        effect size on the exposure "beta.exp" and its standard error "se.exp"; 
#        threshold for focusing "tau.f", threshold for strong IV "tau.s", significance level "alpha";
#         weight.order: allows for weighted median, simple weight with value -1; "Nbs": rounds of bootstrap.
#output: "Vtau": the focused set; "rej.region": the rejection region; "pval": the p-value; 
#         "rej": reject(true) or not; IVW.focus: the value of test statistic.
Focused.med.test<-function(beta.out, beta.exp, se.out, se.exp, tau.f=NULL, tau.s=NULL, alpha=0.95,
                     weight.order=-1, Nbs=300){
  p<-length(beta.out)
  if(is.null(tau.f)){
    tau.f= qnorm(0.9)
  }
  Vtau<-which(abs(beta.out)/se.out <= tau.f)
  if(is.null(tau.s)){
    tau.s= qnorm(1-1/p)
  }
  Vtau<-intersect(Vtau, which(abs(beta.exp)/se.exp>=tau.s))
  if(length(Vtau)<=2){
    return(list(Vtau=Vtau, ctr.bs=NA, 
                rej= TRUE, Q.focus=NA, pval=0))
  }
  wts<-(beta.exp)^weight.order
  Q.focus=median((beta.out[Vtau])*wts[Vtau])
  inf.bs<-Boot.med(wts=wts, b_out=beta.out, se_out=se.out,
                     Nbs=Nbs, tau=tau.f, Q.focus=Q.focus, Vtau=Vtau, alpha=alpha)
  list(Vtau=Vtau, rej.region=inf.bs$ctr.bs, pval=inf.bs$pval,
       rej=(inf.bs$pval<=1-alpha), Med.focus=Q.focus)
}
Boot.med<-function(b_out,se_out, wts, Nbs, tau, Q.focus, Vtau=NULL, alpha) {
  Q.bs <- rep(0, Nbs)
  a.truc<- -tau
  b.truc<- tau
  Z.truc<-pnorm(b.truc)-pnorm(a.truc)
  
  for (i in 1:Nbs) {
    z.tau<-qnorm(runif(length(Vtau))*Z.truc+pnorm(a.truc))-
      (dnorm(a.truc)-dnorm(b.truc))/Z.truc
    z.tau<-z.tau*1.2 #adjusted for remaining bias
    
    Q.bs[i]<- median(wts[Vtau]*se_out[Vtau]*z.tau)
  }
  list(ctr.bs=quantile(Q.bs,c(0.5-alpha/2, 0.5+alpha/2)),pval=mean(abs(Q.bs)>=abs(Q.focus)))
}


