#
mpwgaim.cal<-function(pheno.data,mpInterval,nfounders){
   nf<-nfounders
#
  cal.wald<-function(x){
    x<-as.vector(x)
    a<-x[1:nf]
    a<-a-mean(a)
    va<-x[(1:nf)+nf]
    w<-sum(a^2/va)
    return(w)
  }


  asr0 <- asreml(y ~ 1, random = ~ id, data=pheno.data,family=asreml.gaussian(dispersion=1e-03))
  sim.qtl <- mpwgaim(asr0, pheno.data, mpInterval, merge.by = "id",
                  verboseLev=0, gen.type="marker", na.method.X='include',
                  data.name = "Example")
  
if (is.null(sim.qtl$QTL)){
  newmarkers<-NA
  gamma<-t(rep(NA,nf))
  vgamma<-t(rep(NA,nf))
  w0<-0   
}else{
  parms<-sim.qtl$QTL
  markers<-parms$qtl
  mks<-length(markers)
  newmarkers<-sapply(markers,function(one){
    s0<-4
    s1<-nchar(one)
    newmks<-substr(one,s0,s1)
    return(newmks)
    })  
#      
  gamma<-parms$effects
  gamma<-matrix(gamma,mks,nfounders,byrow=TRUE)
  vgamma<-parms$veffects
  vgamma<-matrix(vgamma,mks,nfounders,byrow=TRUE)
  #
  x0<-cbind(gamma,vgamma)
  w0<-as.numeric(apply(x0,1,cal.wald))
  #
}
  rr<-as.data.frame(cbind(newmarkers,gamma,vgamma,w0))
  names(rr)<-c("markers",paste("gamma",1:nfounders,sep=""),paste("vgamma",1:nfounders,sep=""),"wald")
  return(rr)
}
   
