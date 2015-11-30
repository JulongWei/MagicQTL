# 
simH2<-function(seed=NULL,nfounders=8,GG,addqtl=FALSE,qtl=NULL){  
   set.seed(seed=seed)   
#    
   mu<-10
   sg2<-0.5
   se2<-0.5
   r<-nfounders
   #
   if (is.null(qtl)){
     locus<-NULL
     g.locus<-NULL
   }else{
     locus<-as.numeric(qtl[,1])
     g.locus<-qtl[,-1]
   } 
   gen<-GG[[1]]
   kk<-GG[[2]] 
   n<-ncol(gen)
   n.qtl<-length(locus)
    
##simul-background genetics
   sg<-sqrt(sg2)
   G<-t(chol(kk))%*%matrix(rnorm(n,0,1))*sg
   G<-G*as.numeric(sqrt(sg2/var(G)))   
##simul-qtl    
   gg<-matrix(0,n,1)
   g.qtls<-NULL
   vg<-numeric()
   if (addqtl){
     for ( i in 1:n.qtl){              
       ll<-locus[i]
       g1<-matrix(g.locus[i,],r,1)
       sub.1<-seq((ll-1)*r+1,(ll-1)*r+r)
       gg.1<-t(gen[sub.1,])%*%g1
       vg<-c(vg,var(gg.1))
       gg<-gg+gg.1
       g.qtls<-cbind(g.qtls,gg.1)     
     }  
   }
   vg<-c(vg,var(gg),var(G))
   vg<-c(vg,var(gg+G))
#simul-error
   E<-matrix(rnorm(n,0,1))*sqrt(se2)
   Phe<-mu+G+gg+E
   ve<-var(E)
   vp<-var(Phe)
    
##names Y0, that is "parameter of data"
   Y0<-c(mu,vg,ve,vp)
   if (is.null(g.locus)){
     names(Y0)<-c("mu","v.qtl.all","v.polygen","v.g","ve","vp")  
   }else{
     names(Y0)<-c("mu",paste("v.qtl.",1:n.qtl,sep=""),"v.qtl.all","v.polygen","v.g","ve","vp")
   }
    
##Y1,that is "data of individuals"
   Mu<-rep(mu,n)
   if (is.null(g.locus)){
     Y1<-data.frame(Mu,G,E,Phe)
   }else{
     Y1<-data.frame(Mu,g.qtls,gg,G,E,Phe)
     colnames(Y1)<-c("Mu",paste("qtl.",1:n.qtl,sep=""),"qtl.all","G","E","Phe")
   }
#
   result<-list(Parms=Y0,Datas=Y1)
   return(result)
}
               