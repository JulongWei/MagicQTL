#
##Scan the simulated data from the null model and 
##Produce the threshold value in 100%,99%,95%,90% levels.  
permutScanSim<-function(pathout,gen,map,kk.eigen,model,nfounders,window,step,num,mc){
##
   dir.create(pathout,recursive=TRUE,show=FALSE)
   kk<-kk.eigen[[1]]
   GG<-list(gen[[1]],kk)  
##parallel 
   tmp<-mclapply(1:num,function(i){       
     iseed<-0
     while (TRUE) {    
       set.seed(seed=i+iseed)
       X0<-simulH2(GG=GG,addqtl=FALSE)             
       iseed<-iseed+1
       X1<-X0$Parms
       X2<-X0$Datas
       n<-nrow(X2)
       d<-data.frame(y=X2$Phe,x=rep(1,n))
       Y0<-try(magicScan(dataframe=d,gen=gen,map=map,kk.eigen=kk.eigen,model=model,nfounders=nfounders,window=window,step=step)[[1]])
       if ( !inherits(Y0,"try-error")) break  
     }  
     Y1<-Y0$parr
     Y2<-Y0$blupp
##lrt           
     i0<-4
     Qvalue<-quantile(Y1[,i0],probs=c(1,0.995,0.99,0.98,0.95),na.rm=TRUE,names=FALSE)
     i1<-unique(which.min(abs(Y1[,i0]-Qvalue[1])))
     i2<-unique(which.min(abs(Y1[,i0]-Qvalue[2])))
     i3<-unique(which.min(abs(Y1[,i0]-Qvalue[3])))
     i4<-unique(which.min(abs(Y1[,i0]-Qvalue[4])))
     i5<-unique(which.min(abs(Y1[,i0]-Qvalue[5])))
     YY1<-c(Y1[i1,i0],Y1[i2,i0],Y1[i3,i0],Y1[i4,i0],Y1[i5,i0])
##lrt.log.pvalue             
     i0<-6
     Qvalue<-quantile(Y1[,i0],probs=c(1,0.995,0.99,0.98,0.95),na.rm=TRUE,names=FALSE)
     i1<-unique(which.min(abs(Y1[,i0]-Qvalue[1])))
     i2<-unique(which.min(abs(Y1[,i0]-Qvalue[2])))
     i3<-unique(which.min(abs(Y1[,i0]-Qvalue[3])))
     i4<-unique(which.min(abs(Y1[,i0]-Qvalue[4])))
     i5<-unique(which.min(abs(Y1[,i0]-Qvalue[5])))
     YY2<-c(Y1[i1,i0],Y1[i2,i0],Y1[i3,i0],Y1[i4,i0],Y1[i5,i0])  
     cat(i,"\n")
     XX1<-matrix(X1,1,length(X1))
     YY1<-matrix(YY1,1,length(YY1)
     YY2<-matrix(YY2,1,length(YY2))
     return(list(XX1,YY1,YY2))
   },mc.cores=mc)
   
   ###
   Xparms<-lapply(1:num,function(i) return(tmp[[i]][[1]]) )
   Xparms<-as.data.frame(do.call(rbind,Xparms))
   #
   Yparms1<-lapply(1:num,function(i) return(tmp[[i]][[2]]) )
   Yparms1<-as.data.frame(do.call(rbind,Yparms1))
   #
   Yparms2<-lapply(1:num,function(i) return(tmp[[i]][[3]])
   Yparms2<-as.data.frame(do.call(rbind,Yparms2))
   #
   names(Xparms)<-c("mu","v.qtl.all","v.polygen","v.g","ve","vp")
   names(Yparms1)<-c("max.1","max.0.99","max.0.98","max.0.95")
   names(Yparms2)<-c("max.1","max.0.995","max.0.99","max.0.98","max.0.95")
##
   xfn<-paste(pathout,"dataparms.permut.csv",sep="")     
   write.csv(Xparms,file=xfn,row.names=FALSE)
   yfn1<-paste(pathout,"maxparms.lrt.permut.csv",sep="")     
   write.csv(Yparms1,file=fn,row.names=FALSE)
   yfn2<-paste(pathout,"maxparms.logp.permut.csv",sep="")
   write.csv(Yparms2,file=yfn2,row.names=FALSE)
##
   cat("Permuts",num,"have been done",0,"\n")
   return(0)   
}
