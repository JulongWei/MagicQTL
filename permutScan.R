###
permutScan<-function(pathout,dataframe,gen,map,kk.eigen,model,nfounders,window,step=20,bychr=1,num,mc){
##    
  d<-dataframe
  n0<-nrow(d)
  tmp<-mclapply(1:num,function(i){
    while(TRUE){
      index<-sample(1:n0)
      tmp0<-magicScan(dataframe=d,gen=gen,map=map,kk.eigen=kk.eigen,index=index,model=model,nfounders=nfounders,window=window,step=step,bychr=bychr)
      if (class(tmp0)!="try-error") break
    }
    cat("Num permut:",i,"\n")
    return(tmp0)
  },mc.cores=mc)

   chrnum<-length(gen)
   chr.all<-data.frame()
   parms<-lapply(1:num,function(i){
     Y1<-lapply(1:chrnum,function(ichr){return(tmp[[i]][[ichr]]$parr)})
     Y1<-do.call(rbind,Y1)
     #LRT
     i0<-4
     Qvalue<-quantile(Y1[,i0],probs=c(1,0.995,0.99,0.98,0.95),na.rm=TRUE,names=FALSE)
     i1<-unique(which.min(abs(Y1[,i0]-Qvalue[1])))
     i2<-unique(which.min(abs(Y1[,i0]-Qvalue[2])))
     i3<-unique(which.min(abs(Y1[,i0]-Qvalue[3])))
     i4<-unique(which.min(abs(Y1[,i0]-Qvalue[4])))
     i5<-unique(which.min(abs(Y1[,i0]-Qvalue[5])))
     YY1<-c(Y1[i1,i0],Y1[i2,i0],Y1[i3,i0],Y1[i4,i0],Y1[i5,i0])
     YY1<-matrix(YY1,1,length(YY1)) 
     #logp
     i0<-6
     Qvalue<-quantile(Y1[,i0],probs=c(1,0.995,0.99,0.98,0.95),na.rm=TRUE,names=FALSE)
     i1<-unique(which.min(abs(Y1[,i0]-Qvalue[1])))
     i2<-unique(which.min(abs(Y1[,i0]-Qvalue[2])))
     i3<-unique(which.min(abs(Y1[,i0]-Qvalue[3])))
     i4<-unique(which.min(abs(Y1[,i0]-Qvalue[4])))
     i5<-unique(which.min(abs(Y1[,i0]-Qvalue[5])))
     YY2<-c(Y1[i1,i0],Y1[i2,i0],Y1[i3,i0],Y1[i4,i0],Y1[i5,i0]) 
     YY2<-matrix(YY2,1,length(YY2))
     #
     RR<-list(YY1,YY2)
     return(RR)               
   }
   #
   chrall0<-lapply(1:num,function(i) return(parms[[i]][[1]]))
   chrall0<-as.data.frame(do.call(rbind,chrall0))
   names(chrall0)<-c("max.1","max.0.995","max.0.99","max.0.98","max.0.95")
   logpfile<-paste(pathout,".logp.permut.csv",sep="")
   write.csv(chrall0,file=logpfile,row.names=FALSE)
   #
   chrall1<-lapply(1:num,function(i) return(parms[[i]][[2]]))
   chrall1<-as.data.frame(do.call(rbind,chrall1))
   names(chrall1)<-c("max.1","max.0.995","max.0.99","max.0.98","max.0.95")
   lrtfile<-paste(pathout,".lrt.permut.csv",sep="")
   write.csv(chrall1,file=lrtfile,row.names=FALSE)
   cat("permutations", num, "have been done","\n")
   return(0)
}
