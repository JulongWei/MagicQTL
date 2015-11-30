#scan genome for the simulated data. 
simScan<-function(pathout,gen,map,kk.eigen,model,nfounders,window,step,num=100,mc=10){
  s0<-Sys.time()
  pathdata<-"./Data/Simulation/Simul7qtls/"
  dir.create(pathout,recursive=TRUE,show=FALSE)
##parallel the fucntion
   mc<-mc
   win<-windowsize    
   tmp<-mclapply(1:num,function(i){
     #
     simulfile<-paste(pathdata,"simuldata",i-1,".csv",sep="")
     parmfile<-paste(pathout,model,".",i-1,".parr.csv",sep="")
     bluppfile<-paste(pathout,model,".",i-1,".blupp.csv",sep="")
     #
     X1<-read.csv(file=simulfile)
     n<-nrow(X1)
     X2<-rep(1,n)
     d<-data.frame(y=X1$Phe,x=X2) 
     Y0<-magicScan(dataframe=d,gen=gen,map=map,kk.eigen=kk.eigen,model=model,nfounders=nfounders,window=window,step=step)[[1]]
     Y1<-Y0$parr
     Y2<-Y0$blupp
     write.csv(Y1,file=parmfile,row.names=FALSE)
     write.csv(Y2,file=bluppfile,row.names=FALSE)
   },mc.cores=mc)
   cat("All simulations",num,"have been done","\n")
   return(0)
}     