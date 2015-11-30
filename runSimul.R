#
library(MagicQTL)

rm(list=ls())
#I simulate data 
#generate into the 1000 simulated data.
path<-"./Data/Simulation/Simul7qtls/"
dir.create(path=path,recursive=TRUE,show=FALSE)
qtl<-read.csv(file="simul.7qtl.csv")[,-1]
load(file="chrs.magic.RData")
load(file="chrs.kkship.RData")
#
GG<-list(gen[[1]],kk)
num<-1000
parms<-lapply(1:num,function(ii){
  rr<-simulH2(GG=GG,qtl=qtl)
  datafile<-paste(path,"/simuldata.",ii-1,".csv",sep="")
  write.csv(rr[[2]],file=datafile,row.names=FALSE)
  return(rr[[1]])
})
parms<-do.call(rbind,parms)
parmsfile<-paste(path,"/simulparms.csv",sep="")
write.csv(parms,file=parmsfile,row.names=FALSE)



#II scan
rm(list=ls())
source("simScan.R")
source("simH2.R")
load("chrs.kkship.RData")
load("chrs.magic.RData")
num<-1000
models<-c("Random-A","Random-B","Fixed-A","Fixed-B","IM")
nm<-length(models)
mc<-10
for (i in 1:nm){
  pathout<-paste("./Data/Simulation/",models[i],"/",sep="")
  tmp<-simScan(pathout=pathout,gen=gen,map=map,kk.eigen=kk.eigen,model=models[i],nfounders=8,window=5,step=42,num=num,mc=10)
}
#
models<-"CIM"
steps<-c(42,21,14,10.5,8.4,6.4)
cimpre<-c(".co10",".co20",".co30",".co40",".co50",".co65")
num<-1000
mc<-10
nstep<-length(steps)
for (i in 1:nstep){
  pathout<-paste("./Data/Simulation/",models[i],cimpre,"/",sep="")
  tmp<-simScan(pathout=pathout,gen=gen,map=map,kk.eigen=kk.eigen,model=models,nfounders=8,window=5,step=steps[i],num=num,mc=mc)
}
#

#III  permutation to obtain the threshold value
rm(list=ls())
source("permutScanSim")
source("simH2.R")
load("chrs.kkship.RData")
load("chrs.magic.RData")
models<-c("Random-A","Random-B","Fixed-A","Fixed-B","IM")
nm<-length(models)
num<-1000
mc<-10
for (i in 1:nm){
  pathout<-paste("./Data/Simulation/permutation/",models[i],"/",sep="")
  tmp<-permutScanSim(pathout=pathout,gen=gen,map=map,kk.eigen=kk.eigen,model=models[i],nfounders=8,window=5,step=42,num=num,mc=mc)
}

##
models<-"CIM"
steps<-c(42,21,14,10.5,8.4,6.4)
cimpre<-c(".co10",".co20",".co30",".co40",".co50",".co65")
nstep<-length(steps)
num<-1000
mc<-10
for (i in 1:nstep){
  pathout<-paste("./Data/Simulation/permutation/",models[i],cimpre,"/",sep="")
  tmp<-permutScanSim(pathout=pathout,gen=gen,map=map,kk.eigen=kk.eigen,model=models,nfounders=8,window=5,step=steps[i],num=num,mc=mc)
}
