#Call MPWGAIM functions
##-1
#simulation data analysis
##
rm(list=ls())
library(asreml)
library(mpwgaim)
source("mpwgaim.R")
#
pathout<-"./Data/Simulation/mpwgaim/"    
dir.create(path=pathout,recursive=TRUE,showWarnings=FALSE)
pathdata<-"./Data/Simulation/Simul7qtls/"
load(file="mpInterval.RData")
#
pheno<-factor(1:458)
nfounders<-8
num<-1000
for ( i in 1:num){
  s0<-Sys.time()
  cat("Start:",i,"\n")
  datafile<-paste(pathdata,"simuldata.",i-1,".csv",sep="") 
  y0<-read.csv(file=datafile,header=TRUE)$Phe  
  pheno.data<-data.frame(y=y0,id=pheno)
  parms<-mpwgaim.cal(pheno.data=pheno.data,mpInterval=mpInterval,nfounders=nfounders)
  parmfile<-paste(pathout,"mpwgaim.",i-1,".blupp.csv",sep="")
  write.csv(parms,file=parmfile) 
}

##-2
#permutation-simulating in the null models
rm(list=ls())
library(asreml)
library(mpwgaim)
source("mpwgaim.R")
source("simH2.R")
#
pathout<-"./Data/Simulation/permutation/mpwgaim/"
dir.create(path=pathout,recursive=TRUE,showWarnings=FALSE)
load(file="mpInterval.RData")
load(file="chrs.magic.RData")
load(file="chrs.kkship.RData")
GG<-list(gen[[1]],kk.eigen[[1]]) 
nfounders<-8
num<-1000
wald<-NULL
for (i in 1:1000){
  X0<-simulH2(GG=GG,addqtl=FALSE) 
  X1<-X0$Parms
  X2<-X0$Datas
  indi<-nrow(X2)
  pheno.data<-data.frame(y=X2$Phe,id=factor(1:458))
  parms<-mpwgaim.cal(pheno.data=pheno.data,mpInterval=mpInterval)
  wald.one<-as.numeric(as.character(parms$wald))
  wald<-c(wald,max(wald.one))
  parmfile<-paste(pathout,"mp.",i-1,".blupp.csv",sep="")
  write.csv(parms,file=parmfile,row.names=FALSE)  
}
waldfile<-paste(pathout,"maxparms.wald.permut.csv",sep="")
write.csv(wald,file=waldfile,row.names=FALSE) 


##-3
##Ara
rm(list=ls())
library(asreml)
library(mpwgaim)
source("mpwgaim.R")
#
phe.name<-c("bolt.to.flower","growth.rate")
nphe<-length(phe.name)
Phe<-read.csv(file="Ara.phe.impute.csv",header=TRUE)
indi<-nrow(Phe)
pheno<-factor(1:indi)
nfounders<-19   
load(file="mpInterval.Ara.RData")
names(mpInterval$geno)<-"chr1"
#
for ( i in 1:2){
  pathout<-paste("./Data/Ara/",phe.name[i],"/",sep="")
  ii<-i+1
  y0<-phe[,ii] 
  pheno.data<-data.frame(y=y0,id=pheno)
  parms<-mpwgaim.cal(pheno.data,mpInterval,nfounders=nfounders)
  parmfile<-paste(pathout,"mp.blupp.csv",sep="")
  write.csv(parms,file=parmfile,row.names=FALSE)
}





##-4
rm(list=ls())
library(asreml)
library(mpwgaim)
source("mpwgaim.R")
#
phe.name<-c("PMN","CXCL1")
nphe<-length(phe.name)
Phe<-read.csv(file="CC.phe.impute.csv",header=TRUE)
Phe<-data.frame(subject.id=Phe[,1],PMN=log(Phe[,2]),CXCL1=log(Phe[,3]))
indi<-nrow(Phe)
#
load(file="mpInterval.CC.RData")
pheno<-factor(1:indi)
nfounders<-8
#
for ( i in 1:2){
  #
  pathout<-paste("./Data/CC/",phe.name[i],"/",sep="")
  pheno.data<-data.frame(y=y0,id=pheno)
  parms<-mpwgaim.cal(pheno.data,mpInterval,nfounders=nfounders)
  parmfile<-paste(pathout,"mp.blupp.csv",sep="")
  write.csv(parms,file=parmfile,row.names=FALSE)
}
