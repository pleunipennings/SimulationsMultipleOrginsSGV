##write shell script. 

setwd("/Users/pleunipennings/Documents/Research/HIV/LossAndRecoveryHIVDiversity/Simulations")
system("./make_HIV")	
setwd("/Users/pleunipennings/Documents/Research/HIV/LossAndRecoveryHIVDiversity/Simulations/SimResults")

data<-data.frame(sb=rep(0,1000),sd=0, mu=0,Pfix=0,Psoft2=0)

i=0
#i=length(data$mu)
seed = 7
mut_on=1
mu = 0.00002
WantedNRuns = 5000
Kmain = 10000
sbvalues = c(0.2, 0.1, 0.05, 0.02) 
sdvalues = c(0, 0.01, 0.02, 0.1, 0.2)
#sbvalues=0.2
#sdvalues=0

plot(0:1,0:1,xlim=c(0,1.5*max(sdvalues)),ylim=c(0,1),col=0)
points(0:1,c((mu*Kmain*2)/(mu*Kmain*2+1),(mu*Kmain*2)/(mu*Kmain*2+1)),t="l")

for (sb in sbvalues){
	for (sd in sdvalues){
			i=i+1
		seed = sample(1:10000,1)
			data$Ne[i] = Kmain 
			data$sd[i]=sd
			data$sb[i]=sb
			data$mu[i]=mu
			data$theta[i] = Kmain*mu*2
			data$alphad[i] = data$sd[i]*data$Ne[i]*2
			data$alphab[i] = data$sb[i]*data$Ne[i]*2
			data$predictedPfix[i] = 1- (1+(data$alphab[i]/(data$alphad[i]+1)))^-data$theta[i]
			data$predictedPsoft2[i] = data$theta[i]/(data$theta[i]+1)
			
			nRuns = floor(WantedNRuns/data$predictedPfix[i])
			
			print(paste("sd",sd,"; sb",sb,"; nRuns",nRuns))
			outputfile=	paste("m_seed",seed,"mut_on",mut_on,".csv",sep="")	
			system(paste("rm",outputfile))	
			
			filetowrite="../ShellScriptFixedN.sh"
			write("#!/bin/bash",file=filetowrite)
			write(paste("seed=",seed,"\nnRuns=",nRuns,"\nsd=",sd,"\nsb=",sb,"\nmu=",mu,"\nmu_after_on=",mut_on,"\nKmain=",Kmain,"\n",sep=""),file=filetowrite,append=TRUE)
			
#Dres is the cost of the mutation. Should vary from 0 to Rwt (2 in this case)
			write('echo "
				  $seed 
				  $nRuns
				  $sd
				  $sb
				  $mu
				  $mu_after_on
				  $Kmain
				  " | ../HIVevolution',file=filetowrite,append=TRUE)
			
			system("chmod +x ../ShellScriptFixedN.sh")
			system("../ShellScriptFixedN.sh")	
			read.csv(outputfile,header=FALSE,sep="\t")->X
			data$Pfix[i]=X[,which(X=="Fixprob")+1]
			data$Psoft2[i]=1-X[,which(X=="BenAlleHomozygosity")+1]
#data$PsoftRel[i]=data$Psoft2[i]/data$Pfix[i]
		
		print(data[i,])
		
		points(sd,data$Psoft2[i],pch=16,col=which(sbvalues==sb))

		}}


#data<-data[data$mu>0,]
#print(data[data$mu>0,])
png("Pmultiorigins_SGV.png")

plot(0:1,0:1,xlim=c(0,0.25),ylim=c(0,1),col=0,ylab="P_soft_2",xlab="cost of resistance allele")
points(0:1,c((mu*Kmain*2)/(mu*Kmain*2+1),(mu*Kmain*2)/(mu*Kmain*2+1)),t="l")

for (sb in sbvalues){
	for (sd in sdvalues){
		points(sd+(which(sbvalues==sb)*0.002),data$Psoft2[data$sb==sb&data$sd==sd],pch=16,col=which(sbvalues==sb),cex=2)}}
legend(0.2,0.8,c(sbvalues),col=1:4,pch=16,cex=1,title="sb values")
dev.off()

library(lattice)

#png("ResultsPsoftRel.png")
#levelplot(PsoftRel~sb*sd, data, cuts=8,main="relative prob multiple origin soft sweep")
#dev.off()
#png("ResultsPfix.png")
#levelplot(Pfix~sb*sd, data, cuts=8,main="prob fix")
#dev.off()


