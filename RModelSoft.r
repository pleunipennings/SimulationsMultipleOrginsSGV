##write shell script. 


data<-data.frame(Dres=rep(0,1000),Rwt=0, DDR=0,Pfix=0,Psoft=0)

Rwt=2; 
i=0

for (Rwt in seq(1.1,3,by=0.2)){
	RwtChar=substr(as.character(Rwt+0.00001),1,4)

#for (Dres in seq(0,Rwt,by=0.4)){
	for (Dres in 0.05){
		for (DDR in seq(Dres,1,by=0.1)){
			
			i=i+1
			print(paste("Dres",Dres,"; DDR",DDR))
			
			outputfile=	paste("m_seed7mut_on1Rwt",RwtChar,".csv",sep="")	
			system(paste("rm",outputfile))	
			
			filetowrite="ShellScript.sh"
			write("#!/bin/bash",file=filetowrite)
			write(paste("seed=7\nnRuns=2000\nRwt=",Rwt,"\nDres=",Dres,"\nDDR=",DDR,"\nmu=0.000005\nmu_after_on=1\nmig=0\nKmain=10000\nKrefu=0\nV=1\n",sep=""),file=filetowrite,append=TRUE)
			
#Dres is the cost of the mutation. Should vary from 0 to Rwt (2 in this case)
			write('echo "
				  $seed 
				  $nRuns
				  $Rwt
				  $Dres
				  $DDR
				  $mu
				  $mu_after_on
				  $mig
				  $Kmain
				  $Krefu
				  $V
				  " | ./HIVevolution',file=filetowrite,append=TRUE)
			
			system("chmod +x ./ShellScript.sh")
			system("./ShellScript.sh")	
			read.csv(outputfile,header=FALSE,sep="\t")->X
			data$Rwt[i]=Rwt
			data$Dres[i]=Dres
			data$DDR[i]=DDR
			data$Pfix[i]=X[,which(X=="Fixprob")+1]
			data$Psoft[i]=X[,which(X=="ProbSoft")+1]
			data$PsoftRel[i]=data$Psoft[i]/data$Pfix[i]
			
		}}}


data<-data[data$Rwt>0,]
library(lattice)

png("ResultsPsoftRel.png")
levelplot(PsoftRel~Rwt*DDR, data, cuts=8,main="relative prob multiple origin soft sweep")
dev.off()
png("ResultsPfix.png")
levelplot(Pfix~Rwt*DDR, data, cuts=8,main="prob fix")
dev.off()


if (FALSE){

for (mut_on in 0:1){

read.table(paste("m_seed7mut_on",mut_on,"Rwt2.00.csv",sep=""),stringsAsFactors=FALSE)->Data
#get column names by taking uneven columns of first row

names(Data)[seq(2,length(Data[1,]),by=2)]=as.character(Data[1,seq(1,length(Data[1,]),by=2)])
Data<-Data[,seq(2,length(Data[1,]),by=2)]

png(paste("Picture1_mut_on",mut_on,".png",sep=""))

plot(Data$DDR,Data$Fixprob,t="n",ylab="Pfix",xlab="Strength Drugs",ylim=c(0,1),xlim=c(0,1),main="Prob Fixation from SGV")

polygon(c(0,Data$DDR,1),c(0,Data$Fixprob,0),col=2,)
polygon(c(0,Data$DDR,1),c(0,Data$ProbSoft,0),col=4)

legend("topright", inset=.05, 
c("no adaptation","hard sweep","soft sweep"), fill=c(0,2,4), horiz=FALSE)
dev.off()
png(paste("Picture2_mut_on",mut_on,".png",sep=""))



plot(Data$DDR,Data$Fixprob,t="n",ylab="Pfix",xlab="Strength Drugs",ylim=c(0,1),xlim=c(0,1),main="Prob soft given fixation")

polygon(c(0.12,0.12,1,1),c(0,1,1,0),col=2)
polygon(c(.12,Data$DDR[3:length(Data$DDR)],1),c(0,Data$ProbSoft[3:length(Data$DDR)]/Data$Fixprob[3:length(Data$DDR)],0),col=4)

legend("topright", inset=.05, 
c("no adaptation","hard sweep","soft sweep"), fill=c("white",2,4), horiz=FALSE)

	dev.off()}

}