##write shell script. 

data<-data.frame(Dres=rep(0,1000),Rwt=0, DDR=0,Pfix=0,Psoft2=0)

i=0

for (Rwt in c(2, 4)){
	RwtChar=substr(as.character(Rwt+0.00001),1,4)
	for (Dres in seq(0,0.2,by=0.08)){
		Rres= Rwt*(1-Dres)
		RwtDrug = Rres*(1-0.05)
		DDR = 1-RwtDrug/Rwt
		#choose DDR so that a_b is same all the time!!
 			
			i=i+1
			print(paste("Dres",Dres,"; DDR",DDR))
			
			outputfile=	paste("m_seed7mut_on0Rwt",RwtChar,".csv",sep="")	
			system(paste("rm",outputfile))	
			
			filetowrite="ShellScriptFixedN.sh"
			write("#!/bin/bash",file=filetowrite)
			write(paste("seed=7\nnRuns=500\nRwt=",Rwt,"\nDres=",Dres,"\nDDR=",DDR,"\nmu=0.000005\nmu_after_on=0\nmig=0\nKmain=10000\nKrefu=0\nV=1\n",sep=""),file=filetowrite,append=TRUE)
			
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
			
			system("chmod +x ./ShellScriptFixedN.sh")
			system("./ShellScriptFixedN.sh")	
			read.csv(outputfile,header=FALSE,sep="\t")->X
			data$Rwt[i]=Rwt
			data$Dres[i]=Dres
			data$DDR[i]=DDR
			data$Pfix[i]=X[,which(X=="Fixprob")+1]
			data$Psoft2[i]=1-X[,which(X=="BenAlleHomozygosity")+1]
			data$PsoftRel[i]=data$Psoft2[i]/data$Pfix[i]
			
		}}


data<-data[data$Rwt>0,]
library(lattice)

png("ResultsPsoftRel.png")
levelplot(PsoftRel~Rwt*DDR, data, cuts=8,main="relative prob multiple origin soft sweep")
dev.off()
png("ResultsPfix.png")
levelplot(Pfix~Rwt*DDR, data, cuts=8,main="prob fix")
dev.off()


