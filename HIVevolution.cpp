/***************************************************************************
 *   copyright Pleuni PENNINGS                                 *
 *   pennings@fas.harvard.edu                                  *
 *                                                                         *
 
 ***************************************************************************/


//Plan to make simple model with real growth rates, a poisson distribution of offspring numbers and a real population size, a carrying capacity. 
//At carrying capacity, the pop should behave like a WF model
//The environment changes the R0 of the wild type but not of the resistant type. Resistance comes with a cost. 
//The question is what is the fixation probability, probability of soft and hard sweep from SGVm new mutations or migration. 
//Initially, I will set, Rwt1>1 (initally the WT is fit), Rwt2=0 (the environmental change kills all WT), mig=0 (no migration). We should recover H&P2005 results for Psgv, but in addition determine number of origins (softness). 
//Next will look at Rwt2<1, but Rres2>1, which is more interesting. The prob of escape should be given by Uecker and H 2013. Not sure about Psoft and Phard | Pescape. 
//Note: RES = MUT

//Will change the code to have constant population size and see if the number of origins depends on ad or not. 


// System wide headers:
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sstream> 
#include <iostream> 
#include <fstream>
#include <vector> 
#include <math.h>
#include <string>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
using namespace std;

//Global variables
unsigned int repeat;
ofstream output;
int tspd = 1; //time steps per day I'll do discrete generations now, i.e. tspd = 1; 
double deathrate = 1./tspd; //will add second deathrate due to drugs later
int simlengthBEFORE = 1000; //to reach mut sel balance
int simlength=200+simlengthBEFORE; // total length of simulation

unsigned int popArray[100]; //Array that will hold the WT pop size and the Mut pop size. From 1 to 99 are diff origins of the mutant
unsigned int originGen[100]; //Array hols info of origin of  
unsigned int originType[100]; //Array hols info of origin: 1 = mutation
unsigned int popArrayReplace[100];//the newborns every generation 

//things to be filled during "popAdapt" make space for 5000 generations
int PopSizeWT[5000]; //population size of the wild type
int PopSizeMut[5000]; //pop size of the mutants
int PopTotal[5000]; //sum of WT and Mut popsizes

//descriptive statistics over a couple of runs
double MeanPopSizeWT[5000]; // calculate mean popsize in given generation averaged over all runs
double MeanPopSizeMut[5000]; // mutants mean pop size in all runs
double MeanNewMutants[5000]; //same but averaged over all runs
double MutPresent[5000]; //whether the mutant is present averaged over all runs
double MeanTfix = 0; //calculate time to fixation of mutant
double FixProb = 0; // to calculate the number of runs where the mutant has fixed and thereby the probability that it fixes. 
double PSoft = 0; 
double BenAlleHomozygosity = 0; // to get homozygosity at beneficial allele 1-BenAlleHomozygosity = Psoft2
double Psoft2 = 0; 

//Parameters to be set in "GetParameters"
unsigned int seed, nRuns; 
double Rwt; //fitness of WT in old envi
double Dres; //Addnl death rate of res types
double DDR; //addnl death rate in new envi
double mu;//mutation rate
unsigned int mu_after_on; //whether mutation continues after env change
double mig;//migration rate
unsigned int Kmain; // carrying capacity of main compartment
unsigned int Krefu; // carrying capacity of refugium
double V;

void getParameters(){
	cerr << "enter seed: "; cin >> seed; 
	//seed =1000*seed;
	cerr << "enter number of runs: "; cin >> nRuns; 
	cerr << "enter Rwt ";	cin >> Rwt; 
	cerr << "enter Dres ";	cin >> Dres; 
	cerr << "enter DrugDeathrate"; cin >> DDR; 
	cerr << "enter mu ";	cin >> mu; 
	cerr << "enter mu_after_on ";	cin >> mu_after_on; 
	cerr << "enter mig ";	cin >> mig; 
	cerr << "enter Kmain (carrying capacity) ";	cin >> Kmain;
	cerr << "enter Krefu (carrying capacity) ";	cin >> Krefu;
	cerr << "enter V, the variance in offspr number";cin >> V; 
}

void createOutputfile(){    
	char filename [100];
	sprintf (filename,"%s%u%s%u%s%.2f%s","m_seed", seed, "mut_on",mu_after_on,"Rwt", Rwt,".csv");
	cerr << "the name of the outputfile is: " << filename << "\n" ; 
	output.open(filename , ios::app);
}

void popAdapt(gsl_rng *rng){
	//popAdapt has two parts: initialize (to make a monomorphic population and seed the random number generator) and evolve (a loop that loops until the ancestral allele is lost). 
	
	//initialize
	gsl_rng_set (rng, (unsigned long)repeat+seed+1); 
	//set everything to zero for a new run
	for(int m=0; m < simlength ; m++){PopSizeWT[m]=0; PopSizeMut[m]=0; PopTotal[m]=0;}  
	PopSizeWT[0]=Kmain;PopTotal[0]=Kmain;
	popArray[0]=PopSizeWT[0]; 
	for (int i= 1; i<100; i++){popArray[i]=0;originGen[i]=0; originType[i]=0; }
	
	int wt_from_res = 0; int res_from_wt = 0;
	double ExpNumOffspring = 0; 
	
	//evolve before/after change of env
	for (int t=1; t<simlength;t++)  { //each loop is one generation, continues until end of simlength
		ExpNumOffspring = 1 + (Rwt-1)*(1-(double)(popArray[0]+PopSizeMut[t-1])/(double)Kmain);
		if (t<simlengthBEFORE){
			popArrayReplace[0]=gsl_ran_poisson(rng, (double)popArray[0]*ExpNumOffspring);}
		//1-DDR is the fraction that survives despite drugs 
		if (t>=simlengthBEFORE){
			popArrayReplace[0]=gsl_ran_poisson(rng, (double)popArray[0]*ExpNumOffspring*(1.-DDR));}
		//each of the mutants has poisson num offspring
		for (int i= 1; i<100; i++){
			popArrayReplace[i]=gsl_ran_poisson(rng, (double)popArray[i]*ExpNumOffspring*(1.-Dres));}
		
//		cerr << t << "\t"  << "ExpNU" << "\t" <<ExpNumOffspring << "\t"<< PopSizeWT[t-1]<<"\t"<< PopSizeMut[t-1]<<"\n";

		//replace old array
		for (int j = 0; j<100; j++){
			popArray[j]=popArrayReplace[j]; }
		//count num mutants
		int sumofMutants = 0;
		for (int j = 1; j<100; j++){
			sumofMutants+=popArray[j];}
		
		//mutation only before or before and after the env change
		if (t<simlengthBEFORE){
			wt_from_res=gsl_ran_poisson(rng, mu*sumofMutants);
			res_from_wt=gsl_ran_poisson(rng, mu*popArrayReplace[0]);}
		if (t>=simlengthBEFORE && mu_after_on==1){
			wt_from_res=gsl_ran_poisson(rng, mu*sumofMutants);
			res_from_wt=gsl_ran_poisson(rng, mu*popArrayReplace[0]);}
		if (t>=simlengthBEFORE && mu_after_on==0){
			wt_from_res = 0; res_from_wt = 0;}
		
		//add the mutants to the popArray
		popArray[0]=popArray[0]+wt_from_res; //wt mutants
		for (int i = 0; i<res_from_wt; i++){ // res mutants (more important!)
			//find empty spot in popArray
			int foundspot=0;
			for (int j = 1; foundspot==0; j++){
				if (popArray[j]==0){
					foundspot = 1; popArray[j]++;
					originGen[j] = t; originType[j] = 1; 
				}}}
	
/*		
		for (int j = 0; j<10; j++){
			cerr << popArray[j]<<"\t"; 
		}
		cerr << "\n"; 
		for (int j = 0; j<10; j++){
			cerr << originGen[j]<<"\t"; 
		}
		cerr << "\n"; 
		for (int j = 0; j<10; j++){
			cerr << originType[j]<<"\t"; 
		}
		cerr << "\n\n"; 
*/		
		//Add Migration
		
		//Bookkeeping
		
		PopSizeWT[t]=popArray[0]; 
		for (int j = 1; j<100; j++){
			PopSizeMut[t]+=popArray[j];}
		PopTotal[t]=PopSizeWT[t]+PopSizeMut[t]; 
		
		//bookkeeping for the run
		MeanPopSizeWT[t]+=PopSizeWT[t]; MeanPopSizeMut[t]+=PopSizeMut[t]; if (PopSizeMut[t]>0) MutPresent[t]+=1; // bookkeeping for average over all runs
//		cerr<<t << "\t"<<popArray[0]<<"\t"<<popArray[1]<<"\n";

}

//cerr<<popArray[0]<<"\t"<<PopSizeMut[simlength-1]<<"\n";

if (PopSizeMut[simlength-1]>PopSizeWT[simlength-1]&&PopSizeMut[simlength-1]>0){
	int xx; for (xx=1; PopSizeMut[xx]<PopSizeWT[xx]; xx++); MeanTfix = MeanTfix+xx; }//get time of fixation after start of treatment
if ((double)PopSizeMut[simlength-1]>(double)PopSizeWT[simlength-1]&&PopSizeMut[simlength-1]>0) {FixProb+=1;
//check whether soft sweep (= multiple origins)
	int numOrigins = 0;
	for (int j = 1; j<100; j++){if (popArray[j]>0){numOrigins++;}}
	if (numOrigins>1) PSoft+=1; 
	double BenAlleHomozygosityThisRun=0; 
	for (int j = 1; j<100; j++){BenAlleHomozygosityThisRun+=double(popArray[j])*double(popArray[j])/(double(PopSizeMut[simlength-1])*double(PopSizeMut[simlength-1]));}
	BenAlleHomozygosity += BenAlleHomozygosityThisRun; 
//cerr<<	BenAlleHomozygosityThisRun<<"\t"<<BenAlleHomozygosity<<"\n";
	
} //I assume a mutation will fix if it has more than 0.5 frequency 
}


int main (int argc, char *argv[]){
	getParameters();
	createOutputfile();
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);//choose the random number generator
//	PopSizeWT[0]=Kmain; PopTotal[0]=Kmain; 
	for (int i = 0; i<5000; i++){MeanPopSizeWT[i]=0; MeanPopSizeMut[i]=0; MeanNewMutants[i]=0; MutPresent[i]=0;}
//	MeanPopSizeWT[0]=Kmain; MeanPopSizeMut[0]=0;MeanNewMutants[0]=0; MutPresent[0]=0; 	
	for (repeat = 0; repeat < nRuns; repeat++){
//		cerr << "repeat"<<repeat << "\t"; 
		popAdapt(rng);}//each of these loops is an independent run
		
	//calculate & output means over all generations
//	output<<"t\t" <<"PopSizeWT\t"<<"PopSizeMut\t"<<"MutPresent\t"<<"MeanNewMuts\tMut\tDrug\n";
	for (int t=1; t<simlength;t++)  {
		MeanPopSizeWT[t]/=(double)nRuns; 
		MeanPopSizeMut[t]/=(double)nRuns; 
		MeanNewMutants[t]/=(double)nRuns;
		MutPresent[t]/=(double)nRuns;
//			output<<t<<"\t"<<MeanPopSizeWT[t]<<"\t"<<MeanPopSizeMut[t]<<"\t"<<MutPresent[t]<<"\t"<<"\n";	
	}
	FixProb/=(double)nRuns; //to calculate the probability that adaptation has happened
	PSoft/=(double)nRuns; //to calculate the probability that adaptation has happened
	MeanTfix/=((double)nRuns*FixProb);MeanTfix-=simlengthBEFORE;
	BenAlleHomozygosity/=((double)nRuns*FixProb);
	
	//write summary output
	output << "Nruns\t"<< nRuns; 
	output << "\tSeed\t"<< seed; 
	output << "\tRwt\t"<< Rwt; 
	output << "\tDres\t"<< Dres; 
	output << "\tDDR\t"<< DDR; 
	output <<"\tmu\t"<<mu<<"\tmig\t"<<mig; 
	output <<"\tFixprob\t" <<FixProb<<"\tProbSoft\t"<<PSoft<<"\tV\t"<<V; 
	output <<"\tNMut\t"<<MeanPopSizeMut[simlength-1]; 
	output << "\tKmain\t"<< Kmain<< "\tKrefu\t"<< Krefu ; 
	output <<"\t"<<"\tMeanTfix\t"<<MeanTfix; 
	output <<"\t"<<"\tBenAlleHomozygosity\t"<<BenAlleHomozygosity<<"\n"; 
		
	return 5;
}

