/***************************************************************************
 *   copyright Pleuni PENNINGS                                 *
 *   pennings@fas.harvard.edu                                  *
 *                                                                         *
 
 ***************************************************************************/


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
int simlengthBEFORE = 10000; //to reach mut sel balance
int simlength=1000+simlengthBEFORE; // total length of simulation

unsigned int popArray[100]; //Array that will hold the WT pop size and the Mut pop size. From 1 to 99 are diff origins of the mutant
double weightedArray[100]; //SEP09similar as popArray, but then multiplied with the relative fitness of WT and Mut

unsigned int originGen[100]; //Array holds info of origin of  
unsigned int originType[100]; //Array holds info of origin: 1 = mutation (migration is other option)

//things to be filled during "popAdapt" make space for 11000 generations
int PopSizeWT[11000]; //population size of the wild type
int PopSizeMut[11000]; //pop size of the mutants
double RelFitnessRes[11000]; //SEP09relative fitness of the drug resistant mutant depends on Drug

//descriptive statistics over a couple of runs
double MeanPopSizeWT[11000]; // calculate mean popsize in given generation averaged over all runs
double MeanPopSizeMut[11000]; // mutants mean pop size in all runs
double MeanNewMutants[11000]; //same but averaged over all runs
double MutPresent[11000]; //whether the mutant is present averaged over all runs
double MeanTfix = 0; //calculate time to fixation of mutant
double FixProb = 0; // to calculate the number of runs where the mutant has fixed and thereby the probability that it fixes. 
double PSoftTotal = 0; //whether multiple origins in the population
double BenAlleHomozygosity = 0; // to get homozygosity at beneficial allele 1-BenAlleHomozygosity = Psoft2
double Psoft2 = 0; 

//Parameters to be set in "GetParameters"
unsigned int seed, nRuns; 

double sd; //Cost of mutation
double sb; //benefit of mutation 
double mu;//mutation rate
unsigned int mu_after_on; //whether mutation continues after env change
unsigned int Kmain; // carrying capacity of main compartment

void getParameters(){
	cerr << "enter seed: "; cin >> seed; 
	//seed =1000*seed;
	cerr << "enter number of runs: "; cin >> nRuns; 
	cerr << "enter sd ";	cin >> sd; 
	cerr << "enter sb ";	cin >> sb; 
	cerr << "enter mu ";	cin >> mu; 
	cerr << "enter mu_after_on ";	cin >> mu_after_on; 
	cerr << "enter Kmain (carrying capacity) ";	cin >> Kmain;
}

void createOutputfile(){    
	char filename [100];
//	sprintf (filename,"%s%u%s%u%s%f%s%f%s","m_seed", seed, "mut_on",mu_after_on,"sd", sd,"sb", sb,".csv");
	sprintf (filename,"%s%u%s%u%s","m_seed", seed, "mut_on",mu_after_on,".csv");
	cerr << "the name of the outputfile is: " << filename << "\n" ; 
	output.open(filename , ios::app);
}

void prepareTimeline(){
	int m=0; 
	for(m=0; m < simlengthBEFORE ; m++) {RelFitnessRes[m]=1-sd;}
	for(m=simlengthBEFORE; m < simlength ; m++) {RelFitnessRes[m]=1+sb;}
}

void updateWeightedArray(unsigned int *popArray, double *weightedArray, int t){ 
	weightedArray[0]= popArray[0] *1.0; //the 0 popArray have fitness 1
	for (int i= 1; i<100; i++){
		weightedArray[i] = popArray[i]*(RelFitnessRes[t]); //all 1 popArray have fitness 1+sb or 1-sd depending on whether drugs are taken
	}}


void popAdapt(gsl_rng *rng){
	//popAdapt has two parts: initialize (to make a monomorphic population and seed the random number generator) and evolve (a loop that loops until the ancestral allele is lost). 
	
	//initialize
	gsl_rng_set (rng, (unsigned long)repeat+seed+1); 
	//set everything to zero for a new run
	for(int m=0; m < simlength ; m++){PopSizeWT[m]=0; PopSizeMut[m]=0; }
	//PopTotal will stay the same at Kmain! July 2013
	PopSizeWT[0]=Kmain;
	popArray[0]=PopSizeWT[0]; 
	
	for (int i= 1; i<100; i++){popArray[i]=0;originGen[i]=0; originType[i]=0; }
	
	int wt_from_res = 0; int res_from_wt = 0;
	
	//evolve before/after change of env
	//use multinomial sampling 
	
	updateWeightedArray(popArray, weightedArray,0);
	int lengthofthissimulation = 1; 
	
	for (int t=1; t<simlength && PopSizeMut[t-1]<0.9*Kmain;t++)  { //each loop is one generation, continues until end of simlength
		gsl_ran_multinomial (rng, 100, Kmain, weightedArray, popArray); //reproduction
		//100 is the length of the array		
		//count num mutants in the population
		int sumofMutants = 0;
		for (int j = 1; j<100; j++){
			sumofMutants+=popArray[j];}
		
		//mutation only before or before and after the env change
		//only mutations in one direction is easier
		wt_from_res = 0; 
		if (t<simlengthBEFORE){
			res_from_wt=gsl_ran_poisson(rng, mu*popArray[0]);}
		if (t>=simlengthBEFORE && mu_after_on==1){
			res_from_wt=gsl_ran_poisson(rng, mu*popArray[0]);}
		if (t>=simlengthBEFORE && mu_after_on==0){res_from_wt = 0;}
		
		//add the mutants to the popArray
		popArray[0]=popArray[0]+wt_from_res-res_from_wt; //wt mutants
		for (int i = 0; i<res_from_wt; i++){ // res mutants (more important!)
			//find empty spot in popArray
			int foundspot=0;
			for (int j = 1; foundspot==0; j++){
				if (popArray[j]==0){
					foundspot = 1; popArray[j]++;
					originGen[j] = t; originType[j] = 1; 
				}}}
		
		updateWeightedArray(popArray, weightedArray,t);
					
		////bookkeeping for the run
		PopSizeWT[t]=popArray[0]; 
		for (int j = 1; j<100; j++){
			PopSizeMut[t]+=popArray[j];}
//		if (t>9990){cerr << t << "\t"<<PopSizeMut[t]<<"\n"; }

		// bookkeeping for average over all runs
		MeanPopSizeWT[t]+=PopSizeWT[t]; MeanPopSizeMut[t]+=PopSizeMut[t]; if (PopSizeMut[t]>0) MutPresent[t]+=1; 
		lengthofthissimulation = t; 
	}

	if (PopSizeMut[lengthofthissimulation]>PopSizeWT[lengthofthissimulation]&&PopSizeMut[lengthofthissimulation]>0){
		int xx; for (xx=1; PopSizeMut[xx]<PopSizeWT[xx]; xx++); MeanTfix = MeanTfix+xx; }//get time of fixation after start of treatment
	if ((double)PopSizeMut[lengthofthissimulation]>(double)PopSizeWT[lengthofthissimulation]&&PopSizeMut[lengthofthissimulation]>0) {//I assume a mutation will fix if it has more than 0.5 frequency 
		for (int j = 1; j<100; j++){if (popArray[j]>0){	
//			cerr << "j  " << popArray[j] << "\t" << originGen[j] << "\n";
		}}
		
		int FromSGV=0; 
		for (int j = 1; j<100; j++){if (popArray[j]>0 && originGen[j]<simlengthBEFORE){FromSGV++; 
		}}
		if (FromSGV>0){ // only count if at least one mut from SGV present. 
			FixProb+=1;
			//check whether soft sweep (= multiple origins) if fixation happened
			int numOrigins = 0;
			for (int j = 1; j<100; j++){if (popArray[j]>0){numOrigins++;}}
//			cerr << "numOrigins  "<< numOrigins<< "\n";
			if (numOrigins>1) PSoftTotal+=1; 
			double BenAlleHomozygosityThisRun=0; 
			for (int j = 1; j<100; j++){BenAlleHomozygosityThisRun+=double(popArray[j])*double(popArray[j])/(double(PopSizeMut[lengthofthissimulation])*double(PopSizeMut[lengthofthissimulation]));}
			BenAlleHomozygosity += BenAlleHomozygosityThisRun; 
//			cerr << "BenAlleHomozygosity "<<BenAlleHomozygosity<<"\n"; 
			// check if from SGV
		}}}


int main (int argc, char *argv[]){
	getParameters();
	createOutputfile();
	prepareTimeline();//added july 13
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);//choose the random number generator
//	PopSizeWT[0]=Kmain; 
	for (int i = 0; i<11000; i++){MeanPopSizeWT[i]=0; MeanPopSizeMut[i]=0; MeanNewMutants[i]=0; MutPresent[i]=0;}
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
	PSoftTotal/=(double)nRuns; //to calculate the probability that adaptation has happened
	MeanTfix/=((double)nRuns*FixProb);MeanTfix-=simlengthBEFORE;
	BenAlleHomozygosity/=((double)nRuns*FixProb);
	
	//write summary output
	output << "Nruns\t"<< nRuns; 
	output << "\tSeed\t"<< seed; 
	output << "\tsd\t"<< sd; 
	output << "\tsb\t"<< sb; 
	output <<"\tmu\t"<<mu; 
	output <<"\tFixprob\t" <<FixProb<<"\tProbSoftTotal\t"<<PSoftTotal; 
	output <<"\tNMut\t"<<MeanPopSizeMut[simlength-1]; 
	output << "\tKmain\t"<< Kmain; 
	output <<"\t"<<"\tMeanTfix\t"<<MeanTfix; 
	output <<"\t"<<"\tBenAlleHomozygosity\t"<<BenAlleHomozygosity<<"\n"; 
		
	return 5;
}

