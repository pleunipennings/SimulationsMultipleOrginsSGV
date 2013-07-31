#!/bin/bash
seed=7
nRuns=2000
Rwt=2.9
Dres=0.05
DDR=0.95
mu=0.000005
mu_after_on=1
mig=0
Kmain=10000
Krefu=0
V=1

echo "
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
				  " | ./HIVevolution
