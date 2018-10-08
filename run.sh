#!/bin/bash
N_RUNS=$(seq 1 100)
N_PARALLEL_TASKS=55

N=256
m=100
K=$(seq 15 5 55)
RUN_SESSION=0
PENALIZATION=$(seq 1 3)


for p in $PENALIZATION; do
	for k in $K; do
		for RUN in $N_RUNS; do
			## choose an algorithm:
                	./gorwl $N $k $m $RUN $p > "${k}_${RUN}.txt" &
			##./gorwl $N $k $m $RUN $p 25 > "${k}_${RUN}.txt" &
			##./gorwlist $N $k $m $RUN $p > "${k}_${RUN}.txt" & 
			##./gorwlist $N $k $m $RUN $p 25> "${k}_${RUN}.txt" &       
        	done
	done 



	for k in $K; do
		RSE=$(cat ${k}_*.txt | awk '{T=T+$1} END{print T/NR}')
		SUCCESS=$(cat ${k}_*.txt | awk '{T=T+$2} END{print T/NR}')
		ITER=$(cat ${k}_*.txt | awk '{T=T+$3} END{print T/NR}')
		TIME_SEC=$(cat ${k}_*.txt | awk '{T=T+$4} END{print T/NR}')
		echo "$k $RSE $SUCCESS $ITER $TIME_SEC" | tee -a RWL.dat
	done

	echo " " | tee -a RWL.dat
	echo " " | tee -a RWL.dat
	echo " " | tee -a RWL.dat

done 





