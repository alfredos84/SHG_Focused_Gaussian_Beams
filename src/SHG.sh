#!/bin/bash

# This file contains a set of instructions to run the main file SHG.cpp.
# Use this file to perform simulations in bulk. This means, for example, 
# systematically varying the input power, the beam weist, etc.

TIMEFORMAT='It took %0R seconds.' 
time { 

clear                            # Clear screen
rm *.dat                         # This removes all .dat files in the current folder. Comment this line for safe.
# rm *.txt                       # This removes all .txt files in the current folder. Comment this line for safe. 
rm SHG                           # This removes a previuos executable file (if it exist)


########################################################################
# -----Compilation-----
# Notice that for computing thermal profile user must use -DTHERMAL preprocessor variable in the compilation.
 
# g++ SHG.cpp -I /usr/local/include -L /usr/local/include/fftwf3 -lfftw3 -lfftw3f -lm -o SHG
# FOLDERSIM="Simulations"

g++ SHG.cpp -DTHERMAL -I /usr/local/include -L /usr/local/include/fftwf3 -lfftw3 -lfftw3f -lm -o SHG
FOLDERSIM="Simulations"


########################################################################
# The variables defined below (ARGX) will be passed as arguments to the main file 
# cuSHG.cu on each execution via the argv[X] instruction.

ARG1=(30)     	# Pump power                (ARG1)
ARG2=(29)      	# Beam waist (um)           (ARG2)
ARG3=(55.66) 	# 1st peltier temperature   (ARG3)
ARG4=(55.66) 	# 2nd peltier temperature   (ARG4)

SOLS=1			# Save Only Last Slice (FALSE=0,TRUE=1)
STMP=1			# Save the entire temperature tensor (FALSE=0,TRUE=1)

x=$(awk 'BEGIN{for(i=20; i<=40; i=i+1)print i}')
# n=$(awk 'BEGIN{for(i=1; i<=31; i=i+3)print i}')
# t1=$(awk 'BEGIN{for(i=54; i<=59; i=i+0.2)print i}')
# Each for-loop span over one or more values defined in the previous arguments. 
# 		for N in $x
# 		do
# 
# touch "CPUTime_NZ_100_N_128_T.txt"

# for (( p=0; p<${#ARG1[@]}; p++ ))
# do
# for N in $n
# do
# for (( x=0; x<${#ARG2[@]}; x++ ))
# do
# 	for WA in $x
# 	do
		# for (( t1=0; t1<${#ARG3[@]}; t1++ ))
		# do
# 		for TOVEN1 in $t1
# 		do

			# N=${ARG1[$p]}
			# printf "\nPower			= ${N} W\n" 
			
			# WA=${ARG2[$x]}
			# printf "\nBeam waist	= ${WA} \n" 
			
			# TOVEN1=${ARG3[$t1]}
			# printf "\nT peltier 1		= ${TOVEN1} ºC\n"
			
			# TOVEN2=${TOVEN1}
			# printf "\nT peltier 2	= ${TOVEN2} ºC\n"
			
			printf "\nMaking directory...\n"
			# FOLDER="sPPLT_P_${N}_W_T1_${TOVEN1}_WAIST_${WA}"
			# FILE="sPPLT_P_${N}_W_T1_${TOVEN1}_WAIST_${WA}.txt"
			FOLDER="sPPLT_P_${ARG1}_W_T1_${ARG3}_WAIST_${ARG2}"
			FILE="sPPLT_P_${ARG1}_W_T1_${ARG3}_WAIST_${ARG2}.txt"
			
			printf "Bash execution and writing output file...\n\n"
			# ./SHG $N $WA $TOVEN1 $TOVEN2 $SOLS $STMP| tee -a $FILE
			./SHG $ARG1 $ARG2 $ARG3 $ARG4 $SOLS $STMP| tee -a $FILE

			printf "Bash finished!!\n\n" 
			mkdir $FOLDER
			mv *.dat $FOLDER"/"
			mv FILE.txt $FOLDER"/"

# 		done
# 	done
# done


if [ -d "$FOLDERSIM" ]; then
	echo "Moving simulations in ${FOLDERSIM}..."
	mv sPPLT* $FOLDERSIM"/" 
else

	mkdir $FOLDERSIM
	echo "Creating and moving simulations in ${FOLDERSIM}..."
	mv sPPLT* $FOLDERSIM"/" 
fi

}