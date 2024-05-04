#!/bin/bash

#Paths for the compound:

COMPOUND="Rh1"
PATH_DATA="media/sf_Processing/$COMPOUND"

#Paths of scripts and programs:

PATH_SCRIPTS="media/sf_Processing/ESRF_Scripts/Scripts"
XTAL_DIR="home/piotr/xtal"
LAUE_DIR="home/piotr/Desktop/laueutil/epd-7.0-1-rh5-x86_64/bin"

#remember to have: laser_01.inp, dark_01.inp, .mda, _mono.h5 files ready in a subdirectory named as the processed compound in the main scripts directory (e.g. /SCRIPTS/CuCu/laser_01.inp etc.)

##
## Advanced settings:
##

#Paths of files:

PATH_INP_DARK="${PATH_SCRIPTS}/${COMPOUND}/dark_01.inp"
PATH_INP_LASER="${PATH_SCRIPTS}/${COMPOUND}/laser_01.inp"
PATH_MONO="${PATH_SCRIPTS}/${COMPOUND}/${COMPOUND}_mono.h5"
PATH_BEAM="${PATH_SCRIPTS}/${COMPOUND}/${COMPOUND}.mda"

PATH_OUT="${PATH_DATA}_Processed"
PATH_RESULTS="${PATH_DATA}_Results"

#Specify which full scan directories you wish to omit (i.e. when processing for some data is already done or unimportant).
#For example OMIT_FULL=(1 2 5 7 8) to omit COMP_01, COMP_02, COMP_05, COMP_07, COMP_08 directories in COMP/ directory
OMIT_FULL=(1 2 3 4 5 6 7)

#Specify which subdirectories you wish to omit. 
#For example OMIT_SUB=(3 "dark_01" 8 "laser_01") to omit COMP_03/dark_01, COMP_08/laser_01 directories in COMP/ directory
OMIT_SUB=()

#CAUTION!! Don't use manual exclusions OMIT_FULL or OMIT_SUB if your data is incomplete (e.g. there are COMPOUND_01 and COMPOUND_03 main directories, but no COMPOUND_02 main directory). This exception only applies to main directories, incomplete subdirectories are fine. If you wish to exclude directories anyway, simply count them alphabetically to get their "array" assignment.

#Use options below if you intend to only process dark_xx frames now or if you want to process everything else (without dark_xx). Leave both at 0 otherwise
JUST_DARK=0
ALL_BUT_DARK=0

#Use option below to omit directories which have an underscore ( _ ) sign in the beggining
OMIT_UNDERSCORED=1

#Use option below if you wish to delete exported _frm files after -d subroutine
CLEAN=0

#VERY IMPORTANT:::: (if you are matching pairs for dark sets)
#Make sure that the fake lu_match_pairs.py file in the /SCRIPTS looks like this:
#line101::		ecf = laue_util.expt_filter.ExptClusteringFilter(use_idx_expt=False, min_I=50000)
#line219::   	f= open("PATH_SIZE","w")
#Abovementioned lines can only get corrupt if the match pairs processing routine has been aborted by user (C^)

######################################################################
######################################################################
######################################################################

#Some required repathing

PROCESS_CHOICE=${1?Error: Please define your processing procedure or use ./processing.sh h for help}
MATCH_PAIRS_DIR="${PATH_SCRIPTS}"

#Help regarding choice

if [ "$PROCESS_CHOICE" == "h" ]
then
	echo To specify what would you like to do:
	echo  
	echo c for creating directories and .inp files
	echo e for exporting
	echo s for signal processing
	echo i for integrating
	echo d for dumping to LaueUtil files AND DELETING EXPORTED FILES
	echo a for exporting, signalprocessing, integrating and dumping to LU, taking exemplary frames AND DELETING EXPORTED FILES
	echo p for plotting integrated powerscan data
	echo m for matching pairs in dark sets
	echo o for orientation matrix refinement
	echo hkl for hkl assignment and plotting euler angles
	echo r for hkl index to ratio conversion
	echo sortav for sortav processing
	echo b to save masks
	echo u for model refinement
exit
fi

#Getting the full name of process
PROCESS_NAME=""
case $PROCESS_CHOICE in
	c)PROCESS_NAME="Creating"
	;;
	e)PROCESS_NAME="Exporting"
	;;
	s)PROCESS_NAME="Signal"
	;;
	i)PROCESS_NAME="Integrating"
	;;
	d)PROCESS_NAME="DumpingAndDeleting"
	;;
	a)PROCESS_NAME="All"
	;;
	p)PROCESS_NAME="Powerscans"
	;;
	m)PROCESS_NAME="Matching"
	;;
	o)PROCESS_NAME="OrientationMatrix"
	;;
	hkl)PROCESS_NAME="hklAndPlot"
	;;
	r)PROCESS_NAME="Ratios"
	;;
	sortav)PROCESS_NAME="Sortav"
	;;	
	b)PROCESS_NAME="Masks"
	;;
	u)PROCESS_NAME="Refinement"
	;;
esac

#Dismantling the need to use original data if user has already produced local .inp files
T_TOT=0
if [ "$PROCESS_CHOICE" == "c" ]
	then echo Using raw .mccd data to produce .inp files
	else 
		PATH_DATA=$PATH_OUT
		echo Using .inp files and their local directories for futher processing
fi

#Mapping the main directory into an array

cd /$PATH_DATA/
ARRAY_MAIN=($(ls -d */))

#Echoing the number of viable main directories

echo ${#ARRAY_MAIN[@]} data subdirectories found

#Get into each subdirectory (loop starts here)

for ((j=0; j<${#ARRAY_MAIN[@]}; j++))
do

#Omiting main directories excluded manually
	
	OMIT_CHECK=0
	for ((g=0; g<${#OMIT_FULL[@]}; g++))
	do
		if (( "$j + 1" == "${OMIT_FULL[$g]}" ))
			then OMIT_CHECK=1
		fi
	done
	if [ "$OMIT_CHECK" == 1 ]
		then
			echo ...Omiting ${ARRAY_MAIN[$j]}
			continue
	fi
	
#Map subdirectories into an array
	
	cd /$PATH_DATA/${ARRAY_MAIN[$j]}
	ARRAY=($(ls -d */))
	
#Create an array for potential powerscans
	
	ARRAY_PS=()

#Get into each subdirectory (loop starts here)

	for ((k=0 ; k<${#ARRAY[@]}; k++))
	do
		cd ${ARRAY[$k]}
		FOLDER=${PWD##*/}
		
#Omiting subdirectories according to possible ALL_BUT_DARK or JUST_DARK preferences
		
		if [ "dark" == "${FOLDER:0:4}" ]
			then
			if [ "$ALL_BUT_DARK" == "1" ]
				then 
					cd ..
					continue
			fi
			elif [ "$JUST_DARK" == "1" ] 
				then
					cd ..				
					continue	
		fi
		
#Omiting subdirectories with an underscore as the first sign
		
		if [ "_" == "${FOLDER:0:1}" ]
			then
			if [ "$OMIT_UNDERSCORED" == "1" ]
				then 
					cd ..
					continue
			fi
		fi
		
#Omiting subdirectories excluded manually
		
		OMIT_CHECK=0
		for ((g=0; g<${#OMIT_SUB[@]}; g=g+2))
		do
			if (( "$j + 1" == "${OMIT_SUB[$g]}" ))
				then ((OMIT_CHECK++))
			fi
		done
		for ((g=1; g<${#OMIT_SUB[@]}; g=g+2))
		do
			if [ "$FOLDER" == "${OMIT_SUB[$g]}" ]
				then ((OMIT_CHECK++))
			fi
		done
		if [ "$OMIT_CHECK" == 2 ]
			then 
				echo ...Omiting ${ARRAY_MAIN[$j]}$FOLDER
				cd ..
				continue
		fi
		
#Echoing current directory/subdirectory/ and taking timestamp
		
		echo Processing ${ARRAY_MAIN[$j]}$FOLDER
		T_0=$(date +%s)
		
#Processing the subdirectory
		
		case $PROCESS_CHOICE in
			c)
				mkdir -p /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER
				if [ "dark" == "${FOLDER:0:4}" ]
					then 
						yes | cp  -rf /$PATH_INP_DARK /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER
						mv /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/dark_01.inp /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp 2>/dev/null
					else 
						yes | cp -rf /$PATH_INP_LASER /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER 
						mv /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/laser_01.inp /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp 2>/dev/null
				fi
				sed -i s@PATH_INPUT@$PATH_DATA/${ARRAY_MAIN[$j]}$FOLDER@g /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp
				sed -i s@PATH_LOG@$PATH_DATA/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.log@g /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp
				sed -i s@PATH_OUTPUT@$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER@g /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp
				sed -i s@NAME_LOOK@$FOLDER@g /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp
				if [ ! -f /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/number_of_pairs.txt ]
					then 
						FRAMES=($(ls -l /$PATH_DATA/${ARRAY_MAIN[$j]}${FOLDER}/*001.mccd | wc -l)) # ty J.G.
						echo "$((FRAMES/2))" >> /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/number_of_pairs.txt
				fi
				;;
			s)
				/${XTAL_DIR}/laueproc -r /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp -s
				;;
			e)
				/${XTAL_DIR}/laueproc -r /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp -e
				;;
			i)
				/${XTAL_DIR}/laueproc -r /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp -i			
				;;				
			d)
				/${XTAL_DIR}/laueproc -r /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp -d
				if [ $CLEAN -gt 0 ]
					then
						rm /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}__frm.h5
				fi
				;;
			a)
				/${XTAL_DIR}/laueproc -r /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp -a
				/${XTAL_DIR}/laueproc -r /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp -b -f 0 -l 4 -p "filtered"	
				rm /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}__frm.h5
				;;				
			p)	if [ "${FOLDER:0:9}" == "powerscan" ]
					then
						ARRAY_PS+=("$FOLDER")
					else 
						echo No powerscan detected in this subdirectory
				fi
				if (( "$k + 1" == "${#ARRAY[@]}" ))
					then 
						for ((t=0; t<${#ARRAY_PS[@]}-1; t++))
						do
							cp /$PATH_SCRIPTS/templates/ps_match_pairs.py /$PATH_OUT/${ARRAY_MAIN[$j]}
							sed -i s/POWERSCAN_ONE/${ARRAY_PS[$t]}/g /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
							sed -i -i s@PATH_OUT@$PATH_OUT/${ARRAY_MAIN[$j]}@g /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
							((t++))
							sed -i s/POWERSCAN_TWO/${ARRAY_PS[$t]}/g /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
							sed -i -i s@PATH_OUT@$PATH_OUT/${ARRAY_MAIN[$j]}@g /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
							((t--))
							echo now will plot ${ARRAY_PS[$t]} and one above
							/$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
							rm /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
						done
						
						for ((t=2; t<${#ARRAY_PS[@]}; t++))
						do
							cp /$PATH_SCRIPTS/templates/ps_match_pairs.py /$PATH_OUT/${ARRAY_MAIN[$j]}
							sed -i s/POWERSCAN_ONE/${ARRAY_PS[0]}/g /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
							sed -i -i s@PATH_OUT@$PATH_OUT/${ARRAY_MAIN[$j]}@g /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
							sed -i s/POWERSCAN_TWO/${ARRAY_PS[$t]}/g /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py						
							sed -i -i s@PATH_OUT@$PATH_OUT/${ARRAY_MAIN[$j]}@g /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
							echo now will plot ${ARRAY_PS[0]} and ${ARRAY_PS[$t]}
							/$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
							rm /$PATH_OUT/${ARRAY_MAIN[$j]}match_pairs.py
						done						
				fi
				;;
			m)
				cp /$PATH_SCRIPTS/templates/lu__match_pairs_nogui_backup_.py /$PATH_SCRIPTS/templates/lu__match_pairs_nogui.py
				cp /$PATH_SCRIPTS/templates/lu__match_pairs_nogui_BROKEN_ORIGINAL.py /$PATH_SCRIPTS/templates/lu__match_pairs_nogui_BROKEN.py
				if [ "${FOLDER:0:4}" == "dark" ]
					then 
						echo "${ARRAY_MAIN[$j]}${FOLDER} ${PROCESS_NAME}" >> /$PATH_OUT/matchpairs_log.csv
						sed -i s@PATH_SIZE@/$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/size.txt@g /$PATH_SCRIPTS/templates/lu__match_pairs_nogui_BROKEN.py
						EPSI=150
						INTENSITY=45000
						TRIGGER=0
						while [ ${EPSI#-} -gt 15 ]
						do
							sed -i s/I=50000/I=$INTENSITY/g /$PATH_SCRIPTS/templates/lu__match_pairs_nogui_BROKEN.py							
							timeout 5m /$LAUE_DIR/ipython -wthread /$PATH_SCRIPTS/templates/lu__match_pairs_nogui_BROKEN.py  /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}__expt.h5 /$PATH_MONO
							sleep 10
							read -r SIZE < /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/size.txt
							echo "${INTENSITY} intensity threshold gave ${SIZE} peaks" >> /$PATH_OUT/matchpairs_log.csv
							echo ${INTENSITY} intensity threshold gave ${SIZE} peaks
							EPSI=$((SIZE-100))						
							sed -i s/I=$INTENSITY/I=50000/g /$PATH_SCRIPTS/templates/lu__match_pairs_nogui_BROKEN.py
							rm /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/size.txt
							if [ ${EPSI#-} -lt 15 ]
								then break
							fi
							if [ $INTENSITY -gt 800000 ]
								then break
							fi
							if [ $INTENSITY -le 10000 ]
								then break
							fi
							if [ $EPSI -gt 0 ]
								then 
									if [ $TRIGGER -le 0 ] 
										then
											((INTENSITY=INTENSITY+5000))
											TRIGGER=-1
										else break
									fi
								else 
									if [ $TRIGGER -ge 0 ] 
										then
											((INTENSITY=INTENSITY-5000))
											TRIGGER=1
										else break
									fi
							fi
						done
						sed -i s/I=70000/I=$INTENSITY/g /$PATH_SCRIPTS/templates/lu__match_pairs_nogui.py						
						timeout 30m /$LAUE_DIR/ipython -wthread /$PATH_SCRIPTS/templates/lu__match_pairs_nogui.py  /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}__expt.h5 /$PATH_MONO
						sleep 10						
						sed -i s/I=$INTENSITY/I=70000/g /$PATH_SCRIPTS/templates/lu__match_pairs_nogui.py
					else continue
				fi
				sed -i s@/$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/size.txt@PATH_SIZE@g /$PATH_SCRIPTS/templates/lu__match_pairs_nogui_BROKEN.py
				;;			
			o)
				if [ "${FOLDER:0:5}" == "laser" ]
					then 
						timeout 30m /$LAUE_DIR/ipython -wthread /$PATH_SCRIPTS/templates/lu__refine_orientmatrix_nogui.py /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}__expt.h5 /$PATH_MONO /$PATH_BEAM /$PATH_OUT/${ARRAY_MAIN[$j]}dark_0*/xxx.h5
						sleep 5
				fi
				;;
			hkl)
				if [ "${FOLDER:0:5}" == "laser" ]
					then 
						timeout 30m yes | /$LAUE_DIR/ipython -wthread /$PATH_SCRIPTS/templates/lu__hkl_assign_monolc_modified.py /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}__expt.h5 /$PATH_MONO /$PATH_BEAM /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}__expt_rot.h5
						sleep 5
						timeout 30m /$LAUE_DIR/ipython -wthread /$PATH_SCRIPTS/templates/lu__view_matrix_parameters_nogui.py /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}__expt_rot.h5
						sleep 5
						mv /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/euler_plot.png /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${ARRAY_MAIN[$j]:: -1}_${FOLDER}_euler_plot.png 2>/dev/null
						mkdir -p /$PATH_RESULTS
						cp /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${ARRAY_MAIN[$j]:: -1}_${FOLDER}_euler_plot.png /$PATH_RESULTS/${ARRAY_MAIN[$j]:: -1}_${FOLDER}_euler_plot.png
				
				fi
				;;
			r)
				if [ "${FOLDER:0:5}" == "laser" ]
					then
						cp /$PATH_SCRIPTS/templates/lu_indexing2ratios_ORIGINAL.py /$PATH_SCRIPTS/templates/lu_indexing2ratios.py				
						PAIRS=0
						read -r PAIRS < /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/number_of_pairs.txt
						sed -i s/LASER/${ARRAY_MAIN[$j]:: -1}_${FOLDER}/g /$PATH_SCRIPTS/templates/lu_indexing2ratios.py
						sed -i s/NUM_PAIRS/$PAIRS/g /$PATH_SCRIPTS/templates/lu_indexing2ratios.py
						yes | /$LAUE_DIR/ipython -wthread /$PATH_SCRIPTS/templates/lu_indexing2ratios.py /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/yyy.h5
						sleep 5
				fi
				;;
			sortav)
				if [ "${FOLDER:0:5}" == "laser" ]
					then
						cp /$PATH_SCRIPTS/templates/sortav.inp /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/sortav.inp
						sed -i s@PATH_HKL@${ARRAY_MAIN[$j]:: -1}_${FOLDER}__ratios.hkl@g /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/sortav.inp
						sed -i s/OUTPUT/${ARRAY_MAIN[$j]:: -1}_${FOLDER}/g /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/sortav.inp
						cp /$PATH_SCRIPTS/templates/f_sortav.out /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/f_sortav.out
						/$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/f_sortav.out
						
						NUM_LINES=($(wc -l /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${ARRAY_MAIN[$j]:: -1}_${FOLDER}.sortav | awk '{ print $1 }'))
						echo $((NUM_LINES-7)) reflections
						for ((u=8 ; u<$((NUM_LINES+1)); u++))
						do
							OUT_STRING=""
							CURRENT_LINE="$(cat /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${ARRAY_MAIN[$j]:: -1}_${FOLDER}.sortav | head -$u | tail -1)"
							arr=($CURRENT_LINE)

							ITEM=${arr[0]}
							OUT_STRING="$OUT_STRING $ITEM"
							ITEM=${arr[1]}
							OUT_STRING="$OUT_STRING $ITEM"
							ITEM=${arr[2]}
							OUT_STRING="$OUT_STRING $ITEM"

							ITEM=${arr[3]}
							NUMBER=($(echo $ITEM -n | awk -F"E" 'BEGIN{OFMT="%10.9f"} {print $1 * (10 ^ $2)}'))
							OUT_STRING=" $OUT_STRING $NUMBER"

							ITEM=${arr[4]}
							NUMBER=($(echo $ITEM -n | awk -F"E" 'BEGIN{OFMT="%10.9f"} {print $1 * (10 ^ $2)}'))
							OUT_STRING=" $OUT_STRING $NUMBER"
							mkdir -p /$PATH_RESULTS
							echo $OUT_STRING >> /$PATH_RESULTS/${ARRAY_MAIN[$j]:: -1}_${FOLDER}_final.hkl
						done
						cp /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/sortav.lst /$PATH_RESULTS/${ARRAY_MAIN[$j]:: -1}_${FOLDER}_final.lst
				fi
				;;	
			b)
				/${XTAL_DIR}/laueproc -r /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp -b -f 0 -l 4 -p "filtered"				
				;;	
			u)
				/${XTAL_DIR}/laueproc -r /$PATH_OUT/${ARRAY_MAIN[$j]}$FOLDER/${FOLDER}.inp -u				
				;;					
			*)  echo Not compatible
				;;
		esac
		cd ..
		T_E=$(date +%s)
		T_TOT=$((T_TOT+T_E-T_0))
		echo Processed in $((T_E-T_0)) seconds
		DATE=`date '+%Y-%m-%d %H:%M:%S'`
		echo "${DATE},${ARRAY_MAIN[$j]}${FOLDER},${PROCESS_NAME},$((T_E-T_0))" >> /$PATH_OUT/log.csv		
	done
done
echo Processing successful. Total time was $T_TOT seconds
