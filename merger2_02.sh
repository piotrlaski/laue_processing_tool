#!/bin/bash
LC_NUMERIC=C LC_COLLATE=C
PATH_SORTAV="media/sf_E_DRIVE/OLD/Merged_KW"

cd /$PATH_SORTAV/
ARRAY_SORTAV=($(ls))
printf "%4s %3s %3s %14s %14s %14s %14s %11s %11s %4s\n" "h" "k" "l" "Ioff" "sig(Ioff)" "Ion" "sig(Ion)" "eta" "sig(eta)" "bn" >> merged.m91

for ((b=0 ; b<${#ARRAY_SORTAV[@]}; b++))
do
	NUM_LINES=($(wc -l ${ARRAY_SORTAV[$b]} | awk '{ print $1 }'))
	for ((u=1 ; u<$((NUM_LINES+1)); u++))
	do
		CURRENT_LINE="$(cat ${ARRAY_SORTAV[$b]} | head -$u | tail -1)"
		arr=($CURRENT_LINE)
		ETA=($(echo "scale=4;(${arr[3]}-1)" | bc | awk '{printf "%0.5f", $0}'))
		SIG=($(printf "%0.5f" ${arr[4]}))
		NUM=$((b+1))
		printf "%4s %3s %3s %14s %14s %14s %14s %11s %11s %4s\n" ${arr[0]} ${arr[1]} ${arr[2]} "10.0000" "6000.0000" "9.0000" "5.0000" ${ETA} ${SIG} $NUM >> merged.m91
	done
done
