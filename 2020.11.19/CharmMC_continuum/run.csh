rm -f log/*
rm -f histos/*
rm -f idx_files/*
rm -f root_files/*

#!/bin/bash
filename='Path_CharmMC.txt'
while read line; do
	bsub -q l -o "log/$(basename ${line}).log" "./Run_MyMC.sh ${line}"
done < $filename

