#!/bin/bash

input=./scripts
path=Path_OffResonantData_all.txt
madey=ORD
hbk_path=./histos
filey=$madey-script
sush=./$madey-submit.sh
echo "#!/bin/bash" > $sush
upcount=10
count=0
check=$upcount

while read line
do
    if [ $check -eq $upcount ]
    then
	sfile=$input
	
        ff=$sfile/$filey-$count.sh
	check=0
        echo "#!/bin/sh" >> $ff
	
	echo "###################" >> $ff
	echo "source /sw/belle/local/etc/bashrc_general" >> $ff

	echo "export USE_GRAND_REPROCESS_DATA=1" >> $ff
	echo "export BASF_USER_IF=basfsh.so" >> $ff
	echo "export BASF_USER_INIT=user_init.so" >> $ff
	echo "export BELLE_MESSAGE_LEVEL=ERROR" >> $ff

	echo "#basf << EOF >& ./log/$filey-$count.log" >> $ff
	echo "basf <<EOF" >> $ff
	echo "path create analysis" >> $ff
	echo "module register fix_mdst User_reco user_index" >> $ff
	echo "path add_module main fix_mdst  User_reco user_index" >> $ff
	echo "path add_condition main <:0:KILL" >> $ff
	echo "path add_condition analysis ==:0:KILL" >> $ff
	 
        echo "initialize" >> $ff
	echo "output open" >> $ff
	echo "histogram define  ./histos/$filey-$count.hbk" >> $ff

	echo "process_event $line 0" >> $ff

	check=` expr $check + 1 `
        echo "bsub -q l -o ./log/$filey-$count.log $ff" >> $sush
    else
	echo "process_event $line 0 " >> $ff
        check=`expr $check + 1 `
    fi
    
    if [ $check  -eq $upcount ]
    then
        echo "output close" >> $ff
	echo " " >> $ff
        echo "terminate "    >> $ff
        echo " "    >> $ff
        echo "EOF"  >> $ff
        echo " "    >> $ff
        echo "#exit 0"  >> $ff
        echo " "    >> $ff
	echo "h2root ./histos/$filey-$count.hbk ./root_files/$filey-$count.root"  >> $ff
        echo " "    >> $ff
	count=` expr $count + 1 `
    fi
    
    
done < $path   




echo "#output close" >> $ff
echo " " >> $ff
echo "terminate "    >> $ff
echo " "    >> $ff
echo "EOF"  >> $ff
echo " "    >> $ff
echo "exit 0"  >> $ff
#chmod u+x ./scripts/*.csh
chmod u+x ./scripts/*.sh
chmod u+x *.sh