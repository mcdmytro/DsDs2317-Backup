#!/bin/sh

if test $1
then
    fname=$1
else
    echo fname missing
    exit 1
fi
oname=$(basename ${fname})
oname=${oname%.*}

#####################################
source /sw/belle/local/etc/bashrc_general

export USE_GRAND_REPROCESS_DATA=1
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so

#export BELLE_MESSAGE_LEVEL=DDEBUG
export BELLE_MESSAGE_LEVEL=ERROR
#### case B
# process_url "http://bweb3/mdst.php?bl=caseB&skm=HadronBJ&ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&dv=zfserv&dt=on_resonance  
#unset BELLE_USE_TMP
# for((i=${srun}; i<=${erun}; i++)) do

basf <<EOF
# create analysis path
path create analysis
# register modules, User_reco is your module
module register fix_mdst User_reco user_index
path add_module main fix_mdst  User_reco user_index
# conditions (aka status)
path add_condition main <:0:KILL
path add_condition analysis ==:0:KILL


initialize
# define output histogram, folder histos_MC MUST be createds
# output open idx_files/${oname}.idx
histogram define  histos/${oname}.hbk



# process_url http://bweb3/mdst.php?ex=${exp}&rs=${srun}&re=${erun}&skm=HadronBorJ&dt=on_resonance&bl=caseB
# process_event /group/belle/bdata_b/dstprod/dat/e0000${exp}/HadronBJ/1205/on_resonance/00/e000021r000002-b20021205_1443.mdst 0

# process url, format is self_explanatory, you can get fixed string at bweb3.cc.kekjp
# process_url http://bweb3/montecarlo.php?ex=${exp}&rs=${srun}&re=${erun}&dt=on_resonance&bl=caseB&ty=Any&st=${stm}

process_event ${fname} 0
output close

terminate

EOF

#done

#for((i=${srun};i<=${erun};i++)) do
# convert histogram to ROOT ntuple
h2root histos/${oname}.hbk root_files/${oname}.root 

#remove old histogram
# rm -f histos_MC/${exp}.${srun}.${erun}.hbk

#done

