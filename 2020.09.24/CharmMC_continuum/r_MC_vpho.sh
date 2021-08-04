rm -f log/*
rm -f histos/*
rm -f idx_files/*
rm -f root_files/*

for file in ../../MC/gsim/mdst/evtgen_exp_*_Y4STo3D*; do
	bsub -q s -o "log/$(basename ${file}).log" "./Run_MyMC.sh ${file}"
done
