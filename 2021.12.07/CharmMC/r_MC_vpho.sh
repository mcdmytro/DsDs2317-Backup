rm -f log/*
rm -f histos/*
rm -f idx_files/*
rm -f root_files/*

for file in /group/belle/bdata_b/mcprod/dat/e000061/evtgen/charm/02/all/0127/continuum/07/evtgen-charm-02-all-e000061r000700-b20090127_0910.mdst; do
	bsub -q s -o "log/$(basename ${file}).log" "./Run_MyMC.sh ${file}"
done
