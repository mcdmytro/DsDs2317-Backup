#! /bin/tcsh -f

echo '#! /bin/csh -f'
echo ''
echo 'setenv USE_GRAND_REPROCESS_DATA 1'
echo 'setenv BELLE_LEVEL b20090127_0910'
echo 'setenv BELLE_DEBUG opt'
echo 'setenv CERN_LEVEL  2006'

echo 'source /sw/belle/local/etc/cshrc_general'

set id = $1
set job = $2
set exp = $3
set evts = $4
set type = $5
set basfexp = $6

set dat_file = gen/${id}/'evtgen_exp_'${basfexp}_${id}-${job}.gen
set log_file = log/${id}/'evtgen_exp_'${basfexp}_${id}-${job}.log

echo 'basf <<EOF >&' ${log_file}

echo 'module register evtgen'

echo 'path create main'
echo 'path add_module main evtgen'
echo 'path add_module main genpak'

echo 'module put_parameter evtgen USER_DECAY\\'$id.dec
echo 'module put_parameter evtgen PRINT_FLAG\\2'

echo ''
cat ${type}
echo ''

echo 'module put_parameter evtgen MODEL_TYPE\\0'

echo 'module put_parameter evtgen LER_ENERGY\\3.499218'
echo 'module put_parameter evtgen HER_ENERGY\\7.998213'

echo 'module put_parameter evtgen ExpNo\'${exp}
echo 'module put_parameter evtgen RunNo\\0'

echo 'initialize'

echo 'table savebr gen_info'
echo 'table savebr gen_beam'
echo 'table save gen_pakevt'
echo 'table save gen_pakvtx'
echo 'table save belle_event'
echo 'table save mctype_all'
echo 'table save mctype'
echo 'table save hepevt_all'

echo 'output open '${dat_file}

echo 'generate_event '${evts}

echo 'terminate'
echo 'EOF'

exit
