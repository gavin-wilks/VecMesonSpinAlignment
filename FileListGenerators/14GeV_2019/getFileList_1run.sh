#!/bin/bash
rm pico_1run.list
rm pico_1run_sorted.list
#rm runNumber_1run.list

get_file_list.pl -keys 'path,filename' -cond 'production=P21ic,library=SL21c,trgsetupname=production_19GeV_2019,filetype=daq_reco_picoDst,runnumber~20061036||20080024||20057004||20093036||20076031,filename~st_physics,storage=nfs' -limit 0 -delim '/' -distinct > pico_1run.list

awk -F/ '{print $NF, $0}' pico_1run.list | sort | uniq | cut -f2- -d ' ' > pico_1run_sorted.list
#awk -F/ '{print $(NF-1), $0}' pico_prod.list | cut -d ' ' -f 1 | sort | uniq > runNumber_prod.list
