#!/bin/csh

#rm -rf test.list 
#get_file_list.pl -keys 'runnumber' -cond 'production=P19ib,trgsetupname=27GeV_production_2018,filetype=daq_reco_picoDst,filename~st_physics,storage!=HPSS' -limit 0 -distinct '/' >& Run18_27.runs
get_file_list.pl -keys 'path,filename' -cond 'production=P19ib,trgsetupname=27GeV_production_2018,filetype=daq_reco_picoDst,filename~st_physics,storage!=HPSS,runnumber[]19142020-19142039' -limit 0 -distinct -delim '/' >& Run18_27.list
#get_file_list.pl -keys 'path,filename' -cond 'storage=nfs,filetype=daq_reco_MuDst,filename~st_hlt,production=P16ij,trgsetupname=AuAu_200_production_2016,tpx=1,sanity=1' -limit 0 -distinct -delim '/' >& micro.list

#For picoDst use 'filetype=daq_reco_picoDst'
#For muDst   use 'filetype=daq_reco_MuDst'
