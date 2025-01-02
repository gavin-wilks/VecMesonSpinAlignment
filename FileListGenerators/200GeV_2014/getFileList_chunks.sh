#!/bin/bash

limit=10000

for i in {0..10..1}
do
  rm pico_chunk_${i}.list
  get_file_list.pl -keys 'path,filename' -cond 'production=P21ic,library=SL21c,trgsetupname=production_19GeV_2019,filetype=daq_reco_picoDst,filename~st_physics,storage=nfs' -start `expr ${i} \* ${limit}` -limit ${limit} -delim '/' -distinct > pico_chunk_${i}.list
done
