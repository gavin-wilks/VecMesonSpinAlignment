#!/bin/bash
for i in {2..5}
do 
  root4star -l -b -q toyMcPhiDecay.C\(4,5,${i},0,10,1,0,\"20221209\",10.0,10.0\)
done 
