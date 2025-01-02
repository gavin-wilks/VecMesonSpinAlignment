#!/bin/bash
for i in {0..8}
do 
  root4star -l -b -q toyMcPhiDecay.C\(0,0,${i},0,1,1,0,\"20221209\"\)
done 
