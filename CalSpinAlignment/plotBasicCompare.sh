#!/bin/bash

for i in 0 1 2 3 4
#for i in 0
do 
  root -l -b -q plotKaonTTrees_Compare_BasicVariables.C\(4,0,${i},1,-0.5,0.5,\"n0p5y0p5\"\) 
  root -l -b -q plotKaonTTrees_Compare_BasicVariables.C\(4,0,${i},1,-1.5,-0.5,\"n1p5yn0p5\"\) 
  root -l -b -q plotKaonTTrees_Compare_BasicVariables.C\(4,0,${i},1,0.5,1.5,\"0p5y1p5\"\) 
  root -l -b -q plotKaonTTrees_Compare_BasicVariables.C\(4,0,${i},1,-1.0,-0.8,\"n1yn0p8\"\) 
  root -l -b -q plotKaonTTrees_Compare_BasicVariables.C\(4,0,${i},1,0.8,1.0,\"0p8y1\"\) 
done
