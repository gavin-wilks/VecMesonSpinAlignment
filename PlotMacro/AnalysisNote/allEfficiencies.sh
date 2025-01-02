#!/bin/bash

root -l -b -q plotMcPhiEff.C\(3,0,0,1,\"\"\)
root -l -b -q plotMcPhiEff.C\(3,0,0,2,\"\"\)
root -l -b -q plotMcPhiEff.C\(4,0,0,1,\"\"\)
root -l -b -q plotMcPhiEff.C\(4,0,0,2,\"_FixedRes\"\)

#root -l -b -q plotMcPhiEff_Centrality.C\(3,0,0,1\)
#root -l -b -q plotMcPhiEff_Centrality.C\(3,0,0,2\)
#root -l -b -q plotMcPhiEff_Centrality.C\(4,0,0,1\)
#root -l -b -q plotMcPhiEff_Centrality.C\(4,0,0,2\)
#
#root -l -b -q plotMcPhiEff_Rapidity.C\(3,0,0,1\)
#root -l -b -q plotMcPhiEff_Rapidity.C\(3,0,0,2\)
#root -l -b -q plotMcPhiEff_Rapidity.C\(4,0,0,1\)
#root -l -b -q plotMcPhiEff_Rapidity.C\(4,0,0,2\)
