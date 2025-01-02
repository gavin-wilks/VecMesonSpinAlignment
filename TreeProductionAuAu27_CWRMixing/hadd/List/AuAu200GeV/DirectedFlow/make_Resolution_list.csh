#!/bin/csh
touch DirectedFlow.list
cd /project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/SpinAlignment/ZDCSMD/DirectedFlow/
ls -d file_*.root >! /global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/hadd/List/AuAu200GeV/DirectedFlow/DirectedFlow.list
cd -
