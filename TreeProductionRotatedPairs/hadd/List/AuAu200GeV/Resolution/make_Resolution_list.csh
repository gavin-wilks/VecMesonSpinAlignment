#!/bin/csh
touch Resolution.list
cd /project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/SpinAlignment/ZDCSMD/Resolution/
ls -d file_*.root >! /global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/hadd/List/AuAu200GeV/Resolution/Resolution.list
cd -
