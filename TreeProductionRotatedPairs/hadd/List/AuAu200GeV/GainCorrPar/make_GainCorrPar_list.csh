#!/bin/csh
touch GainCorrPar.list
cd /project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/SpinAlignment/ZDCSMD/GainCorrPar/
ls -d file_*.root >! /global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/hadd/List/AuAu200GeV/GainCorrPar/GainCorrPar.list
cd -
