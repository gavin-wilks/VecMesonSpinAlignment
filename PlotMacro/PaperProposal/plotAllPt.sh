#!/bin/bash

for energy in 0 1 
do 
  for order in 1 2 
  do
    root -l -b -q plotSys_Dca.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSys_NSig.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSys_Normalization.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSys_Poly.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSysCent_Dca.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSysCent_NSig.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSysCent_Normalization.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSysCent_Poly.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSysY_Dca.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSysY_NSig.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSysY_Normalization.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSysY_Poly.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
    root -l -b -q plotSysY_Diff.C\(${energy},${order},\"AccRes\",\"eta1_eta1\"\)
  done
done
