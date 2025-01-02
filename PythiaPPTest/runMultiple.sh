#!/bin/bash

for i in {21..30}
do
  root -l -b -q scan_phi_mesons.C\(10000000,\"${i}\"\)
done 
