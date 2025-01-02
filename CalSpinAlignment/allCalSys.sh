#!/bin/bash

# inputs ==> energy, eta, pid, year, date, random3D
for i in {0..13}
  do
    root -l -b -q subBackGroundPhiEta.C\(4,${i},0,0,\"20240511\",0,\"eta1_eta1\"\)
    root -l -b -q subBackGroundPhiEta.C\(4,${i},0,0,\"20240511\",0,\"eta1_eta1\"\)
done
