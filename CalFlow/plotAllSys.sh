#!/bin/bash

energy=$1

for i in {0..12}
  do
    root -l -b -q calSysErrorKStar.C\($energy,2,$i\)
  done
