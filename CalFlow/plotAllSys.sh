#!/bin/bash

energy=$1

for i in {0..12}
  do
    root -l -b -q calSysError.C\($energy,0,$i\)
  done
