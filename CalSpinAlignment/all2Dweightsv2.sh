#!/bin/bash

for v2 in 0.0333 0.0667 0.1000 0.1333 0.1667 0.2000 0.2333 0.2667 0.3000
do
  root -l -b -q calculate_2D_weights.C\(4,${v2}\) 
done
