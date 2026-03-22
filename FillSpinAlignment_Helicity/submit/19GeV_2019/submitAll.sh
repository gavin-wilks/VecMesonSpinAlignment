#!/bin/bash

for i in 3
do 
  ./submitTreeProduction.sh 0 0 SE ${i} 
  ./submitTreeProduction.sh 0 1 ME ${i} 
done
#./submitEpdTreeProduction.sh 3 0 SE
#./submitEpdTreeProduction.sh 3 1 ME
#./submitEpdTreeProduction.sh 4 0 SE
#./submitEpdTreeProduction.sh 4 1 ME
#./submitEpdTreeProduction.sh 5 0 SE
#./submitEpdTreeProduction.sh 5 1 ME

