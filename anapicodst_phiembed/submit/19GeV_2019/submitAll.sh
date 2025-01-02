#!/bin/bash

#./submitTreeProduction.sh  0
./submitTreeProduction.sh  1

#for i in 0.2333 0.2667 0.3000 0.3333 0.3667 0.4000 0.4333
#do 
#  ./submitTreeProduction.sh ${i} 0.0 0.0 0.0 0.0 0.3333 0.0 0.0 0.0 0.0
#  ./submitTreeProduction.sh 0.3333 0.0 0.0 0.0 0.0 ${i} 0.0 0.0 0.0 0.0
#done 
#
#for i in -0.3 -0.2 -0.1 0.1 0.2 0.3
#do 
#  ./submitTreeProduction.sh 0.3333 ${i} 0.0 0.0 0.0 0.3333 0.0 0.0 0.0 0.0
#  ./submitTreeProduction.sh 0.3333 0.0 ${i} 0.0 0.0 0.3333 0.0 0.0 0.0 0.0
#  ./submitTreeProduction.sh 0.3333 0.0 0.0 ${i} 0.0 0.3333 0.0 0.0 0.0 0.0
#  ./submitTreeProduction.sh 0.3333 0.0 0.0 0.0 ${i} 0.3333 0.0 0.0 0.0 0.0
#  ./submitTreeProduction.sh 0.3333 0.0 0.0 0.0 0.0 0.3333 ${i} 0.0 0.0 0.0
#  ./submitTreeProduction.sh 0.3333 0.0 0.0 0.0 0.0 0.3333 0.0 ${i} 0.0 0.0
#  ./submitTreeProduction.sh 0.3333 0.0 0.0 0.0 0.0 0.3333 0.0 0.0 ${i} 0.0
#  ./submitTreeProduction.sh 0.3333 0.0 0.0 0.0 0.0 0.3333 0.0 0.0 0.0 ${i}
#done 
