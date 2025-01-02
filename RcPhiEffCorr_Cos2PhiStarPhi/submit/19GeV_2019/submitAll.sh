#!/bin/bash

#1st Order
#./submitTreeProduction.sh 0 2 3 0 0 1 
#./submitTreeProduction.sh 0 2 3 0 3 1
#./submitTreeProduction.sh 0 2 3 0 4 1
#./submitTreeProduction.sh 0 2 3 0 5 1
#./submitTreeProduction.sh 1 4 4 0 0 1
#./submitTreeProduction.sh 1 4 4 0 3 1
#./submitTreeProduction.sh 1 4 4 0 4 1
#./submitTreeProduction.sh 1 4 4 0 5 1
#./submitTreeProduction.sh 2 5 5 0 0 1
#./submitTreeProduction.sh 2 5 5 0 3 1
#./submitTreeProduction.sh 2 5 5 0 4 1
#./submitTreeProduction.sh 2 5 5 0 5 1
#
#./submitTreeProduction.sh 0 0 1 1 0 1
#./submitTreeProduction.sh 0 0 1 1 3 1
#./submitTreeProduction.sh 0 0 1 1 4 1
#./submitTreeProduction.sh 0 0 1 1 5 1
#./submitTreeProduction.sh 3 2 2 1 0 1
#./submitTreeProduction.sh 3 2 2 1 3 1
#./submitTreeProduction.sh 3 2 2 1 4 1
#./submitTreeProduction.sh 3 2 2 1 5 1
#
#./submitTreeProduction.sh 0 0 0 2 0 1
#./submitTreeProduction.sh 0 0 0 2 3 1
#./submitTreeProduction.sh 0 0 0 2 4 1
#./submitTreeProduction.sh 0 0 0 2 5 1
#./submitTreeProduction.sh 3 1 1 2 0 1
#./submitTreeProduction.sh 3 1 1 2 3 1
#./submitTreeProduction.sh 3 1 1 2 4 1
#./submitTreeProduction.sh 3 1 1 2 5 1
#
#./submitTreeProduction.sh 4 2 3 -0.3 1000 0.3333 0.0 0.0 0.0 

#./submitTreeProduction.sh 0 2 3 -0.3 1000 0.3333 0.0 0.0 0.0 
#./submitTreeProduction.sh 1 4 4 -0.3 1000 0.3333 0.0 0.0 0.0
#./submitTreeProduction.sh 2 5 5 -0.3 1000 0.3333 0.0 0.0 0.0
#
#./submitTreeProduction.sh 0 2 3 -0.2 1000 0.3333 0.0 0.0 0.0 
#./submitTreeProduction.sh 1 4 4 -0.2 1000 0.3333 0.0 0.0 0.0
#./submitTreeProduction.sh 2 5 5 -0.2 1000 0.3333 0.0 0.0 0.0
#
#./submitTreeProduction.sh 0 2 3 -0.1 1000 0.3333 0.0 0.0 0.0 
#./submitTreeProduction.sh 1 4 4 -0.1 1000 0.3333 0.0 0.0 0.0
#./submitTreeProduction.sh 2 5 5 -0.1 1000 0.3333 0.0 0.0 0.0
#
#./submitTreeProduction.sh 0 2 3  0.1 1000 0.3333 0.0 0.0 0.0 
#./submitTreeProduction.sh 1 4 4  0.1 1000 0.3333 0.0 0.0 0.0
#./submitTreeProduction.sh 2 5 5  0.1 1000 0.3333 0.0 0.0 0.0
#
#./submitTreeProduction.sh 0 2 3  0.2 1000 0.3333 0.0 0.0 0.0 
#./submitTreeProduction.sh 1 4 4  0.2 1000 0.3333 0.0 0.0 0.0
#./submitTreeProduction.sh 2 5 5  0.2 1000 0.3333 0.0 0.0 0.0
#
#./submitTreeProduction.sh 0 2 3  0.3 1000 0.3333 0.0 0.0 0.0 
#./submitTreeProduction.sh 1 4 4  0.3 1000 0.3333 0.0 0.0 0.0
#./submitTreeProduction.sh 2 5 5  0.3 1000 0.3333 0.0 0.0 0.0
#
###2nd Order
for i in 0.20 0.205 0.21 0.215 0.22 0.225 0.23 0.235 0.24 0.245 0.25 0.255 0.26 0.265 0.27 0.275 0.28 0.285 0.29 0.295 0.30 
do 
  ./submitTreeProduction.sh 4 2 5 ${i}
done
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3000 0.0 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3000 0.0 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3000 0.0 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.2667 0.0 0.0 0.0 
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.2667 0.0 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.2667 0.0 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.2333 0.0 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.2333 0.0 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.2333 0.0 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.2000 0.0 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.2000 0.0 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.2000 0.0 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3667 0.0 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3667 0.0 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3667 0.0 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.4000 0.0 0.0 0.0   
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.4000 0.0 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.4000 0.0 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.4333 0.0 0.0 0.0 
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.4333 0.0 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.4333 0.0 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.4667 0.0 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.4667 0.0 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.4667 0.0 0.0 0.0 
#
###########################################################
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333 -0.3 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333 -0.3 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333 -0.3 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333 -0.2 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333 -0.2 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333 -0.2 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333 -0.1 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333 -0.1 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333 -0.1 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.1 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.1 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.1 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.2 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.2 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.2 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.3 0.0 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.3 0.0 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.3 0.0 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0 -0.3 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0 -0.3 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0 -0.3 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0 -0.2 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0 -0.2 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0 -0.2 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0 -0.1 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0 -0.1 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0 -0.1 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0  0.1 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0  0.1 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0  0.1 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0  0.2 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0  0.2 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0  0.2 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0  0.3 0.0  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0  0.3 0.0 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0  0.3 0.0 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0  0.0 -0.3  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0  0.0 -0.3 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0  0.0 -0.3 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0  0.0 -0.2  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0  0.0 -0.2 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0  0.0 -0.2 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0  0.0 -0.1  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0  0.0 -0.1 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0  0.0 -0.1 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0  0.0  0.1  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0  0.0  0.1 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0  0.0  0.1 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0  0.0  0.2  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0  0.0  0.2 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0  0.0  0.2 
#
#./submitTreeProduction.sh 0 2 3 0.0 1000 0.3333  0.0  0.0  0.3  
#./submitTreeProduction.sh 1 4 4 0.0 1000 0.3333  0.0  0.0  0.3 
#./submitTreeProduction.sh 2 5 5 0.0 1000 0.3333  0.0  0.0  0.3 

#################################################################3

#./submitTreeProduction.sh 0 2 3 0 0 2 0.00 0.0625 0.3333 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.00 0.0625 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.00 0.0625 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.00 0.125 0.3333 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.00 0.125 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.00 0.125 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.00 0.25 0.3333 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.00 0.25 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.00 0.25 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.00 0.5 0.3333 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.00 0.5 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.00 0.5 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.00 1 0.3333 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.00 1 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.00 1 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.00 10 0.3333 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.00 10 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.00 10 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.00 100 0.3333 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.00 100 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.00 100 0.3333



#./submitTreeProduction.sh 0 2 3 0 0 2 0.15 1000 0.3333  
#./submitTreeProduction.sh 1 4 4 0 0 2 0.15 1000 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.15 1000 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.30 1000 0.3333
#./submitTreeProduction.sh 1 4 4 0 0 2 0.30 1000 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.30 1000 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.45 1000 0.3333
#./submitTreeProduction.sh 1 4 4 0 0 2 0.45 1000 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.45 1000 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.60 1000 0.3333
#./submitTreeProduction.sh 1 4 4 0 0 2 0.60 1000 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 0.60 1000 0.3333
#
#
#./submitTreeProduction.sh 0 2 3 0 0 2 -0.15 1000 0.3333  
#./submitTreeProduction.sh 1 4 4 0 0 2 -0.15 1000 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 -0.15 1000 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 -0.30 1000 0.3333
#./submitTreeProduction.sh 1 4 4 0 0 2 -0.30 1000 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 -0.30 1000 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 -0.45 1000 0.3333
#./submitTreeProduction.sh 1 4 4 0 0 2 -0.45 1000 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 -0.45 1000 0.3333
#
#./submitTreeProduction.sh 0 2 3 0 0 2 -0.60 1000 0.3333
#./submitTreeProduction.sh 1 4 4 0 0 2 -0.60 1000 0.3333
#./submitTreeProduction.sh 2 5 5 0 0 2 -0.60 1000 0.3333

#./submitTreeProduction.sh 0 2 3 0 0 2 0.0 2 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.0 2 
#./submitTreeProduction.sh 2 5 5 0 0 2 0.0 2 
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.0 1 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.0 1 
#./submitTreeProduction.sh 2 5 5 0 0 2 0.0 1 
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.0 0.5
#./submitTreeProduction.sh 1 4 4 0 0 2 0.0 0.5
#./submitTreeProduction.sh 2 5 5 0 0 2 0.0 0.5 
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.0 0.25
#./submitTreeProduction.sh 1 4 4 0 0 2 0.0 0.25
#./submitTreeProduction.sh 2 5 5 0 0 2 0.0 0.25
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.0 0.125
#./submitTreeProduction.sh 1 4 4 0 0 2 0.0 0.125
#./submitTreeProduction.sh 2 5 5 0 0 2 0.0 0.125
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.0 0.0625
#./submitTreeProduction.sh 1 4 4 0 0 2 0.0 0.0625
#./submitTreeProduction.sh 2 5 5 0 0 2 0.0 0.0625

#./submitTreeProduction.sh 0 2 3 0 0 2 -0.1 
#./submitTreeProduction.sh 1 4 4 0 0 2 -0.1
#./submitTreeProduction.sh 2 5 5 0 0 2 -0.1
#
#./submitTreeProduction.sh 0 2 3 0 0 2 -0.2 
#./submitTreeProduction.sh 1 4 4 0 0 2 -0.2
#./submitTreeProduction.sh 2 5 5 0 0 2 -0.2
#
#./submitTreeProduction.sh 0 2 3 0 0 2 -0.3 
#./submitTreeProduction.sh 1 4 4 0 0 2 -0.3
#./submitTreeProduction.sh 2 5 5 0 0 2 -0.3
#
#./submitTreeProduction.sh 0 2 3 0 0 2 -0.4 
#./submitTreeProduction.sh 1 4 4 0 0 2 -0.4
#./submitTreeProduction.sh 2 5 5 0 0 2 -0.4
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.1 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.1
#./submitTreeProduction.sh 2 5 5 0 0 2 0.1
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.2 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.2
#./submitTreeProduction.sh 2 5 5 0 0 2 0.2
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.3 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.3
#./submitTreeProduction.sh 2 5 5 0 0 2 0.3
#
#./submitTreeProduction.sh 0 2 3 0 0 2 0.4 
#./submitTreeProduction.sh 1 4 4 0 0 2 0.4
#./submitTreeProduction.sh 2 5 5 0 0 2 0.4


#./submitTreeProduction.sh 2 5 5 0 3 2
#./submitTreeProduction.sh 2 5 5 0 4 2
#./submitTreeProduction.sh 2 5 5 0 5 2
#
#./submitTreeProduction.sh 0 0 1 1 0 2
#./submitTreeProduction.sh 0 0 1 1 3 2
#./submitTreeProduction.sh 0 0 1 1 4 2
#./submitTreeProduction.sh 0 0 1 1 5 2
#./submitTreeProduction.sh 3 2 2 1 0 2
#./submitTreeProduction.sh 3 2 2 1 3 2
#./submitTreeProduction.sh 3 2 2 1 4 2
#./submitTreeProduction.sh 3 2 2 1 5 2

#./submitTreeProduction.sh 0 0 0 2 0 2
#./submitTreeProduction.sh 0 0 0 2 3 2
#./submitTreeProduction.sh 0 0 0 2 4 2
#./submitTreeProduction.sh 0 0 0 2 5 2
#./submitTreeProduction.sh 3 1 1 2 0 2
#./submitTreeProduction.sh 3 1 1 2 3 2
#./submitTreeProduction.sh 3 1 1 2 4 2
#./submitTreeProduction.sh 3 1 1 2 5 2
