#!/bin/bash
order=$1
./submitAllCent_EtaMode0.sh ${order} 
./submitAllCent_EtaMode3.sh ${order}
./submitAllCent_EtaMode4.sh ${order}
./submitAllCent_EtaMode5.sh ${order}
