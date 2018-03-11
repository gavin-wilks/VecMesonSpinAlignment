#!/bin/bash
date

#. ./cleanOutPut.sh

if [ $# -eq 0 ]
  then
    Energy=27GeV
    # Energy=200GeV
    SM=SE
    InPutList="./List/AuAu${Energy}/deleteROOT_${Energy}_${SM}.list"
    echo "delete following files from "
    echo $InPutList
    for item in `cat $InPutList`
    do
      echo $item
      rm $item
    done
    rm $InPutList

  else
    echo "Wrong number of parameters"
fi
