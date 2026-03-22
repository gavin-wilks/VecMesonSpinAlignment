#!/bin/bash

#opt=$1
#
#if [ $opt == 0 ]
#then 
##./findnosubmitTree.sh ME SL20d P18ih high
#
#elif [ $opt == 1 ] 
#then
#fi


./findnosubmitTree.sh B128660A47711BBB9AEFF30DF85AD067 SE SL22c P16id low
./preparefilelist.sh  B128660A47711BBB9AEFF30DF85AD067 SL22c P16id low SE

./findnosubmitTree.sh 562143AC09870E4CB2FB8C1E60A64FD6 SE SL22c P16id mid
./preparefilelist.sh  562143AC09870E4CB2FB8C1E60A64FD6 SL22c P16id mid SE

./findnosubmitTree.sh 5E88D972EE0956AA23008B8CB92782E9 ME SL22c P16id low
./preparefilelist.sh  5E88D972EE0956AA23008B8CB92782E9 SL22c P16id low ME

./findnosubmitTree.sh CB8EFFCAAC1F2E1013FB9846F1DAEF16 ME SL22c P16id mid
./preparefilelist.sh  CB8EFFCAAC1F2E1013FB9846F1DAEF16 SL22c P16id mid ME
