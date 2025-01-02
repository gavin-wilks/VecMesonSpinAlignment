#!/bin/tcsh

rm hadd.list
ls -d -1 /gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/PhiSE_RapidityDependence_EtaMode0/* > hadd.list 
./submitEventPlaneHadd.sh SE 0

rm hadd.list
ls -d -1 /gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/PhiSE_RapidityDependence_EtaMode3/* > hadd.list 
./submitEventPlaneHadd.sh SE 3

rm hadd.list
ls -d -1 /gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/PhiSE_RapidityDependence_EtaMode4/* > hadd.list 
./submitEventPlaneHadd.sh SE 4

rm hadd.list
ls -d -1 /gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/PhiSE_RapidityDependence_EtaMode5/* > hadd.list 
./submitEventPlaneHadd.sh SE 5

rm hadd.list
ls -d -1 /gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/PhiME_RapidityDependence_EtaMode0/* > hadd.list 
./submitEventPlaneHadd.sh ME 0

rm hadd.list
ls -d -1 /gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/PhiME_RapidityDependence_EtaMode3/* > hadd.list 
./submitEventPlaneHadd.sh ME 3

rm hadd.list
ls -d -1 /gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/PhiME_RapidityDependence_EtaMode4/* > hadd.list 
./submitEventPlaneHadd.sh ME 4

rm hadd.list
ls -d -1 /gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu19GeV_2019/OutPut/FlowYields/PhiME_RapidityDependence_EtaMode5/* > hadd.list 
./submitEventPlaneHadd.sh ME 5
