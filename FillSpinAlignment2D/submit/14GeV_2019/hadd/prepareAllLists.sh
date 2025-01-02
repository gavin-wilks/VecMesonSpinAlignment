#!/bin/bash

path=/gpfs01/star/scratch/gwilks3/VectorMesonSpinAlignment/AuAu14GeV_2019/OutPut/SpinAlignmentYields

rm hadd_SE_etamode*.list
rm hadd_ME_etamode*.list


ls -d -1 ${path}/PhiSE_FirstOrder_EtaMode0/* > hadd_SE_etamode0.list
ls -d -1 ${path}/PhiME_FirstOrder_EtaMode0/* > hadd_ME_etamode0.list
ls -d -1 ${path}/PhiSE_FirstOrder_EtaMode3/* > hadd_SE_etamode3.list
ls -d -1 ${path}/PhiME_FirstOrder_EtaMode3/* > hadd_ME_etamode3.list
ls -d -1 ${path}/PhiSE_FirstOrder_EtaMode4/* > hadd_SE_etamode4.list 
ls -d -1 ${path}/PhiME_FirstOrder_EtaMode4/* > hadd_ME_etamode4.list 
ls -d -1 ${path}/PhiSE_FirstOrder_EtaMode5/* > hadd_SE_etamode5.list 
ls -d -1 ${path}/PhiME_FirstOrder_EtaMode5/* > hadd_ME_etamode5.list 
