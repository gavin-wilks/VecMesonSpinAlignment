#!/bin/bash
date

#. ./VecMesonTree.sh

if [ $# -eq 0 ]
then
  Energy=27GeV
  # Energy=200GeV
  SM=SE

  OutPutList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/List/AuAu${Energy}/failed_"$Energy"_"$SM".list"
  rm $OutPutList
  touch $OutPutList

  LogDirectory="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/SLURM"
  InPutList="./List/AuAu${Energy}/completed_"$Energy".log"
  grep -l "Work done" $LogDirectory/*${SM}*.log > $InPutList

  TempList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/List/AuAu${Energy}/Temp_"$Energy".list"
  for item in `cat $InPutList`
  do
    # cat $item | grep "list" >> $TempList
    cat $item | grep "Processing VecMesonTree.C" >> $TempList
  done
  # sed -i "s/root4star -b -q -x '//g" $TempList
  sed -i 's/Processing VecMesonTree.C("//g' $TempList
  sed -i 's/\.list.*/.list/' $TempList

  CompletedList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/List/AuAu${Energy}/completed_"$Energy".list"
  cat $TempList | sort > $CompletedList
  rm $TempList

  OriginList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/List/AuAu${Energy}/origin_"$Energy".list"
  sort ./List/AuAu${Energy}/submit_${Energy}.list > $OriginList

  comm -13 $CompletedList $OriginList > $OutPutList
  rm $OriginList
  rm $CompletedList

  OutPutROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/List/AuAu${Energy}/deleteROOT_"$Energy"_"$SM".list"
  rm $OutPutROOT
  touch $OutPutROOT

  CompletedROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/List/AuAu${Energy}/completedROOT_"$Energy".list"
  TempROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/List/AuAu${Energy}/TempCompletedROOT_"$Energy".list"
  cat $InPutList > $TempROOT
  sed -i 's/Log/SpinAlignment/g' $TempROOT
  sed -i 's/log/root/g' $TempROOT
  sed -i 's/phi.E/file_/g' $TempROOT
  # sed -i 's/SLURM/ReCenterParameter/g' $TempROOT # mode = 0
  # sed -i 's/GeV_\(.*\)_/GeV_ReCenterPar_/g' $TempROOT # mode = 0
  # sed -i 's/SLURM/ShiftParameter/g' $TempROOT # mode = 1
  # sed -i 's/GeV_\(.*\)_/GeV_ShiftPar_/g' $TempROOT # mode = 1
  sed -i 's/SLURM/Resolution/g' $TempROOT # mode = 2
  sed -i 's/GeV_\(.*\)_/GeV_Resolution_/g' $TempROOT # mode = 2
  # sed -i 's/SLURM/Phi\/Forest/g' $TempROOT # mode = 3
  # sed -i 's/GeV_\(.*\)_/GeV_Phi_SE_/g' $TempROOT # mode = 3 => SE
  cat $TempROOT | sort > $CompletedROOT
  rm $TempROOT

  OriginROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/List/AuAu${Energy}/originROOT_"$Energy".list"
  OriginForest="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/SpinAlignment"
  # ls -d $OriginForest/ReCenterParameter/*ReCenterPar_*.root | sort > $OriginROOT #mode = 0 => re-center correction
  # ls -d $OriginForest/ShiftParameter/*ShiftPar_*.root | sort > $OriginROOT #mode = 1 => shift correction
  ls -d $OriginForest/Resolution/*Resolution_*.root | sort > $OriginROOT #mode = 2 => shift correction
  # ls -d $OriginForest/Phi/Forest/*${SM}*.root | sort > $OriginROOT

  comm -13 $CompletedROOT $OriginROOT > $OutPutROOT

  rm $InPutList
  rm $OriginROOT
  rm $CompletedROOT
else
  echo "Wrong number of parameters"
fi
