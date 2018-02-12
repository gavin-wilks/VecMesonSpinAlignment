#!/bin/bash
date

if [ $# -eq 0 ]
then
  Energy=200GeV
  SM=ME

  OutPutList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/failed_"$Energy"_"$SM".list"
  rm $OutPutList
  touch $OutPutList

  LogDirectory="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/Log/SLURM"
  InPutList="./completed_"$Energy".log"
  grep -l "Work done" $LogDirectory/*${SM}*.log > $InPutList

  TempList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/Temp_"$Energy".list"
  for item in `cat $InPutList`
  do
    # cat $item | grep "list" >> $TempList
    cat $item | grep "Processing ZdcSmdTree.C" >> $TempList
  done
  # sed -i "s/root4star -b -q -x '//g" $TempList
  sed -i 's/Processing ZdcSmdTree.C("//g' $TempList
  sed -i 's/\.list.*/.list/' $TempList

  CompletedList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/completed_"$Energy".list"
  cat $TempList | sort > $CompletedList
  rm $TempList

  OriginList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/origin_"$Energy".list"
  sort ./submit_${Energy}.list > $OriginList

  comm -13 $CompletedList $OriginList > $OutPutList
  rm $OriginList
  rm $CompletedList

  OutPutROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/deleteROOT_"$Energy"_"$SM".list"

  CompletedROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/completedROOT_"$Energy".list"
  TempROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/TempCompletedROOT_"$Energy".list"
  cat $InPutList > $TempROOT
  sed -i 's/Log/SpinAlignment/g' $TempROOT
  # sed -i 's/SLURM/ZDCSMD\/GainCorrPar/g' $TempROOT # mode = 0
  # sed -i 's/SLURM/ZDCSMD\/ReCenterPar/g' $TempROOT # mode = 1
  # sed -i 's/SLURM/ZDCSMD\/ShiftPar/g' $TempROOT # mode = 2
  # sed -i 's/SLURM/ZDCSMD\/ShiftParFull/g' $TempROOT # mode = 3
  # sed -i 's/SLURM/ZDCSMD\/Resolution/g' $TempROOT # mode = 4
  # sed -i 's/SLURM/ZDCSMD\/DirectedFlow/g' $TempROOT # mode = 5
  sed -i 's/SLURM/ZDCSMD\/Phi\/Forest/g' $TempROOT # mode = 6
  sed -i 's/log/root/g' $TempROOT
  sed -i 's/phi.E/file_/g' $TempROOT
  # sed -i 's/GeV_\(.*\)_/GeV_GainCorrPar_/g' $TempROOT # mode = 0
  # sed -i 's/GeV_\(.*\)_/GeV_ReCenterPar_/g' $TempROOT # mode = 1
  # sed -i 's/GeV_\(.*\)_/GeV_ShiftPar_/g' $TempROOT # mode = 2
  # sed -i 's/GeV_\(.*\)_/GeV_ShiftParFull_/g' $TempROOT # mode = 3
  # sed -i 's/GeV_\(.*\)_/GeV_Resolution_/g' $TempROOT # mode = 4
  # sed -i 's/GeV_\(.*\)_/GeV_DirectedFlow_/g' $TempROOT # mode = 5
  # sed -i 's/GeV_\(.*\)_/GeV_Phi_SE_/g' $TempROOT # mode = 6 => SE
  sed -i 's/GeV_\(.*\)_/GeV_Phi_ME_/g' $TempROOT # mode = 6 => ME
  cat $TempROOT | sort > $CompletedROOT
  rm $TempROOT

  OriginROOT="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/submit/originROOT_"$Energy".list"
  OriginForest="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/SpinAlignment/ZDCSMD"
  # ls -d $OriginForest/GainCorrPar/*GainCorrPar_*.root | sort > $OriginROOT #mode = 0 => gain correction
  # ls -d $OriginForest/ReCenterPar/*ReCenterPar_*.root | sort > $OriginROOT #mode = 1 => re-center correction
  # ls -d $OriginForest/ShiftPar/*ShiftPar_*.root | sort > $OriginROOT #mode = 2 => re-center correction
  # ls -d $OriginForest/ShiftParFull/*ShiftParFull_*.root | sort > $OriginROOT #mode = 3 => re-center correction
  # ls -d $OriginForest/Resolution/*Resolution_*.root | sort > $OriginROOT #mode = 4 => re-center correction
  # ls -d $OriginForest/DirectedFlow/*DirectedFlow_*.root | sort > $OriginROOT #mode = 5 => re-center correction
  ls -d $OriginForest/Phi/Forest/*${SM}*.root | sort > $OriginROOT # mode = 6 => phi production mode

  comm -13 $CompletedROOT $OriginROOT > $OutPutROOT
  # # sed -i 's/Script/SpinAlignment\/Resolution/g' $OutPutROOT # resolution mode
  # # sed -i "s/runPhiSE${Energy}/file_${Energy}_Resolution_/g" $OutPutROOT
  # sed -i 's/Script/SpinAlignment\/Phi\/Forest/g' $OutPutROOT # phi-meson production mode
  # sed -i "s/runPhi${SM}${Energy}/file_${Energy}_Phi_${SM}_/g" $OutPutROOT
  # sed -i "s/csh/root/g" $OutPutROOT

  rm $InPutList
  rm $OriginROOT
  rm $CompletedROOT
else
  echo "Wrong number of parameters"
fi
