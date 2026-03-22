#!/bin/bash
date

#. ./hadd_phi.sh


if [ $# -eq 0 ]
  then
    # Mode=GainCorrPar # mode = 0
    # Mode=ReCenterPar # mode = 1
    # Mode=ShiftPar # mode = 2
    # Mode=ShiftParFull # mode = 3
    # Mode=Resolution # mode = 4
    Mode=DirectedFlow # mode = 5
    Energy=200GeV
    OutDir="/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu$Energy/SpinAlignment/ZDCSMD/$Mode/merged_file/merged_${Mode}_"
    suffix=".root"
    InPutList="/global/homes/x/xusun/STAR/VecMesonSpinAlignment/TreeProduction/hadd/List/AuAu$Energy/$Mode/merged_${Mode}.list"
    counter=0
    for item in `cat $InPutList`
    do
      cp ./run_hadd.csh ./run_hadd_${Mode}_$counter.csh
      echo "cd ./AuAu$Energy/SpinAlignment/ZDCSMD/$Mode" >> run_hadd_${Mode}_$counter.csh
      echo " " >> run_hadd_${Mode}_$counter.csh

      OutName=$OutDir$counter$suffix

      echo "rm $OutName " >> run_hadd_${Mode}_$counter.csh
      echo " " >> run_hadd_${Mode}_$counter.csh
      echo -n "/usr/bin/time -v hadd $OutName " >> run_hadd_${Mode}_$counter.csh

      for yields in `cat $item`
      do
	echo -n "$yields " >> run_hadd_${Mode}_$counter.csh
      done

      echo " " >> run_hadd_${Mode}_$counter.csh
      echo " " >> run_hadd_${Mode}_$counter.csh
      echo "echo 'This is the end of hadd\!\!\!'" >> run_hadd_${Mode}_$counter.csh

      # sbatch -p shared-chos -t 1:00:00 -o /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_${Mode}_$counter.log -e /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Log/hadd/job_${Mode}_$counter.err ./run_hadd_${Mode}_$counter.csh

      mv run_hadd_${Mode}_$counter.csh /global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu${Energy}/Script/hadd/
      let "counter=counter+1"
    done

  else
    echo "Wrong number of parameters"
fi
