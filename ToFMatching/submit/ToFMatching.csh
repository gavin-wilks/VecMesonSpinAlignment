#!/bin/tcsh 
# Note, all code in this scrpt is executed under tcsh
if (  "$SLURM_JOB_PARTITION" =~ *"chos"  ) then  
    echo  task-in-chos
    chosenv
    ls -l /proc/chos/link
else
    echo  task-in-shifter
    echo inShifter:`env|grep  SHIFTER_RUNTIME`
    cat /etc/*release
    #
    # - - - -  D O   N O T  T O U C H  T H I S   S E C T I O N- - - - 
    #
    whoami    
    echo  load STAR enviroment in shifter
    set NCHOS = sl64
    set SCHOS = 64
    set DECHO = 1
    set SCRATCH = $WRK_DIR/out-star1
    setenv GROUP_DIR /common/star/star${SCHOS}/group/
    source $GROUP_DIR/star_cshrc.csh    
    #
    # - - - -   Y O U   C A N   C H A N G E   B E L O W  - - - -
    #    
endif  

    cd ${WRK_DIR} # important
    set daqN = $argv[1]
     
    echo  starting Kaon ToF matching efficiency analysis with PATH_XRD=$PATH_XRD, daqN=$daqN, execName=$EXEC_NAME, workerName=`hostname -f`, startDate=`date`

    echo testing STAR setup $STAR_VER in `pwd`
    starver $STAR_VER 
    env |grep STAR

    echo 'my new STAR ver='$STAR'  test root '
    root4star -b -q 
    if ( $? != 0) then
	echo STAR environment is corrupted, aborting job
	echo $STAR
	which root4star
	exit
    endif
 
    #echo EEEEE ;   exit

    echo `date`" Fire: $EXEC_NAME  [wiat]"
    cd ${WRK_DIR}/VecMesonSpinAlignment/ToFMatching
    echo `pwd`
    /usr/bin/time -v  $EXEC_NAME -b -q ToFMatching.C\(\"$daqN\",${SLURM_ARRAY_TASK_ID},$energy\) >& $LOG_PATH/KaonToFMatching_${Energy}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

