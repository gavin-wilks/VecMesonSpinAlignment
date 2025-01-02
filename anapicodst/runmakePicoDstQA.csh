#!/bin/tcsh
starver SL21c
# void runmakePicoDstQA(TString InputFileList,Int_t nevents,TString OutputFile, TString jobId, Int_t mEnergy, Int_t mGid, Int_t mPid)
# csh runmakePicoDstQA.csh ../FileLists/14GeV_2019/Kplus_embed.list 0 test.histo 1 3 11 0

root4star -l -b -q runmakePicoDstQA.C\(\"$1\",$2,\"$3\",\"$4\",$5,$6,$7\)
#root4star -l -b -q runmakePicoDstQA_Phi.C\(\"$1\",$2,\"$3\",\"$4\",$5,$6,$7\)

