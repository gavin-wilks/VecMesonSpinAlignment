#!/bin/tcsh
starver SL21c
# void makePicoDstQA(TString InputFileList, Int_t nEvents = 0, TString OutputFile = "test.histo", TString jobId = "1", Int_t mEnergy = 4, Int_t mGid = 11);
# csh runmakePicoDstQA.csh ../ 0 test.histo 1 4 11 0

root4star -l -b -q runmakePicoDstQA.C\(\"$1\",$2,\"$3\",\"$4\",$5,$6,$7\)

