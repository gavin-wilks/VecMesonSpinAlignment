#!/bin/tcsh
starver SL21c
# csh runmakePicoDstQA.csh ../FileLists/19GeV_2019/Phi_embed_test.list 0 TestPhiSwitch 1 4 11 0 0 0.3333 0.0 0.0 0.0 0.0 0.3333 0.0 0.0 0.0 0.0

#root4star -l -b -q runmakePicoDstQA.C\(\"$1\",$2,\"$3\",\"$4\",$5,$6,$7\)
#root4star -l -b -q runmakePicoDstQA_Phi.C\(\"$1\",$2,\"$3\",\"$4\",$5,$6,$7\)
#root4star -l -b -q runmakePicoDstQA_Phi.C\(\"$1\",$2,\"$3\",\"$4\",$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18\)
#root4star -l -b -q runmakePicoDstQA_Phi.C\(\"$1\",$2,\"$3\",\"$4\",$5,$6,$7,$8\)
root4star -l -b -q runmakePicoDstQA_Phi_TTrees.C\(\"$1\",$2,\"$3\",\"$4\",$5,$6,$7,$8\)
#root4star -l -b -q makePicoDstQA_Phi.C\(\"$1\",$2,\"$3\",\"$4\",$5,$6,$7\)

