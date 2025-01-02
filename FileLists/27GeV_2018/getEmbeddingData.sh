#!/bin/bash

filename=Phi_embed.list
rm ${filename}

ls /star/embed/embedding/production_19GeV_2019/Phi_2*_20214008/P21ic.SL21c/2019/*/*/*picoDst* > ${filename}
