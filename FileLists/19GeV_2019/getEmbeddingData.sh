#!/bin/bash

filename=Kplus_embed.list
rm ${filename}

ls /star/embed/embedding/production_19GeV_2019/Kplus_2*_20214203/P21ic.SL21c/2019/*/*/*picoDst* > ${filename}
