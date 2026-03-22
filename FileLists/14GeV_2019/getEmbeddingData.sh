#!/bin/bash

particle=$1
request=$2

filename=${particle}_embed.list
rm ${filename}

ls /star/data105/embedding/production_14p5GeV_2019/${particle}_3*_${request}/P21ic.SL21c/2019/*/*/*picoDst* > ${filename}
