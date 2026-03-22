#!/bin/bash

particle=Kminus
filename=${particle}_embed.list
rm ${filename}

ls -d -1 /star/data105/embedding/production_7p7GeV_2021/${particle}_400_20231501/P22ib.SL22b/2021/*/*/*picoDst* > ${filename}

for i in {1..9}
do
  ls -d -1 /star/data105/embedding/production_7p7GeV_2021/${particle}_40${i}_20231501/P22ib.SL22b/2021/*/*/*picoDst* >> ${filename}
done

