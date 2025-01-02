#!/bin/bash

for i in {0..13}
  do
    root -l -b -q subBackGroundY.C\(4,0,0,0,\"20230627\",${i}\)
    root -l -b -q subBackGroundY.C\(4,0,0,3,\"20230627\",${i}\)
    root -l -b -q subBackGroundY.C\(4,0,0,4,\"20230627\",${i}\)
    root -l -b -q subBackGroundY.C\(4,0,0,5,\"20230627\",${i}\)
done
