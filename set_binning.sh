#!/bin/bash

mother_dir=${PWD}
echo $mother_dir

ls -1 systematics/ > folder-list.txt

for line in $(cat folder-list.txt); do
    if [ -d "systematics/${line}/include" ]; then
        echo "Include exists on ${line}. Copying binning settings into ${line}"
        cp include/analysis-binning.h systematics/${line}/include/
    fi
done

rm folder-list.txt