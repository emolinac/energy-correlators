#!/bin/bash

mother_dir=${PWD}
echo $mother_dir

ls -1 systematics/ > folder-list.txt

for line in $(cat ${motherdir}/folder-list.txt); do
    if [ -d "${motherdir}/systematics/${line}/include" ]; then
        echo "Include exists on ${line}. Copying binning settings into ${line}"
        cp ${motherdir}/include/analysis-binning.h ${motherdir}/systematics/${line}/include/

        cd ${motherdir}/systematics/${line}/
        make clean
        make
        cd ${motherdir}
    fi
done

rm folder-list.txt