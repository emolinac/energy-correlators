#!/bin/bash

mother_dir=${PWD}
echo $mother_dir

ls -1 systematics/ > folder-list.txt

for line in $(cat ${mother_dir}/folder-list.txt); do
    if [ -d "${mother_dir}/systematics/${line}/include" ]; then
        echo "Include exists on ${line}. Copying binning settings into ${line}"
        cp ${mother_dir}/include/analysis-binning.h ${mother_dir}/systematics/${line}/include/
        # cp ${mother_dir}/include/analysis-constants.h ${mother_dir}/systematics/${line}/include/

        cd ${mother_dir}/systematics/${line}/
        make clean
        make
        cd ${mother_dir}
    fi
done

rm folder-list.txt