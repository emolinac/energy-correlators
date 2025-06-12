#!/bin/bash

echo "Make sure to execute this code in the home folder."

name_systematic=syst-corr-paradigm

if [ ! -d "./systematics" ]; then
    echo "Systematics folder does not exist. Creating folder ..."
    mkdir ./systematics
fi

if [ ! -d "./systematics/"${name_systematic} ]; then
    echo ${name_systematic}" folder does not exist. Creating folder ..."
    mkdir ./systematics/${name_systematic}

    echo "Copying essential folders ..."
    cp -r bin                  systematics/${name_systematic}
    cp -r include              systematics/${name_systematic}
    cp -r output-files         systematics/${name_systematic}
    cp -r src                  systematics/${name_systematic}
    cp -r src-analysis         systematics/${name_systematic}
    cp -r src-analysis-hadrons systematics/${name_systematic}
    cp -r src-analysis-jets    systematics/${name_systematic}
    cp -r src-analysis-muons   systematics/${name_systematic}

    echo "Copying essential files ..."
    cp Makefile            systematics/${name_systematic}
    cp run_analysis.sh     systematics/${name_systematic}
    cp run_initial_step.sh systematics/${name_systematic}

    echo "Cleaning the house ..."
    cd systematics/${name_systematic}
    make clean
    rm output-files/*.root
    rm src-analysis/plots/*.pdf
    rm src-analysis-hadrons/efficiency/plots/*.pdf
    rm src-analysis-hadrons/purity/plots/*.pdf
    rm src-analysis-hadrons/plots/*.pdf
    rm src-analysis-jets/plots/*.pdf

    echo "Now set all the specs associated to "${name_systematic}
    echo "1.- In directories.h, at the include folder of the systematic, replace the line of the mother folder to:"
    echo "    std::string mother_folder = \"${PWD}/\";"  
    echo "2.- Remember that our default method to do unfolding is 2d"
    
fi