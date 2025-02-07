#!/bin/bash

# Go to src, execute code, and move the respective files
cd ./src
root -q macro_makeclasses.cpp
mv *.C *.h ../include

echo "NOTE : Remember to set the mother and input folders in the directories.h file in the include folder!!"
