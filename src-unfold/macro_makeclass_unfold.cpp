#include <iostream>
#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"

void macro_makeclass_unfold()
{
    // Add the files located in the input folder
    TChain* chain1 = new TChain((name_ntuple_purity).c_str());
    chain1->Add((output_folder+namef_ntuple_e2c_purity).c_str());
    chain1->MakeClass("TZJetsUnfold");
}

