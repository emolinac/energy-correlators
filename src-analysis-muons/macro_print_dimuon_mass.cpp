#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_dimuon_mass()
{
        TFile* f = new TFile((input_folder + "Zjet_Data_2016_MU_04212025.root").c_str());

        gDirectory->cd("StdHltZJets");

        TTree* t = (TTree*) f->Get("DecayTree");

        TH1F* h = new TH1F("h","DiMuon Invariant Mass",100,50,100);

        t->Project("h","Z0_M","");

        h->Draw("E1");
}