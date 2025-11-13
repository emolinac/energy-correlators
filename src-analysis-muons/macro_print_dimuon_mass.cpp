#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.cpp"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

void macro_print_dimuon_mass()
{
        gStyle->SetOptStat(1110);
        
        TFile* f = new TFile((input_folder + "Zjet_Data_2016_MU_04212025.root").c_str());

        gDirectory->cd("StdHltZJets");

        TTree* t = (TTree*) gDirectory->Get("DecayTree");

        TH1F* h = new TH1F("h","DiMuon Invariant Mass",100,55,130);
        set_histogram_style(h, 868, std_line_width, std_marker_style, std_marker_size);

        t->Project("h","Z0_M/1000.","");

        TCanvas* c = new TCanvas("c","",1920,1080);
        h->Draw("E1");
        h->SetTitle(";M_{#mu^{+}#mu^{-}} (GeV);");
        c->Print("./plots/dimuon_mass_data.pdf");
}