#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

std::string variation_folder = "/home/esteban/Documents/umich-work/lhcb/energy-correlator-inclusivetracks/output-files/";

void macro_print_leptonpair_frac()
{
    TFile* fnominal   = new TFile((output_folder+"ntuple_corre2c.root").c_str());
    TFile* fvariation = new TFile((variation_folder+"ntuple_corre2c.root").c_str());

    TNtuple* ntuple_nominal   = (TNtuple*) fnominal->Get(name_ntuple_data.c_str());
    TNtuple* ntuple_variation = (TNtuple*) fvariation->Get(name_ntuple_data.c_str());

    TH1F* hnominal   = new TH1F("hnominal"  ,"",Nbin_rl,rl_binning);
    TH1F* hvariation = new TH1F("hvariation","",Nbin_rl,rl_binning);
    TH1F* hratio     = new TH1F("hratio"    ,"",Nbin_rl,rl_binning);

    ntuple_nominal->Project("hnominal","R_L",pair_cut);
    ntuple_variation->Project("hvariation","R_L",pair_cut);

    hnominal->Scale(1./hnominal->Integral());
    hvariation->Scale(1./hvariation->Integral());
    
    hratio->Divide(hnominal,hvariation,1,1);
    set_histogram_style(hratio, corr_marker_color_jet_pt[0], std_line_width, corr_marker_style_jet_pt[0], std_marker_size+1);

    hratio->Draw();
    gPad->SetLogx(1);
    hratio->SetTitle(";R_{L};Hadrons/Inclusive");
}