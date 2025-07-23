#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_jetpurity()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder + namef_ntuple_jet_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple = (TNtuple*) fpurity->Get((name_ntuple_jetpurity).c_str());

    // Define the necessary histograms to calculate purity
    TH1F* hsig    = new TH1F("hsig"   ,"",nbin_jet_pt+2,unfolding_jet_pt_binning);
    TH1F* hall    = new TH1F("hall"   ,"",nbin_jet_pt+2,unfolding_jet_pt_binning);
    TH1F* hpurity = new TH1F("hpurity","",nbin_jet_pt+2,unfolding_jet_pt_binning);
    hsig->Sumw2();
    hall->Sumw2();
    hpurity->Sumw2();

    set_histogram_style(hpurity, kCyan  , std_line_width, std_marker_style, std_marker_size);
    
    // Project into the histograms
    ntuple->Project("hsig", "jet_pt","jet_pt_truth!=-999");
    ntuple->Project("hall", "jet_pt");
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    hpurity->Divide(hsig,hall,1,1,"B");
    hpurity->Draw();
    hpurity->GetYaxis()->SetRangeUser(0,1.2);
    
    tex->DrawLatexNDC(0.3,0.3,"simulations");

    // c->Print(Form("./plots/jet_purity.pdf"));
}