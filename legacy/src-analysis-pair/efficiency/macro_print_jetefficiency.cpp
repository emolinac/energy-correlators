#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_jetefficiency()
{
    // Open the necessary files
    TFile* fefficiency = new TFile((output_folder + namef_ntuple_jet_efficiency).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple = (TNtuple*) fefficiency->Get((name_ntuple_jetefficiency).c_str());

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig    = new TH1F("hsig"   ,"",nbin_jet_pt+2,unfolding_jet_pt_binning);
    TH1F* hall    = new TH1F("hall"   ,"",nbin_jet_pt+2,unfolding_jet_pt_binning);
    TH1F* hefficiency = new TH1F("hefficiency","",nbin_jet_pt+2,unfolding_jet_pt_binning);
    hsig->Sumw2();
    hall->Sumw2();
    hefficiency->Sumw2();

    set_histogram_style(hefficiency, kCyan  , std_line_width, std_marker_style, std_marker_size);
    
    // Project into the histograms
    ntuple->Project("hsig", "jet_pt_truth","jet_pt!=-999");
    ntuple->Project("hall", "jet_pt_truth");
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    hefficiency->Divide(hsig,hall,1,1,"B");
    hefficiency->Draw();
    hefficiency->GetYaxis()->SetRangeUser(0,1.2);
    
    tex->DrawLatexNDC(0.3,0.3,"simulations");

    // c->Print(Form("./plots/jet_efficiency.pdf"));
}