#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.cpp"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils.cpp"
#include "../../include/utils.h"
#include "../../include/utils-visual.cpp"
#include "../../include/utils-visual.h"

void macro_print_singlepurity_eta_pt(double jet_pt_min_local = jet_pt_min_nom, double jet_pt_max_local = jet_pt_max)
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder + namef_ntuple_reco2truth_singlehadron_match).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_correction_reco).c_str());

    // Define the necessary histograms to calculate purity
    TH2F* hsig    = new TH2F("hsig"   ,"",25,0,50,12,2,4.5);
    TH2F* hall    = new TH2F("hall"   ,"",25,0,50,12,2,4.5);
    TH2F* hpurity = new TH2F("hpurity","",25,0,50,12,2,4.5);
    hsig->Sumw2();
    hall->Sumw2();
    hpurity->Sumw2();

    // Project into the histograms
    ntuple_dtrmatch->Project("hsig","h_eta:h_pt",Form("key_match==1&&jet_pt>%f&&jet_pt<%f",jet_pt_min_local,jet_pt_max_local));
    ntuple_dtrmatch->Project("hall","h_eta:h_pt",Form("jet_pt>%f&&jet_pt<%f",jet_pt_min_local,jet_pt_max_local));
    
    TCanvas* c = new TCanvas("c","",2880,1620);
    c->Draw();

    // PURITY PLOTS
    hpurity->Divide(hsig,hall,1,1,"B");
    
    //set_histogram_style(hpurity, kViolet, std_line_width, std_marker_style, std_marker_size);
    
    hpurity->Draw("coltext");
    //hpurity->GetYaxis()->SetRangeUser(0,1);
    hpurity->SetTitle(";p_{T}(GeV);#eta");

    // gPad->SetLogx(1);

    hpurity->Smooth();
    gStyle->SetPaintTextFormat("4.2f");
    c->Print(Form("./plots/singlehadron_purity_eta_pt_jet_pt%.0fto%.0f.pdf",jet_pt_min_local,jet_pt_max_local));
}