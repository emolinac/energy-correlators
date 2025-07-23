#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_singletrackcorrections(double corr_rel_error_local = corr_rel_error, double jet_pt_min_local = jet_pt_binning[0], double jet_pt_max_local = jet_pt_binning[nbin_jet_pt])
{
    gStyle->SetOptStat(1110);
    // Open the necessary files
    TFile*   fcorr       = new TFile((output_folder + namef_ntuple_e2c_corr).c_str()); 
    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    
    TH2F* h = new TH2F("h","",100,0.2,1,100,0.2,1);
    
    // Project into the histograms
    ntuple_data->Project("h","purity:efficiency",Form("efficiency_relerror<%f&&purity_relerror<%f&&jet_pt>%f&&jet_pt<%f",corr_rel_error_local,corr_rel_error_local,jet_pt_min_local,jet_pt_max_local));
    // ntuple_data->Project("h","purity:efficiency","efficiency>0&&purity>0");
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    h->Draw("col");

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    h->SetTitle(";efficiency;purity");
    h->Smooth();

    c->Print(Form("./plots/eff_purity_distribution_jet_pt%.0fto%.0f.pdf",corr_rel_error_local,jet_pt_min_local,jet_pt_max_local));
}