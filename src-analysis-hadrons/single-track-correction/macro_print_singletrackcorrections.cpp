#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils.h"
#include "../../include/utils-visual.h"

void macro_print_singletrackcorrections(double jet_pt_min_local = jet_pt_binning[0], double jet_pt_max_local = jet_pt_binning[nbin_jet_pt])
{
    gStyle->SetOptStat(1110);
    // Open the necessary files
    TFile*   fcorr       = new TFile((output_folder + namef_ntuple_eec_corr).c_str()); 
    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    
    TH2F* h = new TH2F("h","",150,0.2,1,150,0.2,1);
    
    // Project into the histograms
    ntuple_data->Project("h","purity:efficiency",Form("jet_pt>%f&&jet_pt<%f",jet_pt_min_local,jet_pt_max_local));
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    h->Draw("col");

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    h->SetTitle(";efficiency;purity");
    h->Smooth();

    c->Print(Form("./plots/eff_purity_distribution_jet_pt%.0fto%.0f.pdf",jet_pt_min_local,jet_pt_max_local));
}