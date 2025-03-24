#include "../../include/analysis-constants.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_singleefficiency_eta_momentum(double jet_pt_min_local = jet_pt_min_nom, double jet_pt_max_local = jet_pt_max)
{
// Open the necessary files
    TFile* fefficiency = new TFile((output_folder+namef_ntuple_e2c_efficiency).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_mcreco = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_reco).c_str());
    TNtuple* ntuple_mc     = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_mc).c_str());

    // Define the necessary histograms to calculate efficiency
    //TH2F* hsig        = new TH2F("hsig"       ,"",ndim_corr,track_p_min,track_p_max,ndim_corr,eta_min,eta_max);
    //TH2F* hall        = new TH2F("hall"       ,"",ndim_corr,track_p_min,track_p_max,ndim_corr,eta_min,eta_max);
    //TH2F* hefficiency = new TH2F("hefficiency","",ndim_corr,track_p_min,track_p_max,ndim_corr,eta_min,eta_max);
    TH2F* hsig        = new TH2F("hsig"       ,"",ic_p_nbins,ic_p_binning,12,2,4.5);
    TH2F* hall        = new TH2F("hall"       ,"",ic_p_nbins,ic_p_binning,12,2,4.5);
    TH2F* hefficiency = new TH2F("hefficiency","",ic_p_nbins,ic_p_binning,12,2,4.5);
    hsig->Sumw2();
    hall->Sumw2();
    hefficiency->Sumw2();

    // Project into the histograms
    ntuple_mcreco->Project("hsig","h_eta_truth:h_p_truth",single_signal_cut+Form("jet_pt>%f&&jet_pt<%f",jet_pt_min_local,jet_pt_max_local));
    ntuple_mc->Project("hall","h_eta:h_p",pair_cut             +Form("jet_pt>%f&&jet_pt<%f",jet_pt_min_local,jet_pt_max_local));

    TCanvas* c = new TCanvas("c","",2880,1620);
    c->Draw();

    // efficiency PLOTS
    hefficiency->Divide(hsig,hall,1,1,"B");

    for(int bin_x = 1 ; bin_x <= ndim_corr ; bin_x++)
    {
        for(int bin_y = 1 ; bin_y <= ndim_corr ; bin_y++)
        {
            int bin = hefficiency->GetBin(bin_x,bin_y);
            if(hefficiency->GetBinContent(bin)>1) hefficiency->SetBinContent(bin,0);
        }
    }
        
    hefficiency->Draw("coltext");
    hefficiency->SetTitle(";p(GeV);#eta");

    gPad->SetLogx(1);

    hefficiency->Smooth();
    gStyle->SetPaintTextFormat("4.2f");
    c->Print(Form("./plots/singlehadron_efficiency_eta_momentum_jetpt%.0fto%.0f.pdf",jet_pt_min_local,jet_pt_max_local));
}