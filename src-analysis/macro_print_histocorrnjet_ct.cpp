#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_histocorrnjet_ct(int niter = nominal_niter)
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + Form("histos_njet_niter%i_ct.root",niter)).c_str(),"RECREATE");
        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_3dpaircorr_rl_jetpt_weightpt_histos_ct).c_str()); 
        if (fcorr->IsZombie()) 
                return;

        TH1F* h_njet           = (TH1F*) fcorr->Get("h_njet");
        TH1F* h_njet_wmuoneff  = (TH1F*) fcorr->Get("h_njet_wmuoneff");
        TH1F* h_efficiency_jet = (TH1F*) fcorr->Get("hefficiency_jet");
        TH1F* h_purity_jet     = (TH1F*) fcorr->Get("hpurity_jet");

        TH1F* h_njet_truth     = (TH1F*) fcorr->Get("h_njet_truth");
        
        // Correct the jets
        TFile* f = new TFile((output_folder + namef_ntuple_jet_purity).c_str());

        TNtuple* ntuple_jet_unfolding = (TNtuple*) f->Get(name_ntuple_jetpurity.c_str());
        
        float jet_pt_unfolding_reco, jet_pt_unfolding_truth;
        set_unfolding_jet_ntuple_branches(ntuple_jet_unfolding, &jet_pt_unfolding_reco, &jet_pt_unfolding_truth);
        
        TH1D* hmeas_jet = new TH1D("hmeas_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1D* htrue_jet = new TH1D("htrue_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2D* hresp_jet = new TH2D("hresp_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        
        for (int evt = 0 ; evt < ntuple_jet_unfolding->GetEntries() ; evt++) {
                ntuple_jet_unfolding->GetEntry(evt);

                if (jet_pt_unfolding_truth != -999)
                        hresp_jet->Fill(jet_pt_unfolding_reco, jet_pt_unfolding_truth);
        }

        RooUnfoldResponse* response_jet = new RooUnfoldResponse(hmeas_jet, htrue_jet, hresp_jet, "response_jet");
        
        TH1F* h_njet_purity_corrected = new TH1F("h_njet_purity_corrected","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        h_njet_purity_corrected->Multiply(h_njet,h_purity_jet,1,1);

        RooUnfoldBayes unfold_jet(response_jet, h_njet_purity_corrected, nominal_niter);

        TH1D* h_njet_unfolded = (TH1D*) unfold_jet.Hunfold();
        
        h_njet_unfolded->Divide(h_efficiency_jet);

        TH1F* h_muon_eff = new TH1F("h_muon_eff","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        h_muon_eff->Divide(h_njet, h_njet_wmuoneff);

        h_njet_unfolded->Divide(h_muon_eff);
        
        fout->cd();
        h_efficiency_jet->Write("hefficiency_jet");
        h_purity_jet->Write("hpurity_jet");
        h_muon_eff->Write("hefficiency_muon");
        h_njet->Write("hnjet_uncorr");
        h_njet_unfolded->Write("hnjet_unfolded");
        h_njet_truth->Write("hnjet_truth");
        gROOT->cd();       

        h_njet_truth->Divide(h_njet_unfolded);

        fout->cd();
        h_njet_truth->Write("truth_to_pseudodata_njet");
        gROOT->cd();       
}
