#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_histocorreec_rl_jetpt_weightpt_ct(int niter = 8, int niter_jet = 4)
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niter%i_niterjet%i_ct.root",niter, niter_jet)).c_str(),"RECREATE");
        TFile* f    = new TFile((output_folder + namef_ntuple_eec_paircorrections_ct).c_str());
        TFile* fjet = new TFile((output_folder + namef_ntuple_jet_efficiency_ct).c_str());

        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_3dpaircorr_rl_jetpt_weightpt_histos_ct).c_str());
        
        TH3D* h_npair       = (TH3D*) fcorr->Get("h_npair");
        TH3D* h_npair_wmuon = (TH3D*) fcorr->Get("h_npair_wmuon");
        TH3D* h_efficiency  = (TH3D*) fcorr->Get("hefficiency");
        TH3D* h_purity      = (TH3D*) fcorr->Get("hpurity");
        
        TH1F* h_njet           = (TH1F*) fcorr->Get("h_njet");
        TH1F* h_njet_wmuoneff  = (TH1F*) fcorr->Get("h_njet_wmuoneff");
        TH1F* h_efficiency_jet = (TH1F*) fcorr->Get("hefficiency_jet");
        TH1F* h_purity_jet     = (TH1F*) fcorr->Get("hpurity_jet");

        TH3D* h_npair_truth = (TH3D*) fcorr->Get("h_npair_truth");
        TH1F* h_njet_truth  = (TH1F*) fcorr->Get("h_njet_truth");

        // Correct the jets
        TRandom3* rndm = new TRandom3();
        
        TNtuple* ntuple_jet_unfolding = (TNtuple*) fjet->Get(name_ntuple_jetefficiency.c_str());
        
        float jet_pt_unfolding_reco, jet_pt_unfolding_truth;
        set_unfolding_jet_ntuple_branches(ntuple_jet_unfolding, &jet_pt_unfolding_reco, &jet_pt_unfolding_truth);
        
        TH1D* hmeas_jet = new TH1D("hmeas_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1D* htrue_jet = new TH1D("htrue_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2D* hresp_jet = new TH2D("hresp_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        
        for (int evt = 0 ; evt < ntuple_jet_unfolding->GetEntries() ; evt++) {
                ntuple_jet_unfolding->GetEntry(evt);

                if (jet_pt_unfolding_reco != -999)
                        hresp_jet->Fill(jet_pt_unfolding_reco, jet_pt_unfolding_truth);
        }

        RooUnfoldResponse* response_jet = new RooUnfoldResponse(hmeas_jet, htrue_jet, hresp_jet, "response_jet");
        
        TH1F* h_njet_purity_corrected = new TH1F("h_njet_purity_corrected","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        h_njet_purity_corrected->Multiply(h_njet_wmuoneff,h_purity_jet,1,1);

        RooUnfoldBayes unfold_jet(response_jet, h_njet_purity_corrected, niter_jet);

        TH1D* h_njet_unfolded = (TH1D*) unfold_jet.Hreco();

        h_njet_unfolded->Divide(h_efficiency_jet);

        // Correct the npairs
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());
        
        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        
        TH3D* hmeas_npair = new TH3D("hmeas_npair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* htrue_npair = new TH3D("htrue_npair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);

        RooUnfoldResponse* response_npair = new RooUnfoldResponse(hmeas_npair, htrue_npair, "response_npair");
        
        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                ntuple->GetEntry(evt);

                if (R_L_truth != -999)
                        response_npair->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
        }

        TH3D* h_npair_purity_corrected = new TH3D("h_npair_purity_corrected","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        h_npair_purity_corrected->Multiply(h_npair_wmuon,h_purity,1,1);
        
        RooUnfoldBayes unfold_npair(response_npair, h_npair_purity_corrected, niter);

        TH3D* h_npair_unfolded = (TH3D*) unfold_npair.Hreco();
        
        h_npair_unfolded->Divide(h_efficiency);

        for (int i = 1 ; i <= h_npair_unfolded->GetNbinsX(); i++) {
                for (int j = 1 ; j <= h_npair_unfolded->GetNbinsY(); j++) {
                        for (int k = 1 ; k <= h_npair_unfolded->GetNbinsZ() ; k++) {
                                double reweight = h_purity_jet->GetBinContent(j)/h_efficiency_jet->GetBinContent(j);
                                
                                h_npair_unfolded->SetBinContent(i, j, k, h_npair_unfolded->GetBinContent(i, j, k) * reweight);
                                h_npair_unfolded->SetBinError(i, j, k, h_npair_unfolded->GetBinError(i, j, k) * reweight);
                        }
                }
        }

        TH2D* h_eec_unfolded_2d    = new TH2D("h_eec_unfolded_2d"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* h_eec_truth_2d       = new TH2D("h_eec_truth_2d"      ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        apply_unfolded_weights(h_npair_unfolded, h_eec_unfolded_2d);
        apply_unfolded_weights(h_npair_truth   , h_eec_truth_2d);

        TH2D* h_npair_unfolded_2d  = (TH2D*) h_npair_unfolded->Project3D("yx");
        TH2D* h_npair_truth_2d = (TH2D*) h_npair_truth->Project3D("yx");

        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_npair[nbin_jet_pt]; 
        TH1F* hcorr_eec_truth[nbin_jet_pt]; 
        TH1F* hcorr_npair_truth[nbin_jet_pt]; 
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_eec[bin]    = new TH1F(Form("hcorr_eec%i",bin)   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hcorr_npair[bin]  = new TH1F(Form("hcorr_npair%i",bin) , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hcorr_eec[bin]   , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_npair[bin] , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                hcorr_eec_truth[bin]   = new TH1F(Form("hcorr_eec_truth%i",bin)   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hcorr_npair_truth[bin] = new TH1F(Form("hcorr_npair_truth%i",bin) , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hcorr_eec_truth[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_npair_truth[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                for (int i = 1 ; i <= hcorr_eec[bin]->GetNbinsX(); i++) {
                        if (i == 1 || i == hcorr_eec[bin]->GetNbinsX()) {
                                hcorr_eec[bin]->SetBinContent(i, 0);
                                hcorr_eec[bin]->SetBinError(i, 0);
                                hcorr_npair[bin]->SetBinContent(i, 0);
                                hcorr_npair[bin]->SetBinError(i, 0);        

                                hcorr_eec_truth[bin]->SetBinContent(i, 0);
                                hcorr_eec_truth[bin]->SetBinError(i, 0);
                                hcorr_npair_truth[bin]->SetBinContent(i, 0);
                                hcorr_npair_truth[bin]->SetBinError(i, 0);        

                                continue;
                        }

                        hcorr_eec[bin]->SetBinContent(i, h_eec_unfolded_2d->GetBinContent(i, bin + 3));
                        hcorr_eec[bin]->SetBinError(i, h_eec_unfolded_2d->GetBinError(i, bin + 3));
                        hcorr_npair[bin]->SetBinContent(i, h_npair_unfolded_2d->GetBinContent(i, bin + 3));
                        hcorr_npair[bin]->SetBinError(i, h_npair_unfolded_2d->GetBinError(i, bin + 3));

                        hcorr_eec_truth[bin]->SetBinContent(i, h_eec_truth_2d->GetBinContent(i, bin + 3));
                        hcorr_eec_truth[bin]->SetBinError(i, h_eec_truth_2d->GetBinError(i, bin + 3));
                        hcorr_npair_truth[bin]->SetBinContent(i, h_npair_truth_2d->GetBinContent(i, bin + 3));
                        hcorr_npair_truth[bin]->SetBinError(i, h_npair_truth_2d->GetBinError(i, bin + 3));
                }

                hcorr_eec[bin]->Scale(1./h_njet_unfolded->GetBinContent(bin + 3),"width");
                hcorr_npair[bin]->Scale(1./h_njet_unfolded->GetBinContent(bin + 3),"width");

                hcorr_eec_truth[bin]->Scale(1./h_njet_truth->GetBinContent(bin + 3),"width");
                hcorr_npair_truth[bin]->Scale(1./h_njet_truth->GetBinContent(bin + 3),"width");

                fout->cd();
                hcorr_eec[bin]->Write();
                hcorr_npair[bin]->Write();
                hcorr_eec_truth[bin]->Write();
                hcorr_npair_truth[bin]->Write();
                gROOT->cd();

                hcorr_eec[bin]->Divide(hcorr_eec_truth[bin]);
                hcorr_eec[bin]->SetTitle("EEC(pseudodata/truth)");
                hcorr_npair[bin]->Divide(hcorr_npair_truth[bin]);
                hcorr_npair[bin]->SetTitle("Npair(pseudodata/truth)");
                
                fout->cd();
                hcorr_eec[bin]->Write(Form("pseudodata_to_truth_eec%i",bin));
                hcorr_npair[bin]->Write(Form("pseudodata_to_truth_npair%i",bin));
                gROOT->cd();       
        }
}
