#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_histocorreec_rl_jetpt_weightpt_shapect(int niter = 4, int niter_jet = 4)
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niter%i_niterjet%i_shapect.root",niter, niter_jet)).c_str(),"RECREATE");
        TFile* f    = new TFile((output_folder + namef_ntuple_reco2truth_match_ct).c_str());
        TFile* fjet = new TFile((output_folder + namef_ntuple_truth2reco_match_ct).c_str());

        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_3dpaircorr_rl_jetpt_weightpt_histos_ct).c_str());
        TFile* fcorr_data = new TFile((output_folder + namef_3dpaircorr_rl_jetpt_weightpt_histos).c_str());

        //  Pseudodata
        TH3D* h_npair      = (TH3D*) fcorr->Get("h_npair_wmuon");
        TH3D* h_eqchnpair  = (TH3D*) fcorr->Get("h_eqchnpair_wmuon");
        TH3D* h_neqchnpair = (TH3D*) fcorr->Get("h_neqchnpair_wmuon");
        TH3D* h_efficiency = (TH3D*) fcorr->Get("hefficiency");
        TH3D* h_purity     = (TH3D*) fcorr->Get("hpurity");
        
        TH1F* h_njet           = (TH1F*) fcorr->Get("h_njet");
        TH1F* h_efficiency_jet = (TH1F*) fcorr->Get("hefficiency_jet");
        TH1F* h_purity_jet     = (TH1F*) fcorr->Get("hpurity_jet");

        // Truth
        TH3D* h_npair_truth      = (TH3D*) fcorr->Get("h_npair_truth");
        TH3D* h_eqchnpair_truth  = (TH3D*) fcorr->Get("h_eqchnpair_truth");
        TH3D* h_neqchnpair_truth = (TH3D*) fcorr->Get("h_neqchnpair_truth");
        TH1F* h_njet_truth       = (TH1F*) fcorr->Get("h_njet_truth");

        TH2D* h_eec_truth_2d      = new TH2D("h_eec_truth_2d"     ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* h_eqcheec_truth_2d  = new TH2D("h_eqcheec_truth_2d" ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* h_neqcheec_truth_2d = new TH2D("h_neqcheec_truth_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        apply_unfolded_weights(h_npair_truth     , h_eec_truth_2d);
        apply_unfolded_weights(h_eqchnpair_truth , h_eqcheec_truth_2d);
        apply_unfolded_weights(h_neqchnpair_truth, h_neqcheec_truth_2d);

        // Real data
        TH3D* h_npair_data_wmuon      = (TH3D*) fcorr_data->Get("h_npair_wmuon");
        TH3D* h_eqchnpair_data_wmuon  = (TH3D*) fcorr_data->Get("h_eqchnpair_wmuon");
        TH3D* h_neqchnpair_data_wmuon = (TH3D*) fcorr_data->Get("h_neqchnpair_wmuon");

        // Correct the jets
        TRandom3* rndm = new TRandom3(0);
        
        TNtuple* ntuple_jet_unfolding = (TNtuple*) fjet->Get(name_ntuple_jet_truth2reco_match.c_str());
        
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
        h_njet_purity_corrected->Multiply(h_njet,h_purity_jet,1,1);

        RooUnfoldBayes unfold_jet(response_jet, h_njet_purity_corrected, niter_jet);

        TH1D* h_njet_unfolded = (TH1D*) unfold_jet.Hreco();

        h_njet_unfolded->Divide(h_efficiency_jet);

        // Correct the npairs
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());
        
        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth, eq_charge_reco;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        ntuple->SetBranchAddress("eq_charge", &eq_charge_reco);

        TH3D* hmeas_npair = new TH3D("hmeas_npair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* htrue_npair = new TH3D("htrue_npair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);

        RooUnfoldResponse* response_npair = new RooUnfoldResponse(hmeas_npair, htrue_npair, "response_npair");
        
        TH3D* hmeas_eqchnpair = new TH3D("hmeas_eqchnpair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* htrue_eqchnpair = new TH3D("htrue_eqchnpair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        RooUnfoldResponse* response_eqchnpair = new RooUnfoldResponse(hmeas_eqchnpair, htrue_eqchnpair, "response_eqchnpair");

        TH3D* hmeas_neqchnpair = new TH3D("hmeas_neqchnpair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* htrue_neqchnpair = new TH3D("htrue_neqchnpair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        RooUnfoldResponse* response_neqchnpair = new RooUnfoldResponse(hmeas_neqchnpair, htrue_neqchnpair, "response_neqchnpair");

        // Create the reweighting matrix
        TH3D* h_npair_reweight      = new TH3D("h_npair_reweight" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* h_eqchnpair_reweight  = new TH3D("h_eqchnpair_reweight" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* h_neqchnpair_reweight = new TH3D("h_neqchnpair_reweight" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);

        TH3D* h_npair_norm           = (TH3D*) h_npair->Clone("h_npair_norm");
        TH3D* h_eqchnpair_norm       = (TH3D*) h_eqchnpair->Clone("h_eqchnpair_norm");
        TH3D* h_neqchnpair_norm      = (TH3D*) h_neqchnpair->Clone("h_neqchnpair_norm");

        TH3D* h_npair_data_wmuon_norm      = (TH3D*) h_npair_data_wmuon->Clone("h_npair_data_wmuon_norm");
        TH3D* h_eqchnpair_data_wmuon_norm  = (TH3D*) h_eqchnpair_data_wmuon->Clone("h_eqchnpair_data_wmuon_norm");
        TH3D* h_neqchnpair_data_wmuon_norm = (TH3D*) h_neqchnpair_data_wmuon->Clone("h_neqchnpair_data_wmuon_norm");

        h_npair_norm->Scale(1./h_npair_norm->Integral());
        h_eqchnpair_norm->Scale(1./h_eqchnpair_norm->Integral());
        h_neqchnpair_norm->Scale(1./h_neqchnpair_norm->Integral());
        h_npair_data_wmuon_norm->Scale(1./h_npair_data_wmuon_norm->Integral());
        h_eqchnpair_data_wmuon_norm->Scale(1./h_eqchnpair_data_wmuon_norm->Integral());
        h_neqchnpair_data_wmuon_norm->Scale(1./h_neqchnpair_data_wmuon_norm->Integral());

        h_npair_reweight->Divide(h_npair_norm,h_npair_data_wmuon_norm);
        h_eqchnpair_reweight->Divide(h_eqchnpair_norm,h_eqchnpair_data_wmuon_norm);
        h_neqchnpair_reweight->Divide(h_neqchnpair_norm,h_neqchnpair_data_wmuon_norm);

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        gStyle->SetOptStat("");
        gStyle->SetPaintTextFormat("4.2f");        
        
        TLatex latex;
        latex.SetTextAlign(22); // center alignment
        latex.SetTextSize(text_size_correction_plots);
        latex.SetTextColor(kBlack);

        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                ntuple->GetEntry(evt);

                double reweigh_npair      = h_npair_reweight->GetBinContent(h_npair_reweight->FindBin(R_L_reco, jet_pt_reco, weight_pt_reco));
                double reweigh_echnpair   = h_eqchnpair_reweight->GetBinContent(h_eqchnpair_reweight->FindBin(R_L_reco, jet_pt_reco, weight_pt_reco));
                double reweigh_neqchnpair = h_neqchnpair_reweight->GetBinContent(h_neqchnpair_reweight->FindBin(R_L_reco, jet_pt_reco, weight_pt_reco));

                if (R_L_truth != -999)
                        response_npair->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth, reweigh_npair);
                if (R_L_truth != -999 && eq_charge_reco > 0)
                        response_eqchnpair->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth, reweigh_echnpair);
                if (R_L_truth != -999 && eq_charge_reco < 0)
                        response_neqchnpair->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth, reweigh_neqchnpair);
        }

        TH3F* h_npair_purity_corrected      = new TH3F("h_npair_purity_corrected",     "",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_eqchnpair_purity_corrected  = new TH3F("h_eqchnpair_purity_corrected", "",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_neqchnpair_purity_corrected = new TH3F("h_neqchnpair_purity_corrected","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        
        h_npair_purity_corrected->Multiply(h_npair,h_purity,1,1);
        h_eqchnpair_purity_corrected->Multiply(h_eqchnpair,h_purity,1,1);
        h_neqchnpair_purity_corrected->Multiply(h_neqchnpair,h_purity,1,1);
        
        RooUnfoldBayes unfold_npair(response_npair, h_npair_purity_corrected, niter);
        RooUnfoldBayes unfold_eqchnpair(response_eqchnpair, h_eqchnpair_purity_corrected, niter);
        RooUnfoldBayes unfold_neqchnpair(response_neqchnpair, h_neqchnpair_purity_corrected, niter);

        TH3D* h_npair_unfolded      = (TH3D*) unfold_npair.Hreco();
        TH3D* h_eqchnpair_unfolded  = (TH3D*) unfold_eqchnpair.Hreco();
        TH3D* h_neqchnpair_unfolded = (TH3D*) unfold_neqchnpair.Hreco();

        h_npair_unfolded->Divide(h_efficiency);
        h_eqchnpair_unfolded->Divide(h_efficiency);
        h_neqchnpair_unfolded->Divide(h_efficiency);

        apply_jet_weight_to_npairs(h_npair_unfolded, h_purity_jet, h_efficiency_jet);
        apply_jet_weight_to_npairs(h_eqchnpair_unfolded, h_purity_jet, h_efficiency_jet);
        apply_jet_weight_to_npairs(h_neqchnpair_unfolded, h_purity_jet, h_efficiency_jet);

        TH2D* h_eec_unfolded_2d      = new TH2D("h_eec_unfolded_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* h_eqcheec_unfolded_2d  = new TH2D("h_eqcheec_unfolded_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* h_neqcheec_unfolded_2d = new TH2D("h_neqcheec_unfolded_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        apply_unfolded_weights(h_npair_unfolded, h_eec_unfolded_2d);
        apply_unfolded_weights(h_eqchnpair_unfolded, h_eqcheec_unfolded_2d);
        apply_unfolded_weights(h_neqchnpair_unfolded, h_neqcheec_unfolded_2d);
        
        TH2D* h_npair_unfolded_2d  = (TH2D*) h_npair_unfolded->Project3D("yx");
        TH2D* h_npair_truth_2d = (TH2D*) h_npair_truth->Project3D("yx");

        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_eqcheec[nbin_jet_pt]; 
        TH1F* hcorr_neqcheec[nbin_jet_pt]; 
        TH1F* hcorr_tau[nbin_jet_pt]; 
        
        TH1F* hcorr_eec_truth[nbin_jet_pt]; 
        TH1F* hcorr_eqcheec_truth[nbin_jet_pt]; 
        TH1F* hcorr_neqcheec_truth[nbin_jet_pt]; 
        TH1F* hcorr_tau_truth[nbin_jet_pt]; 

        // Define unity
        TH1F* h_unity_uounderflow = new TH1F("h_unity_uounderflow","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
        h_unity_uounderflow->Sumw2();
        set_unity_content(h_unity_uounderflow);

        TH1F* h_unity_taubinning[nbin_jet_pt];

        double tau_binning[nbin_jet_pt][nbin_rl_nominal + 1];
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                double avge_pt2_jet = (jet_pt_binning[bin + 1] + jet_pt_binning[bin])/2.;
                
                get_tau_binning_from_eec_binning(tau_binning[bin], rl_nominal_binning, avge_pt2_jet);
                
                int nominal_jet_pt_bin = bin + 3;

                h_unity_taubinning[bin] = new TH1F(Form("h_unity%i",bin),"",nbin_rl_nominal,tau_binning[bin]);
                h_unity_taubinning[bin]->Sumw2();
                set_unity_content(h_unity_taubinning[bin]);

                // Pseudodata operations
                hcorr_eec[bin]      = new TH1F(Form("hcorr_eec%i",bin)     , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hcorr_eqcheec[bin]  = new TH1F(Form("hcorr_eqcheec%i",bin) , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hcorr_neqcheec[bin] = new TH1F(Form("hcorr_neqcheec%i",bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hcorr_tau[bin]      = new TH1F(Form("hcorr_tau%i",bin)     , "", nbin_rl_nominal,tau_binning[bin]);
                
                set_histogram_style(hcorr_eec[bin]     , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_eqcheec[bin] , corr_marker_color_jet_pt[bin], std_line_width-1, std_marker_style_jet_pt[bin] , std_marker_size+1);
                set_histogram_style(hcorr_neqcheec[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_tau[bin]   , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h_eec_unfolded_2d     , hcorr_eec[bin]     , nominal_jet_pt_bin);
                project_nominal_phase_space(h_eqcheec_unfolded_2d , hcorr_eqcheec[bin] , nominal_jet_pt_bin);
                project_nominal_phase_space(h_neqcheec_unfolded_2d, hcorr_neqcheec[bin], nominal_jet_pt_bin);
                
                hcorr_eqcheec[bin]->Divide(hcorr_eec[bin]);
                hcorr_neqcheec[bin]->Divide(hcorr_eec[bin]);

                hcorr_eec[bin]->Scale(1./h_njet_unfolded->GetBinContent(bin + 3),"width");
                
                get_tau_from_uoflow_eec(hcorr_eec[bin], hcorr_tau[bin], avge_pt2_jet);
                
                // Truth operations
                hcorr_eec_truth[bin]      = new TH1F(Form("hcorr_eec_truth%i",bin)     , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hcorr_eqcheec_truth[bin]  = new TH1F(Form("hcorr_eqcheec_truth%i",bin) , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hcorr_neqcheec_truth[bin] = new TH1F(Form("hcorr_neqcheec_truth%i",bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hcorr_tau_truth[bin]      = new TH1F(Form("hcorr_tau_truth%i",bin)     , "", nbin_rl_nominal,tau_binning[bin]);;

                project_nominal_phase_space(h_eec_truth_2d     , hcorr_eec_truth[bin]     , nominal_jet_pt_bin);
                project_nominal_phase_space(h_eqcheec_truth_2d , hcorr_eqcheec_truth[bin] , nominal_jet_pt_bin);
                project_nominal_phase_space(h_neqcheec_truth_2d, hcorr_neqcheec_truth[bin], nominal_jet_pt_bin);

                hcorr_eqcheec_truth[bin]->Divide(hcorr_eec_truth[bin]);
                hcorr_neqcheec_truth[bin]->Divide(hcorr_eec_truth[bin]);

                hcorr_eec_truth[bin]->Scale(1./h_njet_truth->GetBinContent(bin + 3),"width");
                
                get_tau_from_uoflow_eec(hcorr_eec_truth[bin], hcorr_tau_truth[bin], avge_pt2_jet);

                fout->cd();
                hcorr_eec[bin]->Write();
                hcorr_tau[bin]->Write();
                hcorr_eec_truth[bin]->Write();
                hcorr_tau_truth[bin]->Write();
                hcorr_eqcheec[bin]->Write();
                hcorr_neqcheec[bin]->Write();
                hcorr_eqcheec_truth[bin]->Write();
                hcorr_neqcheec_truth[bin]->Write();
                gROOT->cd();

                hcorr_eec[bin]->Divide(hcorr_eec_truth[bin]);
                hcorr_eec[bin]->SetTitle("EEC(pseudodata/truth)");
                hcorr_tau[bin]->Divide(hcorr_tau_truth[bin]);
                hcorr_tau[bin]->SetTitle("Npair(pseudodata/truth)");
                hcorr_eqcheec[bin]->Divide(hcorr_eqcheec_truth[bin]);
                hcorr_eqcheec[bin]->SetTitle("Eq. Charged EEC(pseudodata/truth)");
                hcorr_neqcheec[bin]->Divide(hcorr_neqcheec_truth[bin]);
                hcorr_neqcheec[bin]->SetTitle("Op. Charged EEC(pseudodata/truth)");
                
                hcorr_eec[bin]->Add(h_unity_uounderflow, -1);
                hcorr_tau[bin]->Add(h_unity_taubinning[bin], -1);
                hcorr_eqcheec[bin]->Add(h_unity_uounderflow, -1);
                hcorr_neqcheec[bin]->Add(h_unity_uounderflow, -1);

                hcorr_eec[bin]->Multiply(hcorr_eec[bin]);
                hcorr_tau[bin]->Multiply(hcorr_tau[bin]);
                hcorr_eqcheec[bin]->Multiply(hcorr_eqcheec[bin]);
                hcorr_neqcheec[bin]->Multiply(hcorr_neqcheec[bin]);

                square_root_bins(hcorr_eec[bin]);
                square_root_bins(hcorr_tau[bin]);
                square_root_bins(hcorr_eqcheec[bin]);
                square_root_bins(hcorr_neqcheec[bin]);

                fout->cd();
                hcorr_eec[bin]->Write(Form("relerror_eec%i",bin));
                hcorr_tau[bin]->Write(Form("relerror_tau%i",bin));
                hcorr_eqcheec[bin]->Write(Form("relerror_eqcheec%i",bin));
                hcorr_neqcheec[bin]->Write(Form("relerror_neqcheec%i",bin));
                gROOT->cd();       
        }

        // Print all the relevant ratios
        c = new TCanvas("c","",1800,600);
        c->Draw();
        c->Divide(3,1);
        
        TLatex* tex = new TLatex();
        tex->SetTextColorAlpha(16,0.3);
        tex->SetTextSize(0.1991525);
        tex->SetTextAngle(26.15998);
        tex->SetLineWidth(2);

        THStack* s[3];
        TLegend* l[3];

        TLine* line = new TLine(unfolding_rl_nominal_binning[1], 1, unfolding_rl_nominal_binning[nbin_rl_nominal_unfolding-1], 1);

        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                c->cd(bin+1);
                s[bin] = new THStack();
                l[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,gPad->GetBottomMargin()+0.01,1-gPad->GetRightMargin()-0.01,gPad->GetBottomMargin()+0.12);

                s[bin]->Add(hcorr_eec[bin],"PE1");
                
                s[bin]->Draw("NOSTACK");
                
                s[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};Corr. Pseudodata / Truth",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s[bin]->SetMaximum(0.4);
                s[bin]->SetMinimum(0.);
                s[bin]->GetXaxis()->SetRangeUser(unfolding_rl_nominal_binning[1],unfolding_rl_nominal_binning[nbin_rl_nominal_unfolding-1]);

                l[bin]->AddEntry(hcorr_eec[bin], "EEC", "p");
                
                gPad->SetLogx(1);
                
                l[bin]->Draw("SAME");    
                line->Draw("SAME");
        }

        c->Print("./plots/closure-test-shape-eec.pdf");

        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                c->cd(bin+1);
                s[bin] = new THStack();
                l[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,gPad->GetBottomMargin()+0.01,1-gPad->GetRightMargin()-0.01,gPad->GetBottomMargin()+0.12);

                s[bin]->Add(hcorr_eqcheec[bin],"PE1");
                s[bin]->Add(hcorr_neqcheec[bin],"PE1");
                
                s[bin]->Draw("NOSTACK");
                
                s[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};Corr. Pseudodata / Truth",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s[bin]->SetMaximum(0.4);
                s[bin]->SetMinimum(0.);
                s[bin]->GetXaxis()->SetRangeUser(unfolding_rl_nominal_binning[1],unfolding_rl_nominal_binning[nbin_rl_nominal_unfolding-1]);

                l[bin]->AddEntry(hcorr_eqcheec[bin], "Eq. Ch. EEC", "p");
                l[bin]->AddEntry(hcorr_neqcheec[bin], "Op. Ch. EEC", "p");
                
                gPad->SetLogx(1);
                
                l[bin]->Draw("SAME");    
                line->Draw("SAME");
        }

        c->Print("./plots/closure-test-shape-chargedeec.pdf");
}

