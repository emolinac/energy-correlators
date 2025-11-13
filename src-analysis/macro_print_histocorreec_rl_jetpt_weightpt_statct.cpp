#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.cpp"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

void macro_print_histocorreec_rl_jetpt_weightpt_statct(int niter = 4, int niter_jet = 4)
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niter%i_niterjet%i_statct_niterct%i.root",niter, niter_jet,niter_ct)).c_str(),"RECREATE");
        TFile* f    = new TFile((output_folder + namef_ntuple_reco2truth_match_ct).c_str());
        TFile* fjet = new TFile((output_folder + namef_ntuple_truth2reco_match_ct).c_str());

        gROOT->cd();

        TFile* fcorr      = new TFile((output_folder + namef_3dpaircorr_rl_jetpt_weightpt_histos_ct).c_str());
        TFile* fcorr_data = new TFile((output_folder + namef_3dpaircorr_rl_jetpt_weightpt_histos).c_str());

        // Pseudodata
        TH3D* h_npair            = (TH3D*) fcorr->Get("h_npair");
        TH3D* h_npair_wmuon      = (TH3D*) fcorr->Get("h_npair_wmuon");
        TH3D* h_eqchnpair_wmuon  = (TH3D*) fcorr->Get("h_eqchnpair_wmuon");
        TH3D* h_neqchnpair_wmuon = (TH3D*) fcorr->Get("h_neqchnpair_wmuon");
        TH3D* h_efficiency       = (TH3D*) fcorr->Get("hefficiency");
        TH3D* h_purity           = (TH3D*) fcorr->Get("hpurity");
        
        TH1F* h_njet           = (TH1F*) fcorr->Get("h_njet");
        TH1F* h_njet_wmuoneff  = (TH1F*) fcorr->Get("h_njet_wmuoneff");
        TH1F* h_efficiency_jet = (TH1F*) fcorr->Get("hefficiency_jet");
        TH1F* h_purity_jet     = (TH1F*) fcorr->Get("hpurity_jet");

        // Make N copies in order to not spoil the CT
        TH3D* h_npair_wmuon_sample[niter_ct];
        TH3D* h_eqchnpair_wmuon_sample[niter_ct];
        TH3D* h_neqchnpair_wmuon_sample[niter_ct];
        
        // Truth
        TH3D* h_npair_truth      = (TH3D*) fcorr->Get("h_npair_truth");
        TH3D* h_eqchnpair_truth  = (TH3D*) fcorr->Get("h_eqchnpair_truth");
        TH3D* h_neqchnpair_truth = (TH3D*) fcorr->Get("h_neqchnpair_truth");
        TH1F* h_njet_truth       = (TH1F*) fcorr->Get("h_njet_truth");

        // Real data
        TH3D* h_npair_data_wmuon      = (TH3D*) fcorr_data->Get("h_npair_wmuon");
        TH3D* h_eqchnpair_data_wmuon  = (TH3D*) fcorr_data->Get("h_eqchnpair_wmuon");
        TH3D* h_neqchnpair_data_wmuon = (TH3D*) fcorr_data->Get("h_neqchnpair_wmuon");
        
        // Truth operations
        TH2D* h_eec_truth_2d         = new TH2D("h_eec_truth_2d"     ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* h_eqcheec_truth_2d     = new TH2D("h_eqcheec_truth_2d" ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* h_neqcheec_truth_2d    = new TH2D("h_neqcheec_truth_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);

        apply_unfolded_weights(h_npair_truth     , h_eec_truth_2d);
        apply_unfolded_weights(h_eqchnpair_truth , h_eqcheec_truth_2d);
        apply_unfolded_weights(h_neqchnpair_truth, h_neqcheec_truth_2d);
        
        TH2D* h_npair_truth_2d = (TH2D*) h_npair_truth->Project3D("yx");
        TH1F* hcorr_eec_truth[nbin_jet_pt]; 
        TH1F* hcorr_eqcheec_truth[nbin_jet_pt]; 
        TH1F* hcorr_neqcheec_truth[nbin_jet_pt]; 
        TH1F* hcorr_tau_truth[nbin_jet_pt]; 

        double tau_binning[nbin_jet_pt][nbin_rl_nominal + 1];
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                double avge_pt2_jet = (jet_pt_binning[bin + 1] + jet_pt_binning[bin])/2.;
                
                get_tau_binning_from_eec_binning(tau_binning[bin], rl_nominal_binning, avge_pt2_jet);
                
                int nominal_jet_pt_bin = bin + 3;
                
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
                
        }
        
        // Pseudodata operations
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
        h_njet_purity_corrected->Multiply(h_njet_wmuoneff,h_purity_jet,1,1);

        RooUnfoldBayes unfold_jet(response_jet, h_njet_purity_corrected, niter_jet);

        TH1D* h_njet_unfolded = (TH1D*) unfold_jet.Hreco();

        h_njet_unfolded->Divide(h_efficiency_jet);

        // Correct the npairs
        TH3D* hmeas_npair      = new TH3D("hmeas_npair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* htrue_npair      = new TH3D("htrue_npair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* hmeas_eqchnpair  = new TH3D("hmeas_eqchnpair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* htrue_eqchnpair  = new TH3D("htrue_eqchnpair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* hmeas_neqchnpair = new TH3D("hmeas_neqchnpair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* htrue_neqchnpair = new TH3D("htrue_neqchnpair" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        
        TH3F* h_npair_purity_corrected      = new TH3F("h_npair_purity_corrected",     "",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_eqchnpair_purity_corrected  = new TH3F("h_eqchnpair_purity_corrected", "",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_neqchnpair_purity_corrected = new TH3F("h_neqchnpair_purity_corrected","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);

        TH3D* h_npair_unfolded[niter_ct];
        TH3D* h_eqchnpair_unfolded[niter_ct];
        TH3D* h_neqchnpair_unfolded[niter_ct];

        TH2D* h_eec_unfolded_2d      = new TH2D("h_eec_unfolded_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* h_eqcheec_unfolded_2d  = new TH2D("h_eqcheec_unfolded_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* h_neqcheec_unfolded_2d = new TH2D("h_neqcheec_unfolded_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        TH2D* h_npair_unfolded_2d[niter_ct];

        TH1F* hcorr_eec[nbin_jet_pt][niter_ct]; 
        TH1F* hcorr_eqcheec[nbin_jet_pt][niter_ct]; 
        TH1F* hcorr_neqcheec[nbin_jet_pt][niter_ct]; 
        TH1F* hcorr_tau[nbin_jet_pt][niter_ct]; 
        
        TH1F* hcorr_eec_total[nbin_jet_pt]; 
        TH1F* hcorr_eqcheec_total[nbin_jet_pt]; 
        TH1F* hcorr_neqcheec_total[nbin_jet_pt]; 
        TH1F* hcorr_tau_total[nbin_jet_pt]; 

        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                hcorr_eec_total[i]      = new TH1F(Form("hcorr_eec_total_%i",i)     ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hcorr_eqcheec_total[i]  = new TH1F(Form("hcorr_eqcheec_total_%i",i) ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hcorr_neqcheec_total[i] = new TH1F(Form("hcorr_neqcheec_total_%i",i),"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hcorr_tau_total[i]      = new TH1F(Form("hcorr_tau_total_%i",i)     ,"",nbin_rl_nominal,tau_binning[i]);
        }
        
        RooUnfoldResponse* response_npair[niter_ct];
        RooUnfoldResponse* response_eqchnpair[niter_ct];
        RooUnfoldResponse* response_neqchnpair[niter_ct];
        RooUnfoldBayes* unfold_npair[niter_ct];
        RooUnfoldBayes* unfold_eqchnpair[niter_ct];
        RooUnfoldBayes* unfold_neqchnpair[niter_ct];
                        
        TH1F* h_unity_uounderflow = new TH1F("h_unity_uounderflow","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
        TH1F* h_unity_taubinning[nbin_jet_pt];
        h_unity_taubinning[0] = new TH1F(Form("h_unity%i",0),"",nbin_rl_nominal,tau_binning[0]);
        h_unity_taubinning[1] = new TH1F(Form("h_unity%i",1),"",nbin_rl_nominal,tau_binning[1]);
        h_unity_taubinning[2] = new TH1F(Form("h_unity%i",2),"",nbin_rl_nominal,tau_binning[2]);
        
        h_unity_uounderflow->Sumw2();
        h_unity_taubinning[0]->Sumw2();
        h_unity_taubinning[1]->Sumw2();
        h_unity_taubinning[2]->Sumw2();

        set_unity_content(h_unity_uounderflow);
        set_unity_content(h_unity_taubinning[0]);
        set_unity_content(h_unity_taubinning[1]);
        set_unity_content(h_unity_taubinning[2]);

        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());
        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth, eq_charge_reco;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        ntuple->SetBranchAddress("eq_charge", &eq_charge_reco);
        // Pseudodata operations
        for (int ct_iter = 0 ; ct_iter < niter_ct ; ct_iter++) {
                response_npair[ct_iter]      = new RooUnfoldResponse(hmeas_npair     , htrue_npair     , Form("response_npair_%i"     ,ct_iter));
                response_eqchnpair[ct_iter]  = new RooUnfoldResponse(hmeas_eqchnpair , htrue_eqchnpair , Form("response_eqchnpair_%i" ,ct_iter));
                response_neqchnpair[ct_iter] = new RooUnfoldResponse(hmeas_neqchnpair, htrue_neqchnpair, Form("response_neqchnpair_%i",ct_iter));
        
                for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                        ntuple->GetEntry(evt);

                        if (R_L_truth != -999)
                                response_npair[ct_iter]->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
                        if (R_L_truth != -999 && eq_charge_reco > 0)
                                response_eqchnpair[ct_iter]->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
                        if (R_L_truth != -999 && eq_charge_reco < 0)
                                response_neqchnpair[ct_iter]->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
                }

                h_npair_wmuon_sample[ct_iter]      = (TH3D*) h_npair_wmuon->Clone(Form("h_npair_wmuon_sample%i",ct_iter));
                h_eqchnpair_wmuon_sample[ct_iter]  = (TH3D*) h_eqchnpair_wmuon->Clone(Form("h_eqchnpair_wmuon_sample%i",ct_iter));
                h_neqchnpair_wmuon_sample[ct_iter] = (TH3D*) h_neqchnpair_wmuon->Clone(Form("h_neqchnpair_wmuon_sample%i",ct_iter));

                // Smear the pseudodata
                smear_pseudodata(h_npair_wmuon_sample[ct_iter]     , h_npair_data_wmuon     , rndm);
                smear_pseudodata(h_eqchnpair_wmuon_sample[ct_iter] , h_eqchnpair_data_wmuon , rndm);
                smear_pseudodata(h_neqchnpair_wmuon_sample[ct_iter], h_neqchnpair_data_wmuon, rndm);
                
                // Correct the pseudodata
                h_npair_purity_corrected->Multiply(h_npair_wmuon_sample[ct_iter],h_purity,1,1);
                h_eqchnpair_purity_corrected->Multiply(h_eqchnpair_wmuon_sample[ct_iter],h_purity,1,1);
                h_neqchnpair_purity_corrected->Multiply(h_neqchnpair_wmuon_sample[ct_iter],h_purity,1,1);
                
                unfold_npair[ct_iter]      = new RooUnfoldBayes(response_npair[ct_iter], h_npair_purity_corrected, niter);
                unfold_eqchnpair[ct_iter]  = new RooUnfoldBayes(response_eqchnpair[ct_iter], h_eqchnpair_purity_corrected, niter);
                unfold_neqchnpair[ct_iter] = new RooUnfoldBayes(response_neqchnpair[ct_iter], h_neqchnpair_purity_corrected, niter);

                h_npair_unfolded[ct_iter]      = (TH3D*) unfold_npair[ct_iter]->Hreco();
                h_eqchnpair_unfolded[ct_iter]  = (TH3D*) unfold_eqchnpair[ct_iter]->Hreco();
                h_neqchnpair_unfolded[ct_iter] = (TH3D*) unfold_neqchnpair[ct_iter]->Hreco();

                h_npair_unfolded[ct_iter]->Divide(h_efficiency);
                h_eqchnpair_unfolded[ct_iter]->Divide(h_efficiency);
                h_neqchnpair_unfolded[ct_iter]->Divide(h_efficiency);

                apply_jet_weight_to_npairs(h_npair_unfolded[ct_iter]     , h_purity_jet, h_efficiency_jet);
                apply_jet_weight_to_npairs(h_eqchnpair_unfolded[ct_iter] , h_purity_jet, h_efficiency_jet);
                apply_jet_weight_to_npairs(h_neqchnpair_unfolded[ct_iter], h_purity_jet, h_efficiency_jet);

                apply_unfolded_weights(h_npair_unfolded[ct_iter]     , h_eec_unfolded_2d);
                apply_unfolded_weights(h_eqchnpair_unfolded[ct_iter] , h_eqcheec_unfolded_2d);
                apply_unfolded_weights(h_neqchnpair_unfolded[ct_iter], h_neqcheec_unfolded_2d);
                
                h_npair_unfolded_2d[ct_iter] = (TH2D*) h_npair_unfolded[ct_iter]->Project3D("yx");

                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        double avge_pt2_jet = (jet_pt_binning[bin + 1] + jet_pt_binning[bin])/2.;
                        int nominal_jet_pt_bin = bin + 3;

                        // Pseudodata operations
                        hcorr_eec[bin][ct_iter]      = new TH1F(Form("hcorr_eec%i%i",bin, ct_iter)     , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                        hcorr_eqcheec[bin][ct_iter]  = new TH1F(Form("hcorr_eqcheec%i%i",bin, ct_iter) , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                        hcorr_neqcheec[bin][ct_iter] = new TH1F(Form("hcorr_neqcheec%i%i",bin, ct_iter), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                        hcorr_tau[bin][ct_iter]      = new TH1F(Form("hcorr_tau%i%i",bin, ct_iter)     , "", nbin_rl_nominal,tau_binning[bin]);
                        
                        project_nominal_phase_space(h_eec_unfolded_2d      , hcorr_eec[bin][ct_iter]     , nominal_jet_pt_bin);
                        project_nominal_phase_space(h_eqcheec_unfolded_2d  , hcorr_eqcheec[bin][ct_iter] , nominal_jet_pt_bin);
                        project_nominal_phase_space(h_neqcheec_unfolded_2d , hcorr_neqcheec[bin][ct_iter], nominal_jet_pt_bin);
                        
                        hcorr_eqcheec[bin][ct_iter]->Divide(hcorr_eec[bin][ct_iter]);
                        hcorr_neqcheec[bin][ct_iter]->Divide(hcorr_eec[bin][ct_iter]);

                        hcorr_eec[bin][ct_iter]->Scale(1./h_njet_unfolded->GetBinContent(bin + 3),"width");

                        get_tau_from_uoflow_eec(hcorr_eec[bin][ct_iter], hcorr_tau[bin][ct_iter], avge_pt2_jet);
                        
                        fout->cd();
                        hcorr_eec[bin][ct_iter]->Write();
                        hcorr_tau[bin][ct_iter]->Write();
                        hcorr_eqcheec[bin][ct_iter]->Write();
                        hcorr_neqcheec[bin][ct_iter]->Write();
                        gROOT->cd();

                        hcorr_eec[bin][ct_iter]->Divide(hcorr_eec_truth[bin]);
                        hcorr_eec[bin][ct_iter]->SetTitle("EEC(pseudodata/truth)");
                        hcorr_tau[bin][ct_iter]->Divide(hcorr_tau_truth[bin]);
                        hcorr_tau[bin][ct_iter]->SetTitle("Npair(pseudodata/truth)");
                        hcorr_eqcheec[bin][ct_iter]->Divide(hcorr_eqcheec_truth[bin]);
                        hcorr_eqcheec[bin][ct_iter]->SetTitle("Eq. Charged EEC(pseudodata/truth)");
                        hcorr_neqcheec[bin][ct_iter]->Divide(hcorr_neqcheec_truth[bin]);
                        hcorr_neqcheec[bin][ct_iter]->SetTitle("Op. Charged EEC(pseudodata/truth)");
                        
                        hcorr_eec[bin][ct_iter]->Add(h_unity_uounderflow, -1);
                        hcorr_tau[bin][ct_iter]->Add(h_unity_taubinning[bin], -1);
                        hcorr_eqcheec[bin][ct_iter]->Add(h_unity_uounderflow, -1);
                        hcorr_neqcheec[bin][ct_iter]->Add(h_unity_uounderflow, -1);

                        hcorr_eec[bin][ct_iter]->Multiply(hcorr_eec[bin][ct_iter]);
                        hcorr_tau[bin][ct_iter]->Multiply(hcorr_tau[bin][ct_iter]);
                        hcorr_eqcheec[bin][ct_iter]->Multiply(hcorr_eqcheec[bin][ct_iter]);
                        hcorr_neqcheec[bin][ct_iter]->Multiply(hcorr_neqcheec[bin][ct_iter]);

                        hcorr_eec_total[bin]->Add(hcorr_eec[bin][ct_iter]);  
                        hcorr_tau_total[bin]->Add(hcorr_tau[bin][ct_iter]); 
                        hcorr_eqcheec_total[bin]->Add(hcorr_eqcheec[bin][ct_iter]); 
                        hcorr_neqcheec_total[bin]->Add(hcorr_neqcheec[bin][ct_iter]);

                        
                        fout->cd();
                        hcorr_eec[bin][ct_iter]->Write(Form("squared_relerror_eec%i%i",bin,ct_iter));
                        hcorr_tau[bin][ct_iter]->Write(Form("squared_relerror_npair%i%i",bin,ct_iter));
                        hcorr_eqcheec[bin][ct_iter]->Write(Form("squared_relerror_eqcheec%i%i",bin,ct_iter));
                        hcorr_neqcheec[bin][ct_iter]->Write(Form("squared_relerror_neqcheec%i%i",bin,ct_iter));
                        gROOT->cd();       
                }

                // Recycling
                hmeas_npair->Reset();
                htrue_npair->Reset();
                hmeas_eqchnpair->Reset();
                htrue_eqchnpair->Reset();
                hmeas_neqchnpair->Reset();
                htrue_neqchnpair->Reset();
                
                h_npair_purity_corrected->Reset();
                h_eqchnpair_purity_corrected->Reset();
                h_neqchnpair_purity_corrected->Reset();

                h_eec_unfolded_2d->Reset();
                h_eqcheec_unfolded_2d->Reset();
                h_neqcheec_unfolded_2d->Reset();
        }

        // Estimate the average difference
        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                square_root_bins(hcorr_eec_total[i]);
                square_root_bins(hcorr_tau_total[i]);
                square_root_bins(hcorr_eqcheec_total[i]);
                square_root_bins(hcorr_neqcheec_total[i]);

                hcorr_eec_total[i]->Scale(1./std::sqrt(niter_ct));
                hcorr_tau_total[i]->Scale(1./std::sqrt(niter_ct));
                hcorr_eqcheec_total[i]->Scale(1./std::sqrt(niter_ct));
                hcorr_neqcheec_total[i]->Scale(1./std::sqrt(niter_ct));

                set_histogram_style(hcorr_eec_total[i]     , corr_marker_color_jet_pt[i], std_line_width-1, corr_marker_style_jet_pt[i], std_marker_size);
                set_histogram_style(hcorr_eqcheec_total[i] , corr_marker_color_jet_pt[i], std_line_width-1, std_marker_style_jet_pt[i] , std_marker_size);
                set_histogram_style(hcorr_neqcheec_total[i], corr_marker_color_jet_pt[i], std_line_width-1, corr_marker_style_jet_pt[i], std_marker_size);

                fout->cd();
                hcorr_eec_total[i]->Write(Form("relerror_eec%i",i));
                hcorr_tau_total[i]->Write(Form("relerror_tau%i",i));
                hcorr_eqcheec_total[i]->Write(Form("relerror_eqcheec%i",i));
                hcorr_neqcheec_total[i]->Write(Form("relerror_neqcheec%i",i));
                gROOT->cd();
        }

        // Print all the relevant ratios
        TCanvas* c = new TCanvas("c","",1800,600);
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

                s[bin]->Add(hcorr_eec_total[bin],"PE1");
                
                s[bin]->Draw("NOSTACK");
                
                s[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};Corr. Pseudodata / Truth",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s[bin]->SetMaximum(0.5);
                s[bin]->SetMinimum(-0.5);
                s[bin]->GetXaxis()->SetRangeUser(unfolding_rl_nominal_binning[1],unfolding_rl_nominal_binning[nbin_rl_nominal_unfolding-1]);

                l[bin]->AddEntry(hcorr_eec_total[bin], "EEC", "p");
                
                gPad->SetLogx(1);
                
                l[bin]->Draw("SAME");    
                line->Draw("SAME");
        }

        c->Print(Form("./plots/closure-test-stat-eec-nsamples%i.pdf",niter_ct));

        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                c->cd(bin+1);
                s[bin] = new THStack();
                l[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,gPad->GetBottomMargin()+0.01,1-gPad->GetRightMargin()-0.01,gPad->GetBottomMargin()+0.12);

                s[bin]->Add(hcorr_eqcheec_total[bin],"PE1");
                s[bin]->Add(hcorr_neqcheec_total[bin],"PE1");
                
                s[bin]->Draw("NOSTACK");
                
                s[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};Corr. Pseudodata / Truth",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s[bin]->SetMaximum(0.5);
                s[bin]->SetMinimum(-0.5);
                s[bin]->GetXaxis()->SetRangeUser(unfolding_rl_nominal_binning[1],unfolding_rl_nominal_binning[nbin_rl_nominal_unfolding-1]);

                l[bin]->AddEntry(hcorr_eqcheec_total[bin], "Eq. Ch. EEC", "p");
                l[bin]->AddEntry(hcorr_neqcheec_total[bin], "Op. Ch. EEC", "p");
                
                gPad->SetLogx(1);
                
                l[bin]->Draw("SAME");    
                line->Draw("SAME");
        }

        c->Print(Form("./plots/closure-test-stat-chargedeec-nsamples%i.pdf",niter_ct));
}
