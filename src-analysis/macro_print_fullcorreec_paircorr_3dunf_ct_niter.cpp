#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"
#include "../include/TZJets2016Data.h"
#include "../include/TZJets2016Data.C"
#include "TRandom3.h"

void macro_print_fullcorreec_paircorr_3dunf_ct_niter(int niter = nominal_niter, int ct_niter = 10, bool do_print = true, bool compare_to_truth = true)
{
        std::string systematic = available_systematics[1]; // choose CT systematic
        
        if (systematic != "ct") {
                std::cout<<"Make sure to set the CT systematic"<<std::endl; 
                return;
        }

        TFile* fout_dev = new TFile((output_folder + devfromnom_namef[systematic]).c_str(),"RECREATE");
        TFile* fout     = new TFile((output_folder + namef_histos_paircorr_eec_ct).c_str(),"RECREATE");
        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_ntuple_eec_paircorr_ct).c_str()); 
        if (fcorr->IsZombie()) 
                return;

        TNtuple* ntuple_pseudodata = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet        = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        TNtuple* ntuple_mc         = (TNtuple*) fcorr->Get((name_ntuple_mc).c_str());
        TNtuple* ntuple_mc_jet     = (TNtuple*) fcorr->Get((name_ntuple_mc_jet).c_str());

        // Set the branches of data
        float R_L, jet_pt, weight_pt, event_weight, efficiency, purity, efficiency_relerror, purity_relerror, eq_charge, h1_pt, h2_pt;
        set_data_ntuple_branches(ntuple_pseudodata, &event_weight, &R_L, &jet_pt, &weight_pt, &efficiency, &purity, &efficiency_relerror, &purity_relerror, &eq_charge);
        ntuple_pseudodata->SetBranchAddress("h1_pt",&h1_pt);
        ntuple_pseudodata->SetBranchAddress("h2_pt",&h2_pt);
        
        // Set the response matrix to unfold down the road
        TFile* f = new TFile((output_folder + namef_ntuple_eec_paircorrections_ct).c_str());
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth, h1_pt_reco, h1_pt_truth, h2_pt_reco, h2_pt_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        ntuple->SetBranchAddress("h1_pt",&h1_pt_reco);
        ntuple->SetBranchAddress("h1_pt_truth",&h1_pt_truth);
        ntuple->SetBranchAddress("h2_pt",&h2_pt_reco);
        ntuple->SetBranchAddress("h2_pt_truth",&h2_pt_truth);
        
        TH3D* hpurcorr = new TH3D("hpurcorr","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_ptprod, ptprod_binning);
        TH3D* hmeas    = new TH3D("hmeas"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_ptprod, ptprod_binning);
        TH3D* htrue    = new TH3D("htrue"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_ptprod, ptprod_binning);
        RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");
        
        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                ntuple->GetEntry(evt);

                if (abs(R_L_truth - R_L_reco) < rl_resolution) 
                        response->Fill(R_L_reco, jet_pt_reco, h1_pt_reco*h2_pt_reco, R_L_truth, jet_pt_truth, h1_pt_truth*h2_pt_truth);
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        // Adding content with errors
        TLatex latex;
        latex.SetTextAlign(22); // center alignment
        latex.SetTextSize(text_size_correction_plots);
        latex.SetTextColor(kBlack);

        gStyle->SetPaintTextFormat("4.2f");
        
        THStack* s_data     = new THStack();
        TLegend* l_data     = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);
        TLegend* l_data_tau = new TLegend(gPad->GetLeftMargin()+0.01,1-gPad->GetTopMargin()-0.01,gPad->GetLeftMargin()+0.21,1-gPad->GetTopMargin()-0.21);

        TH1F* hcorr_jet[nbin_jet_pt];
        TH1F* hcorr_jet_centroid[nbin_jet_pt];
        TH1F* hcorr_npair[nbin_jet_pt]; 
        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_tau[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_total[nbin_jet_pt]; // necessary due to the difference in the type of binning
        
        TH1F* hcorr_ratio_eec[nbin_jet_pt][ct_niter]; 
        TH1F* hcorr_ratio_eec_total[nbin_jet_pt]; 
        TH1F* hcorr_ratio_npair[nbin_jet_pt][ct_niter]; 
        TH1F* hcorr_ratio_npair_total[nbin_jet_pt]; 
        TH1F* hcorr_ratio_tau[nbin_jet_pt][ct_niter]; 
        TH1F* hcorr_ratio_tau_total[nbin_jet_pt]; 
        TH1F* hcorr_ratio_eec_eqcharge[nbin_jet_pt][ct_niter]; 
        TH1F* hcorr_ratio_eec_eqcharge_total[nbin_jet_pt]; 
        TH1F* hcorr_ratio_eec_neqcharge[nbin_jet_pt][ct_niter]; 
        TH1F* hcorr_ratio_eec_neqcharge_total[nbin_jet_pt]; 
        
        TH1F* hnominal[nbin_jet_pt];

        TH1F* htruth_jet[nbin_jet_pt];
        TH1F* htruth_eec[nbin_jet_pt]; 
        TH1F* htruth_tau[nbin_jet_pt]; 
        TH1F* htruth_npair[nbin_jet_pt];     
        TH1F* htruth_eec_eqcharge[nbin_jet_pt]; 
        TH1F* htruth_eec_neqcharge[nbin_jet_pt];     
        TH1F* htruth_eec_total[nbin_jet_pt];

        // Initialize the histograms you will use in order to avoid memory leaks
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)         ,"",1  ,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                hcorr_jet_centroid[bin] = new TH1F(Form("hcorr_jet_centroid%i" ,bin),"",200,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                htruth_jet[bin]         = new TH1F(Form("htruth_jet%i" ,bin)        ,"",200,jet_pt_binning[bin],jet_pt_binning[bin + 1]);     

                ntuple_jet->Project(Form("hcorr_jet%i" ,bin)         , "jet_pt",jet_full_corr[bin]);
                ntuple_jet->Project(Form("hcorr_jet_centroid%i" ,bin), "jet_pt",jet_full_corr[bin]);
                ntuple_mc_jet->Project(Form("htruth_jet%i" ,bin), "jet_pt");
                
                hcorr_npair[bin]         = new TH1F(Form("hcorr_npair%i",bin)        ,"",nbin_rl_nominal,rl_nominal_binning );
                hcorr_eec[bin]           = new TH1F(Form("hcorr_eec%i",bin)          ,"",nbin_rl_nominal,rl_nominal_binning );
                hcorr_tau[bin]           = new TH1F(Form("hcorr_tau%i",bin)          ,"",nbin_rl_nominal,tau_nominal_binning);
                hcorr_eec_eqcharge[bin]  = new TH1F(Form("hcorr_eec_eqcharge%i",bin) ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                hcorr_eec_neqcharge[bin] = new TH1F(Form("hcorr_eec_neqcharge%i",bin),"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                hcorr_eec_total[bin]     = new TH1F(Form("hcorr_eec_total%i",bin)    ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                
                htruth_eec[bin]           = new TH1F(Form("htruth_eec%i",bin)           ,"",nbin_rl_nominal,rl_nominal_binning);
                htruth_tau[bin]           = new TH1F(Form("htruth_tau%i",bin)           ,"",nbin_rl_nominal,tau_nominal_binning);
                htruth_npair[bin]         = new TH1F(Form("htruth_npair%i",bin)         ,"",nbin_rl_nominal,rl_nominal_binning);
                htruth_eec_eqcharge[bin]  = new TH1F(Form("htruth_eec_eqcharge%i",bin)  ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                htruth_eec_neqcharge[bin] = new TH1F(Form("htruth_eec_neqcharge%i",bin) ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                htruth_eec_total[bin]     = new TH1F(Form("htruth_eec_total%i",bin)     ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                
                hcorr_ratio_eec_total[bin]           = new TH1F(Form("hcorr_ratio_eec_total%i",bin)  ,"",nbin_rl_nominal,rl_nominal_binning);
                hcorr_ratio_npair_total[bin]         = new TH1F(Form("hcorr_ratio_npair_total%i",bin),"",nbin_rl_nominal,rl_nominal_binning);
                hcorr_ratio_tau_total[bin]           = new TH1F(Form("hcorr_ratio_tau_total%i",bin)  ,"",nbin_rl_nominal,tau_nominal_binning);
                hcorr_ratio_eec_eqcharge_total[bin]  = new TH1F(Form("hcorr_ratio_eec_eqcharge_total%i",bin) ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                hcorr_ratio_eec_neqcharge_total[bin] = new TH1F(Form("hcorr_ratio_eec_neqcharge_total%i",bin),"",nbin_chargedeec_nominal,rl_chargedeec_binning);

                ntuple_mc->Project(Form("htruth_eec%i",bin),"R_L",eec_jet_pt_cut[bin]);
                ntuple_mc->Project(Form("htruth_npair%i",bin),"R_L",pair_jet_pt_cut[bin]);
                ntuple_mc->Project(Form("htruth_tau%i",bin),Form("R_L*%f",htruth_jet[bin]->GetMean()),eec_jet_pt_cut[bin]);
                ntuple_mc->Project(Form("htruth_eec_eqcharge%i",bin),"R_L",eec_eqcharge_jet_pt_cut[bin]);
                ntuple_mc->Project(Form("htruth_eec_neqcharge%i",bin),"R_L",eec_neqcharge_jet_pt_cut[bin]);
                ntuple_mc->Project(Form("htruth_eec_total%i",bin),"R_L",eec_jet_pt_cut[bin]);

                htruth_npair[bin]->Scale(1./htruth_jet[bin]->Integral());
                htruth_eec[bin]->Scale(1./htruth_jet[bin]->Integral());
                htruth_tau[bin]->Scale(1./htruth_jet[bin]->Integral());

                htruth_eec_total[bin]->Add(htruth_eec_eqcharge[bin],htruth_eec_neqcharge[bin],1,1);
                htruth_eec_eqcharge[bin]->Divide(htruth_eec_total[bin]);
                htruth_eec_neqcharge[bin]->Divide(htruth_eec_total[bin]);
        }

        TFile* fdata = new TFile((output_folder + namef_ntuple_eec_paircorr).c_str());

        TNtuple* ntuple_data = (TNtuple*) fdata->Get(name_ntuple_data.c_str());

        TRandom3* rndm = new TRandom3(0);
        for (int ct_iter = 0 ; ct_iter < ct_niter ; ct_iter++) {
                std::cout<<"Iteration "<<ct_iter<<std::endl;

                TH3D* hdataunf_ref = new TH3D("hdataunf_ref","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_ptprod, ptprod_binning);
                TH3D* hdatashift   = new TH3D("hdatashift"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_ptprod, ptprod_binning);    
                
                TH3D* hunfolded_ratio      = new TH3D("hunfolded_ratio"     ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_ptprod, ptprod_binning);
                TH3D* hpuritycorrected     = new TH3D("hpuritycorrected"    ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_ptprod, ptprod_binning);
                TH3D* hpuritycorrected_ref = new TH3D("hpuritycorrected_ref","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_ptprod, ptprod_binning);
                
                ntuple_data->Project("hdataunf_ref", "h1_pt*h2_pt:jet_pt:R_L");
                set_shift_histo(hdataunf_ref,hdatashift,rndm);

                ntuple_pseudodata->Project("hpuritycorrected"    ,"h1_pt*h2_pt:jet_pt:R_L","purity");
                ntuple_pseudodata->Project("hpuritycorrected_ref","h1_pt*h2_pt:jet_pt:R_L","purity");

                hpuritycorrected->Multiply(hdatashift);
                hpuritycorrected_ref->Multiply(hdatashift);
                
                RooUnfoldBayes* unfold = new RooUnfoldBayes(response, hpuritycorrected, niter);

                TH3D* hunfolded_bayes = (TH3D*) unfold->Hunfold();
                
                hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected_ref,1,1);
                
                TH2D* hunfolded_ratio_jet_pt_rl = (TH2D*) hunfolded_ratio->Project3D("yx");
                
                hunfolded_ratio_jet_pt_rl->Draw("col");

                for (int i = 2; i < hunfolded_ratio_jet_pt_rl->GetNbinsX(); ++i) {
                        for (int j = 2; j < hunfolded_ratio_jet_pt_rl->GetNbinsY(); ++j) {
                                double x = hunfolded_ratio_jet_pt_rl->GetXaxis()->GetBinCenter(i);
                                double y = hunfolded_ratio_jet_pt_rl->GetYaxis()->GetBinCenter(j);
                                double content = hunfolded_ratio_jet_pt_rl->GetBinContent(i, j);
                                double error = hunfolded_ratio_jet_pt_rl->GetBinError(i, j);
                                // Draw content and error in the format "content ± error"
                                latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                        }
                }
                
                hunfolded_ratio_jet_pt_rl->SetTitle("Unf3D: Purity Corrected Unfolded/Purity Corrected;R_{L};p_{T,jet} (GeV)");
                hunfolded_ratio_jet_pt_rl->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                hunfolded_ratio_jet_pt_rl->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                
                if (do_print) 
                        c->Print(Form("./plots/unfolded3d_3dunf-niter%i_ratio_ct.pdf",niter));

                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hcorr_ratio_eec[bin][ct_iter]   = new TH1F(Form("hcorr_ratio_eec%i%i",bin,ct_iter),"",nbin_rl_nominal,rl_nominal_binning);
                        hcorr_ratio_tau[bin][ct_iter]   = new TH1F(Form("hcorr_ratio_tau%i%i",bin,ct_iter),"",nbin_rl_nominal,tau_nominal_binning);
                        hcorr_ratio_npair[bin][ct_iter] = new TH1F(Form("hcorr_ratio_npair%i%i",bin,ct_iter),"",nbin_rl_nominal,rl_nominal_binning);
                        hcorr_ratio_eec_eqcharge[bin][ct_iter] = new TH1F(Form("hcorr_ratio_eec_eqcharge%i%i",bin,ct_iter),"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                        hcorr_ratio_eec_neqcharge[bin][ct_iter] = new TH1F(Form("hcorr_ratio_eec_neqcharge%i%i",bin,ct_iter),"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                        
                        double jet_pt_centroid = hcorr_jet_centroid[bin]->GetMean();

                        for (int entry = 0 ; entry < ntuple_pseudodata->GetEntries() ; entry++) {
                                ntuple_pseudodata->GetEntry(entry);
                
                                if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) 
                                        continue;
                                
                                double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt,h1_pt*h2_pt));
                                
                                if (unfolding_weight <= 0) 
                                        unfolding_weight = 1;
                
                                double content_shift = hdatashift->GetBinContent(hdatashift->FindBin(R_L,jet_pt,h1_pt*h2_pt));
                
                                hcorr_npair[bin]->Fill(R_L,event_weight*purity*unfolding_weight*content_shift/efficiency);
                                hcorr_eec[bin]->Fill(R_L,event_weight*purity*unfolding_weight*weight_pt*content_shift/efficiency);
                                hcorr_tau[bin]->Fill(R_L*jet_pt_centroid,event_weight*purity*unfolding_weight*weight_pt*content_shift/efficiency);

                                if (eq_charge > 0)
                                        hcorr_eec_eqcharge[bin]->Fill(R_L,event_weight*purity*unfolding_weight*weight_pt/efficiency);
                                else if (eq_charge < 0)
                                        hcorr_eec_neqcharge[bin]->Fill(R_L,event_weight*purity*unfolding_weight*weight_pt/efficiency);
                        }

                        hcorr_eec_total[bin]->Add(hcorr_eec_eqcharge[bin],hcorr_eec_neqcharge[bin],1,1);
                        hcorr_eec_eqcharge[bin]->Divide(hcorr_eec_total[bin]);
                        hcorr_eec_neqcharge[bin]->Divide(hcorr_eec_total[bin]);

                        hcorr_npair[bin]->Scale(1./hcorr_jet[bin]->Integral());
                        hcorr_eec[bin]->Scale(1./hcorr_jet[bin]->Integral());
                        hcorr_tau[bin]->Scale(1./hcorr_jet[bin]->Integral());
                
                        // Get the delta corresponding to the ith iteration
                        hcorr_ratio_eec[bin][ct_iter]->Divide(hcorr_eec[bin],htruth_eec[bin],1,1);
                        hcorr_ratio_eec_total[bin]->Add(hcorr_ratio_eec[bin][ct_iter],1);

                        hcorr_ratio_npair[bin][ct_iter]->Divide(hcorr_npair[bin],htruth_npair[bin],1,1);
                        hcorr_ratio_npair_total[bin]->Add(hcorr_ratio_npair[bin][ct_iter],1);

                        hcorr_ratio_tau[bin][ct_iter]->Divide(hcorr_tau[bin],htruth_tau[bin],1,1);
                        hcorr_ratio_tau_total[bin]->Add(hcorr_ratio_tau[bin][ct_iter],1);

                        hcorr_ratio_eec_eqcharge[bin][ct_iter]->Divide(hcorr_eec_eqcharge[bin],htruth_eec_eqcharge[bin],1,1);
                        hcorr_ratio_eec_eqcharge_total[bin]->Add(hcorr_ratio_eec_eqcharge[bin][ct_iter],1);

                        hcorr_ratio_eec_neqcharge[bin][ct_iter]->Divide(hcorr_eec_neqcharge[bin],htruth_eec_neqcharge[bin],1,1);
                        hcorr_ratio_eec_neqcharge_total[bin]->Add(hcorr_ratio_eec_neqcharge[bin][ct_iter],1);

                        // Reset histograms to use again
                        hcorr_eec[bin]->Reset();
                        hcorr_npair[bin]->Reset();
                        hcorr_tau[bin]->Reset();
                        hcorr_eec_eqcharge[bin]->Reset();
                        hcorr_eec_neqcharge[bin]->Reset();
                        hcorr_eec_total[bin]->Reset();
                }

                hdataunf_ref->Delete();
                hdatashift->Delete();
                hunfolded_ratio->Delete();
                hpuritycorrected->Delete();
                hpuritycorrected_ref->Delete();
                hunfolded_bayes->Delete();
                unfold->Delete();
        }
        
        // Scale the relevant histograms
        if (ct_niter > 1) {
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hcorr_ratio_eec_total[bin]->Scale(1./ct_niter);
                        hcorr_ratio_npair_total[bin]->Scale(1./ct_niter);
                        hcorr_ratio_tau_total[bin]->Scale(1./ct_niter);
                        hcorr_ratio_eec_eqcharge_total[bin]->Scale(1./ct_niter);
                        hcorr_ratio_eec_neqcharge_total[bin]->Scale(1./ct_niter);
                        
                        fout->cd();
                        hcorr_ratio_eec_total[bin]->Write();
                        hcorr_ratio_npair_total[bin]->Write();
                        hcorr_ratio_tau_total[bin]->Write();
                        hcorr_ratio_eec_eqcharge_total[bin]->Write();
                        hcorr_ratio_eec_neqcharge_total[bin]->Write();
                        gROOT->cd();
                }
        }
        
        if (compare_to_truth) {
                TH2F* hct_ratio        = new TH2F("hct_ratio"       ,"",nbin_rl_nominal,rl_nominal_binning,nbin_jet_pt,jet_pt_binning);
                TH2F* hct_npairs_ratio = new TH2F("hct_npairs_ratio","",nbin_rl_nominal,rl_nominal_binning,nbin_jet_pt,jet_pt_binning);
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {   
                        for (int bin_rl = 1 ; bin_rl <= hcorr_npair[bin]->GetNbinsX() ; bin_rl++) {
                                hct_ratio->SetBinContent(bin_rl, bin + 1, hcorr_ratio_eec_total[bin]->GetBinContent(bin_rl));
                                hct_ratio->SetBinError(bin_rl, bin + 1, hcorr_ratio_eec_total[bin]->GetBinError(bin_rl));

                                hct_npairs_ratio->SetBinContent(bin_rl, bin + 1, hcorr_ratio_npair_total[bin]->GetBinContent(bin_rl));
                                hct_npairs_ratio->SetBinError(bin_rl, bin + 1, hcorr_ratio_npair_total[bin]->GetBinError(bin_rl));
                        } 

                        fout_dev->cd();
                        hcorr_ratio_eec_total[bin]->Write(Form("h_deviations_eec%i",bin));
                        hcorr_ratio_npair_total[bin]->Write(Form("h_deviations_npair%i",bin));
                        hcorr_ratio_tau_total[bin]->Write(Form("h_deviations_tau%i",bin));
                        hcorr_ratio_eec_eqcharge_total[bin]->Write(Form("h_deviations_eec_eqcharge%i",bin));
                        hcorr_ratio_eec_neqcharge_total[bin]->Write(Form("h_deviations_eec_neqcharge%i",bin));
                        gROOT->cd();
                }
                
                hct_ratio->Draw("col");
                
                for (int i = 1; i <= hct_ratio->GetNbinsX(); ++i) {
                        for (int j = 1; j <= hct_ratio->GetNbinsY(); ++j) {
                                double x = hct_ratio->GetXaxis()->GetBinCenter(i);
                                double y = hct_ratio->GetYaxis()->GetBinCenter(j);
                                double content = hct_ratio->GetBinContent(i, j);
                                double error = hct_ratio->GetBinError(i, j);
                                // Draw content and error in the format "content ± error"
                                latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                        }
                }

                hct_ratio->SetTitle("Closure Test: Norm. Corr. Pseudodata / Norm. Truth ;R_{L};p_{T,jet} (GeV)");
                hct_ratio->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                hct_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                if (do_print) 
                        c->Print(Form("./plots/closuretest_3dunf-niter%i_ratio_ctniter%i.pdf",niter,ct_niter));

                hct_npairs_ratio->Draw("col");
                
                for (int i = 1; i <= hct_npairs_ratio->GetNbinsX(); ++i) {
                        for (int j = 1; j <= hct_npairs_ratio->GetNbinsY(); ++j) {
                                double x = hct_npairs_ratio->GetXaxis()->GetBinCenter(i);
                                double y = hct_npairs_ratio->GetYaxis()->GetBinCenter(j);
                                double content = hct_npairs_ratio->GetBinContent(i, j);
                                double error = hct_npairs_ratio->GetBinError(i, j);
                                // Draw content and error in the format "content ± error"
                                latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                        }
                }

                hct_npairs_ratio->SetTitle("Closure Test Npairs: Norm. Corr. Pseudodata / Norm. Truth ;R_{L};p_{T,jet} (GeV)");
                hct_npairs_ratio->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                hct_npairs_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                // hct_npairs_ratio->Smooth();
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                if (do_print) 
                        c->Print(Form("./plots/closuretest_3dunf-niter%i_ratio_ctniter%i_npairs.pdf",niter,ct_niter));
        }
}
