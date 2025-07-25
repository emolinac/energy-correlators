#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"
#include "TRandom3.h"

void macro_print_fullcorre2c_paircorr_2dunf_ct_niter(int niter = nominal_niter, int ct_niter = 10, bool do_print = true, bool compare_to_truth = true)
{
        std::string systematic = available_systematics[1]; // choose CT systematic
        
        if (systematic != "ct") {
                std::cout<<"Make sure to set the CT systematic"<<std::endl; 
                return;
        }

        TFile* fout_dev = new TFile((output_folder + devfromnom_namef[systematic]).c_str(),"RECREATE");
        TFile* fout     = new TFile((output_folder + namef_histos_paircorr_e2c_ct).c_str(),"RECREATE");
        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_ntuple_e2c_paircorr_ct).c_str()); 
        if (fcorr->IsZombie()) 
                return;

        TFile* fnominal = new TFile((output_folder + namef_histos_paircorr_e2c).c_str());
        if (fnominal->IsZombie()) 
                return;
        
        TNtuple* ntuple_data   = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet    = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        TNtuple* ntuple_mc     = (TNtuple*) fcorr->Get((name_ntuple_mc).c_str());
        TNtuple* ntuple_mc_jet = (TNtuple*) fcorr->Get((name_ntuple_mc_jet).c_str());
        
        // Set the branches of data
        float R_L, jet_pt, weight_pt, efficiency, purity, efficiency_relerror, purity_relerror;
        set_data_ntuple_branches(ntuple_data, &R_L, &jet_pt, &weight_pt, &efficiency, &purity, &efficiency_relerror, &purity_relerror);
        
        // Unfold the purity corrected data
        TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections_ct).c_str());
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        
        TH2D* hpurcorr = new TH2D("hpurcorr","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hmeas    = new TH2D("hmeas"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* htrue    = new TH2D("htrue"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");
        
        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                ntuple->GetEntry(evt);

                if (abs(R_L_truth - R_L_reco) < rl_resolution) 
                        response->Fill(R_L_reco, jet_pt_reco, R_L_truth, jet_pt_truth);
        }

        TH2D* hunfolded_ratio     = new TH2D("hunfolded_ratio"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hpuritycorrected    = new TH2D("hpuritycorrected" ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hpuritycorrected2   = new TH2D("hpuritycorrected2","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        ntuple_data->Project("hpuritycorrected" , "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
        ntuple_data->Project("hpuritycorrected2", "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
        
        RooUnfoldBayes unfold(response, hpuritycorrected, niter);
        TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
        hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected2,1,1);

        TH1F* hcorr_jet[nbin_jet_pt];
        TH1F* hcorr_jet_centroid[nbin_jet_pt];
        TH1F* hcorr_e2c_nonorm[nbin_jet_pt]; 
        TH1F* hcorr_e2c[nbin_jet_pt]; 
        TH1F* hcorr_e2c_nounf[nbin_jet_pt]; 
        TH1F* hcorr_tau[nbin_jet_pt]; 
        TH1F* hcorr_tau_nounf[nbin_jet_pt]; 

        TH1F* hcorr_ratio_e2c[nbin_jet_pt][ct_niter]; 
        TH1F* hcorr_ratio_e2c_total[nbin_jet_pt]; 
        TH1F* hcorr_ratio_tau[nbin_jet_pt][ct_niter]; 
        TH1F* hcorr_ratio_tau_total[nbin_jet_pt]; 
        TH1F* hnominal[nbin_jet_pt];

        TH1F* htruth_jet[nbin_jet_pt];
        TH1F* htruth[nbin_jet_pt]; 
        TH1F* htruth_tau[nbin_jet_pt];     

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        // Adding content with errors
        TLatex latex;
        latex.SetTextAlign(22); // center alignment
        latex.SetTextSize(rl_resolution);
        latex.SetTextColor(kBlack);

        gStyle->SetPaintTextFormat("4.2f");
        hunfolded_ratio->Draw("col");

        for (int i = 2; i < hunfolded_ratio->GetNbinsX(); ++i) {
                for (int j = 2; j < hunfolded_ratio->GetNbinsY(); ++j) {
                        double x = hunfolded_ratio->GetXaxis()->GetBinCenter(i);
                        double y = hunfolded_ratio->GetYaxis()->GetBinCenter(j);
                        double content = hunfolded_ratio->GetBinContent(i, j);
                        double error = hunfolded_ratio->GetBinError(i, j);
                        // Draw content and error in the format "content ± error"
                        latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                }
        }

        hunfolded_ratio->SetTitle("Purity Corrected Unfolded/Purity Corrected;R_{L};p_{T,jet} (GeV)");
        hunfolded_ratio->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
        hunfolded_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        if (do_print) c->Print(Form("./plots/unfolded2d_niter%i_ratio_ctniter%i.pdf",niter,ct_niter));

        THStack* s_data     = new THStack();
        TLegend* l_data     = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);
        TLegend* l_data_tau = new TLegend(gPad->GetLeftMargin()+0.01,1-gPad->GetTopMargin()-0.01,gPad->GetLeftMargin()+0.21,1-gPad->GetTopMargin()-0.21);

        // Construct the histograms you will use in order to avoid memory leaks
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)         ,"",1  ,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                hcorr_jet_centroid[bin] = new TH1F(Form("hcorr_jet_centroid%i" ,bin),"",200,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                hcorr_e2c_nonorm[bin]   = new TH1F(Form("hcorr_e2c_nonorm%i",bin)   ,"",nbin_rl_nominal,rl_nominal_binning );
                hcorr_e2c[bin]          = new TH1F(Form("hcorr_e2c%i",bin)          ,"",nbin_rl_nominal,rl_nominal_binning );
                hcorr_e2c_nounf[bin]    = new TH1F(Form("hcorr_e2c_nounf%i",bin)    ,"",nbin_rl_nominal,rl_nominal_binning );
                hcorr_tau[bin]          = new TH1F(Form("hcorr_tau%i",bin)          ,"",nbin_rl_nominal,tau_nominal_binning);
                hcorr_tau_nounf[bin]    = new TH1F(Form("hcorr_tau_nounf%i",bin)    ,"",nbin_rl_nominal,tau_nominal_binning);

                // set_histogram_style(hcorr_e2c[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                // set_histogram_style(hcorr_tau[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        
                ntuple_jet->Project(Form("hcorr_jet%i" ,bin)         , "jet_pt",jet_full_corr[bin]);
                ntuple_jet->Project(Form("hcorr_jet_centroid%i" ,bin), "jet_pt",jet_full_corr[bin]);

                
                hnominal[bin] = (TH1F*) fnominal->Get(Form("hcorr_e2c%i",bin));

                hcorr_ratio_e2c_total[bin] = new TH1F(Form("hcorr_ratio_e2c_total%i",bin),"",nbin_rl_nominal,rl_nominal_binning);
                htruth[bin]                = new TH1F(Form("htruth%i",bin)               ,"",nbin_rl_nominal,rl_nominal_binning);

                htruth_tau[bin]            = new TH1F(Form("htruth_tau%i",bin)           ,"",nbin_rl_nominal,tau_nominal_binning);
                hcorr_ratio_tau_total[bin] = new TH1F(Form("hcorr_ratio_tau_total%i",bin),"",nbin_rl_nominal,tau_nominal_binning);

                htruth_jet[bin]            = new TH1F(Form("htruth_jet%i" ,bin)          ,"",200,jet_pt_binning[bin],jet_pt_binning[bin + 1]);     

                ntuple_mc_jet->Project(Form("htruth_jet%i" ,bin), "jet_pt");
                ntuple_mc->Project(Form("htruth%i",bin),"R_L",e2c_jet_pt_cut[bin]);
                ntuple_mc->Project(Form("htruth_tau%i",bin),Form("R_L*%f",htruth_jet[bin]->GetMean()),e2c_jet_pt_cut[bin]);

                // Normalize the truth distributions to unity for further comparison with corr pseudodata
                htruth[bin]->Scale(1./htruth[bin]->Integral(),"width");
                htruth_tau[bin]->Scale(1./htruth_tau[bin]->Integral(),"width");
        }

        // Create a data histograma to know how much you have to vary its "content"
        TH2D* hdataunf_ref = new TH2D("hdataunf_ref","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hdatashift   = new TH2D("hdatashift"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);    
        ntuple_data->Project("hdataunf_ref", "jet_pt:R_L");

        TRandom3* rndm = new TRandom3();
        for (int ct_iter = 0 ; ct_iter < ct_niter ; ct_iter++) {
                std::cout<<"Iteration "<<ct_iter<<std::endl;

                set_shift_histo(hdataunf_ref,hdatashift,rndm);

                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hcorr_ratio_e2c[bin][ct_iter] = new TH1F(Form("hcorr_ratio_e2c%i%i",bin,ct_iter),"",nbin_rl_nominal,rl_nominal_binning);
                        hcorr_ratio_tau[bin][ct_iter] = new TH1F(Form("hcorr_ratio_tau%i%i",bin,ct_iter),"",nbin_rl_nominal,tau_nominal_binning);
                        
                        double jet_pt_centroid = hcorr_jet_centroid[bin]->GetMean();

                        for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++) {
                                ntuple_data->GetEntry(entry);
                
                                if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) 
                                        continue;
                                
                                double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt));
                                if (unfolding_weight <= 0) 
                                        unfolding_weight = 1;
                
                                double content_shift = hdatashift->GetBinContent(hdatashift->FindBin(R_L,jet_pt));
                
                                hcorr_e2c_nonorm[bin]->Fill(R_L,purity*unfolding_weight*weight_pt*content_shift/efficiency);
                                hcorr_e2c[bin]->Fill(R_L,purity*unfolding_weight*weight_pt*content_shift/efficiency);
                                hcorr_e2c_nounf[bin]->Fill(R_L,purity*weight_pt*content_shift/efficiency);
                                hcorr_tau[bin]->Fill(R_L*jet_pt_centroid,purity*unfolding_weight*weight_pt*content_shift/efficiency);
                                hcorr_tau_nounf[bin]->Fill(R_L*jet_pt_centroid,purity*weight_pt*content_shift/efficiency);
                        }
                
                        // Normalize the distributions to unity
                        // Log binning
                        hcorr_e2c[bin]->Scale(1./hcorr_e2c[bin]->Integral(),"width");
                        hcorr_e2c_nounf[bin]->Scale(1./hcorr_e2c_nounf[bin]->Integral(),"width");

                        hcorr_tau[bin]->Scale(1./hcorr_tau[bin]->Integral(),"width");
                        hcorr_tau_nounf[bin]->Scale(1./hcorr_tau_nounf[bin]->Integral(),"width");
                        
                        // Get the delta corresponding to the ith iteration
                        hcorr_ratio_e2c[bin][ct_iter]->Divide(hcorr_e2c[bin],htruth[bin],1,1);
                        hcorr_ratio_e2c_total[bin]->Add(hcorr_ratio_e2c[bin][ct_iter],1);

                        hcorr_ratio_tau[bin][ct_iter]->Divide(hcorr_tau[bin],htruth_tau[bin],1,1);
                        hcorr_ratio_tau_total[bin]->Add(hcorr_ratio_tau[bin][ct_iter],1);

                        // Reset histograms to use again
                        hcorr_e2c[bin]->Reset();
                        hcorr_e2c_nounf[bin]->Reset();

                        hcorr_tau[bin]->Reset();
                        hcorr_tau_nounf[bin]->Reset();
                }
        }
        
        // Scale the relevant histograms
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_ratio_e2c_total[bin]->Scale(1./ct_niter);
                hcorr_ratio_tau_total[bin]->Scale(1./ct_niter);
                
                fout->cd();
                hcorr_ratio_e2c_total[bin]->Write();
                hcorr_ratio_tau_total[bin]->Write();
                gROOT->cd();
        }

        if (compare_to_truth) {
                TH2F* hct_ratio = new TH2F("hct_ratio","",nbin_rl_nominal,rl_nominal_binning,nbin_jet_pt,jet_pt_binning);
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {   
                        // Normalize both to unity such that we can compare the shapes
                        for (int bin_rl = 1 ; bin_rl <= hcorr_e2c_nonorm[bin]->GetNbinsX() ; bin_rl++) {
                                hct_ratio->SetBinContent(bin_rl, bin + 1, hcorr_ratio_e2c_total[bin]->GetBinContent(bin_rl));
                                hct_ratio->SetBinError(bin_rl, bin + 1, hcorr_ratio_e2c_total[bin]->GetBinError(bin_rl));
                        } 

                        // // For clarity in drawing
                        // set_histo_null_errors(hcorr_ratio_e2c_total[bin]);
                        // set_histo_null_errors(hcorr_ratio_tau_total[bin]);

                        fout_dev->cd();
                        hcorr_ratio_e2c_total[bin]->Write(Form("h_deviations%i",bin));
                        hcorr_ratio_tau_total[bin]->Write(Form("h_deviations_tau%i",bin));
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
                if (do_print) c->Print(Form("./plots/closuretest_niter%i_ratio_ctniter%i.pdf",niter,ct_niter));
        }
}
