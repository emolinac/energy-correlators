#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorre2c_paircorr_2dunf_inclusive(int niter = 4, bool do_print = true, bool compare_to_nominal = false)
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + namef_histos_paircorr_e2c_logbin_inclusive).c_str(),"RECREATE");
        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_ntuple_e2c_paircorr_inclusive).c_str()); 
        if (fcorr->IsZombie()) 
                return;
        
        TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        // Set the branches of data
        float R_L, jet_pt, weight_pt, efficiency, purity, efficiency_relerror, purity_relerror;
        set_data_ntuple_branches(ntuple_data, &R_L, &jet_pt, &weight_pt, &efficiency, &purity, &efficiency_relerror, &purity_relerror);
        
        // UNFOLDING FIRST
        TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections_inclusive).c_str());
        
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        
        // Create histograms with different types of binning
        TH2D* hpurcorr = new TH2D("hpurcorr","",Nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        TH2D* hmeas    = new TH2D("hmeas"   ,"",Nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        TH2D* htrue    = new TH2D("htrue"   ,"",Nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");
        
        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                ntuple->GetEntry(evt);

                if (abs(R_L_truth - R_L_reco) < rl_resolution) 
                        response->Fill(R_L_reco, jet_pt_reco, R_L_truth, jet_pt_truth);
        }

        // Fill the purity corrected distributions
        TH2D* hunfolded_ratio     = new TH2D("hunfolded_ratio"  ,"",Nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        TH2D* hpuritycorrected    = new TH2D("hpuritycorrected" ,"",Nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        TH2D* hpuritycorrected2   = new TH2D("hpuritycorrected2","",Nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        
        ntuple_data->Project("hpuritycorrected" , "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
        ntuple_data->Project("hpuritycorrected2", "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
        
        // Unfold the purity corrected pairs
        RooUnfoldBayes unfold(response, hpuritycorrected, niter);
        TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
        hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected2,1,1);

        TH1F* hcorr_jet[Nbin_jet_pt];
        TH1F* hcorr_jet_centroid[Nbin_jet_pt];
        TH1F* hcorr_e2c_nonorm[Nbin_jet_pt]; 
        TH1F* hcorr_e2c[Nbin_jet_pt]; 
        TH1F* hcorr_e2c_nounf[Nbin_jet_pt]; 
        TH1F* hcorr_tau[Nbin_jet_pt]; 
        TH1F* hcorr_tau_nounf[Nbin_jet_pt]; 
        
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

        hunfolded_ratio->SetTitle("Purity Corrected Unfolded/Purity Corrected;R_{L};p^{jet}_{T}GeV");
        hunfolded_ratio->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[Nbin_rl_nominal]);
        hunfolded_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        
        if (do_print) 
                c->Print(Form("./plots/unfolded2d_niter%i_ratio_logbinning_inclusive.pdf",niter));

        // Fill the histograms
        for (int bin = 0 ; bin < Nbin_jet_pt ; bin++) {
                hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)         ,"", 1  ,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                hcorr_jet_centroid[bin] = new TH1F(Form("hcorr_jet_centroid%i" ,bin),"", 200,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 

                hcorr_e2c_nonorm[bin]        = new TH1F(Form("hcorr_e2c_nonorm%i",bin)       ,"", Nbin_rl_nominal,rl_nominal_binning );
                hcorr_e2c[bin]          = new TH1F(Form("hcorr_e2c%i",bin)         ,"", Nbin_rl_nominal,rl_nominal_binning );
                hcorr_e2c_nounf[bin]    = new TH1F(Form("hcorr_e2c_nounf%i",bin)   ,"", Nbin_rl_nominal,rl_nominal_binning );
                hcorr_tau[bin]          = new TH1F(Form("hcorr_tau%i",bin)         ,"", Nbin_rl_nominal,tau_nominal_binning);
                hcorr_tau_nounf[bin]    = new TH1F(Form("hcorr_tau_nounf%i",bin)   ,"", Nbin_rl_nominal,tau_nominal_binning);
                
                set_histogram_style(hcorr_e2c[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_tau[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        
                ntuple_jet->Project(Form("hcorr_jet%i" ,bin)         , "jet_pt",jet_full_corr[bin]);
                ntuple_jet->Project(Form("hcorr_jet_centroid%i" ,bin), "jet_pt",jet_full_corr[bin]);

                double jet_pt_centroid = hcorr_jet_centroid[bin]->GetMean();
                for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++) {
                        ntuple_data->GetEntry(entry);

                        if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) 
                                continue;
                        
                        double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt));
                        if (unfolding_weight <= 0) 
                                unfolding_weight = 1;

                        hcorr_e2c_nonorm[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_e2c[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_e2c_nounf[bin]->Fill(R_L,purity*weight_pt/efficiency);
                        hcorr_tau[bin]->Fill(R_L*jet_pt_centroid,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_tau_nounf[bin]->Fill(R_L*jet_pt_centroid,purity*weight_pt/efficiency);
                }

                hcorr_e2c[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_e2c_nounf[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");

                fout->cd();
                hcorr_e2c[bin]->Write();
                hcorr_e2c_nounf[bin]->Write();
                hcorr_tau[bin]->Write();
                hcorr_tau_nounf[bin]->Write();
                gROOT->cd();
        }

        THStack* s_data     = new THStack();
        TLegend* l_data     = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);
        TLegend* l_data_tau = new TLegend(gPad->GetLeftMargin()+0.01,1-gPad->GetTopMargin()-0.01,gPad->GetLeftMargin()+0.21,1-gPad->GetTopMargin()-0.21);

        for (int bin = 0 ; bin < Nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_tau[bin],"E");
                l_data_tau->AddEntry(hcorr_tau[bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L} #LT p^{jet}_{t} #GT(GeV);#Sigma_{EEC}(R_{L})");
        l_data_tau->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

        if (do_print) 
                c->Print(Form("./plots/paircorrtau_niter%i_logbinning_2dunf_inclusive.pdf",niter));

        s_data = new THStack();
        for (int bin = 0 ; bin < Nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_e2c[bin],"E");
                l_data->AddEntry(hcorr_e2c[bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

        if (do_print) 
                c->Print(Form("./plots/paircorre2c_niter%i_logbinning_2dunf_inclusive.pdf",niter));

        if (compare_to_nominal) {
                TFile* fnominal = new TFile((output_folder + namef_histos_paircorr_e2c_logbin).c_str());
                if (fnominal->IsZombie()) 
                        return;

                TH2F* hct_ratio = new TH2F("hct_ratio","",Nbin_rl_nominal,rl_nominal_binning,Nbin_jet_pt,jet_pt_binning);
                TH1F* hnominal[Nbin_jet_pt]; 

                for (int bin = 0 ; bin < Nbin_jet_pt ; bin++) {
                        hnominal[bin] = (TH1F*) fnominal->Get(Form("hcorr_e2c%i",bin));

                        // Normalize both to unity sucha that we can compare the shapes
                        hcorr_e2c[bin]->Scale(1./hcorr_e2c[bin]->Integral());
                        hnominal[bin]->Scale(1./hnominal[bin]->Integral());

                        hcorr_e2c[bin]->Divide(hnominal[bin]);
                        for (int bin_rl = 1 ; bin_rl <= hcorr_e2c[bin]->GetNbinsX() ; bin_rl++) {
                                hct_ratio->SetBinContent(bin_rl, bin + 1, hcorr_e2c[bin]->GetBinContent(bin_rl));
                                hct_ratio->SetBinError(bin_rl, bin + 1, hcorr_e2c[bin]->GetBinError(bin_rl));
                        } 
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

                hct_ratio->SetTitle("Norm. Corr. Pseudodata / Norm. Corr. Data ;R_{L};p^{jet}_{T}GeV");
                hct_ratio->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[Nbin_rl_nominal]);
                hct_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                if (do_print) c->Print(Form("./plots/nom_inclusive_comp_niter%i_ratio_logbinning.pdf",niter));
        }
}
