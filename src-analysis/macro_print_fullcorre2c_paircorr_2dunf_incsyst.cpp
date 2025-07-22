#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorre2c_paircorr_2dunf_incsyst(int niter = 4, bool do_print = true)
{
        gStyle->SetPadTopMargin(0.08);

        TFile* fcorr = new TFile((output_folder + namef_ntuple_e2c_paircorr).c_str()); 
        if (fcorr->IsZombie()) 
                return;
        
        TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        // Set the branches of data
        float R_L, jet_pt, weight_pt, efficiency, purity, efficiency_relerror, purity_relerror;
        set_data_ntuple_branches(ntuple_data, &R_L, &jet_pt, &weight_pt, &efficiency, &purity, &efficiency_relerror, &purity_relerror);
        
        // UNFOLDING FIRST
        TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections).c_str());
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        
        // Create histograms with different types of binning
        TH2D* hpurcorr = new TH2D("hpurcorr","",Nbin_R_L_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        TH2D* hmeas    = new TH2D("hmeas"   ,"",Nbin_R_L_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        TH2D* htrue    = new TH2D("htrue"   ,"",Nbin_R_L_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");
        
        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
        {
                // Access entry of ntuple
                ntuple->GetEntry(evt);

                if (abs(R_L_truth-R_L_reco)>0.015) 
                        continue;
        
                response->Fill(R_L_reco, jet_pt_reco, R_L_truth, jet_pt_truth);
        }

        // Fill the purity corrected distributions
        TH2D* hunfolded_ratio     = new TH2D("hunfolded_ratio"  ,"",Nbin_R_L_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        TH2D* hpuritycorrected    = new TH2D("hpuritycorrected" ,"",Nbin_R_L_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        TH2D* hpuritycorrected2   = new TH2D("hpuritycorrected2","",Nbin_R_L_nominal_unfolding,unfolding_rl_nominal_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
        
        ntuple_data->Project("hpuritycorrected" , "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
        ntuple_data->Project("hpuritycorrected2", "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
        
        // Unfold the purity corrected pairs
        RooUnfoldBayes unfold(response, hpuritycorrected, niter);
        TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
        hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected2,1,1);

        TH1F* hcorr_jet[Nbin_jet_pt];
        TH1F* hcorr_jet_centroid[Nbin_jet_pt];
        TH1F* hcorr_e2c[Nbin_jet_pt]; 
        TH1F* hcorr_e2c_syst[Nbin_jet_pt]; 
        TH1F* hcorr_e2c_nounf[Nbin_jet_pt]; 
        TH1F* hcorr_tau[Nbin_jet_pt]; 
        TH1F* hcorr_tau_syst[Nbin_jet_pt]; 
        TH1F* hcorr_tau_nounf[Nbin_jet_pt]; 
        
        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        // Adding content with errors
        TLatex latex;
        latex.SetTextAlign(22); // center alignment
        latex.SetTextSize(0.015);
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
        hunfolded_ratio->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[Nbin_R_L_nominal]);
        hunfolded_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        
        if (do_print) 
                c->Print(Form("./plots/unfolded2d_initer%i_ratio_logbinning_incsyst.pdf",niter));

        // Fill the NOIMNAL histograms
        for (int bin = 0 ; bin < Nbin_jet_pt ; bin++) {
                hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)         ,"", 1  ,jet_pt_binning[bin],jet_pt_binning[bin+1]); 
                hcorr_jet_centroid[bin] = new TH1F(Form("hcorr_jet_centroid%i" ,bin),"", 200,jet_pt_binning[bin],jet_pt_binning[bin+1]); 

                hcorr_e2c[bin]          = new TH1F(Form("hcorr_e2c%i",bin)         ,"", Nbin_R_L_nominal,rl_nominal_binning );
                hcorr_e2c_syst[bin]     = new TH1F(Form("hcorr_e2c_syst%i",bin)    ,"", Nbin_R_L_nominal,rl_nominal_binning );
                hcorr_e2c_nounf[bin]    = new TH1F(Form("hcorr_e2c_nounf%i",bin)   ,"", Nbin_R_L_nominal,rl_nominal_binning );
                hcorr_tau[bin]          = new TH1F(Form("hcorr_tau%i",bin)         ,"", Nbin_R_L_nominal,tau_nominal_binning);
                hcorr_tau_syst[bin]     = new TH1F(Form("hcorr_tau_syst%i",bin)    ,"", Nbin_R_L_nominal,tau_nominal_binning );
                hcorr_tau_nounf[bin]    = new TH1F(Form("hcorr_tau_nounf%i",bin)   ,"", Nbin_R_L_nominal,tau_nominal_binning);

                set_histogram_style(hcorr_e2c[bin]     , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_e2c_syst[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size  );
                set_histogram_style(hcorr_tau[bin]     , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_tau_syst[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size  );
        
                hcorr_e2c[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_e2c_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_tau[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_tau_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                
                ntuple_jet->Project(Form("hcorr_jet%i" ,bin)         , "jet_pt",jet_full_corr[bin]);
                ntuple_jet->Project(Form("hcorr_jet_centroid%i" ,bin), "jet_pt",jet_full_corr[bin]);

                double jet_pt_centroid = hcorr_jet_centroid[bin]->GetMean();
                for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++) {
                        ntuple_data->GetEntry(entry);

                        if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin+1]) 
                                continue;
                        
                        if (efficiency <= 0 || efficiency > 1) 
                                efficiency = 1;
                        
                        if (purity <= 0 || purity > 1) 
                                purity = 1;
                        
                        double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt));
                        if (unfolding_weight <= 0) 
                                unfolding_weight = 1;

                        hcorr_e2c[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_e2c_syst[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_e2c_nounf[bin]->Fill(R_L,purity*weight_pt/efficiency);
                        hcorr_tau[bin]->Fill(R_L*jet_pt_centroid,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_tau_syst[bin]->Fill(R_L*jet_pt_centroid,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_tau_nounf[bin]->Fill(R_L*jet_pt_centroid,purity*weight_pt/efficiency);
                }

                hcorr_e2c[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_e2c_syst[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_e2c_nounf[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                
                hcorr_tau[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_tau_syst[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                
        }

        // Include the systematics in the whole deal
        const int nsyst = sizeof(available_systematics)/sizeof(available_systematics[0]);
        TFile* fsyst[nsyst];
        TH1F* hdev[Nbin_jet_pt];
        TH1F* hdev_tau[Nbin_jet_pt];

        std::cout<<"Source & $20<p_T(jet)<30$ & $30<p_T(jet)<50$ & $50<p_T(jet)<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                fsyst[syst_index] = new TFile((output_folder+devfromnom_namef[available_systematics[syst_index]]).c_str());
                
                if (fsyst[syst_index]->IsZombie()) 
                        continue;
                
                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                
                for (int bin = 0 ; bin < Nbin_jet_pt ; bin++) {
                        hdev[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_deviations%i",bin));

                        set_histo_with_systematics(hdev[bin], hcorr_e2c[bin], hcorr_e2c_syst[bin], systematic_errtype[available_systematics[syst_index]]);

                        if (bin!=Nbin_jet_pt-1) 
                                std::cout<<" & ";
                        else
                                std::cout<<" \\\\ ";
                }

                std::cout<<std::endl;

                delete fsyst[syst_index];
        }

        std::cout<<"Source & $20<p_T(jet)<30$ & $30<p_T(jet)<50$ & $50<p_T(jet)<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                fsyst[syst_index] = new TFile((output_folder+devfromnom_namef[available_systematics[syst_index]]).c_str());
                
                if (fsyst[syst_index]->IsZombie()) 
                        continue;

                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                for (int bin = 0 ; bin < Nbin_jet_pt ; bin++) {
                        hdev_tau[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_deviations_tau%i",bin));

                        set_histo_with_systematics(hdev_tau[bin], hcorr_tau[bin], hcorr_tau_syst[bin], systematic_errtype[available_systematics[syst_index]]);

                        if (bin!=Nbin_jet_pt-1) 
                                std::cout<<" & ";
                        else
                                std::cout<<" \\\\ ";
                }

                std::cout<<std::endl;

                delete fsyst[syst_index];
        }

        THStack* s_data     = new THStack();
        TLegend* l_data     = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);
        TLegend* l_data_tau = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);

        TLatex* lhcbprint = new TLatex();
        
        for (int bin = 0 ; bin < Nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_tau[bin],"E X0");
                s_data->Add(hcorr_tau_syst[bin],"E2");
                l_data_tau->AddEntry(hcorr_tau[bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]),"lpf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L} #LT p^{jet}_{t} #GT(GeV);#Sigma_{EEC}(R_{L})");
        l_data_tau->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        draw_lhcb_tag(lhcbprint);

        if (do_print) 
                c->Print(Form("./plots/paircorrtau_niter%i_logbinning_2dunf_incsyst.pdf",niter));

        s_data = new THStack();
        for (int bin = 0 ; bin < Nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_e2c[bin],"E X0");
                s_data->Add(hcorr_e2c_syst[bin],"E2");
                l_data->AddEntry(hcorr_e2c[bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]),"lpf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        s_data->SetMaximum(0.67);
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);

        draw_lhcb_tag(lhcbprint);

        if (do_print) 
                c->Print(Form("./plots/paircorre2c_niter%i_logbinning_2dunf_incsyst.pdf",niter));
}
