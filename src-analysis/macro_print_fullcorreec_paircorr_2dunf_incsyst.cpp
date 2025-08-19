#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorreec_paircorr_2dunf_incsyst(int niter = nominal_niter, bool do_print = true)
{
        gStyle->SetPadTopMargin(0.08);

        TFile* fcorr = new TFile((output_folder + namef_ntuple_eec_paircorr).c_str()); 
        if (fcorr->IsZombie()) 
                return;
        
        TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        // Set the branches of data
        float R_L, jet_pt, weight_pt, event_weight, efficiency, purity, efficiency_relerror, purity_relerror, eq_charge;
        set_data_ntuple_branches(ntuple_data, &event_weight, &R_L, &jet_pt, &weight_pt, &efficiency, &purity, &efficiency_relerror, &purity_relerror, &eq_charge);
        
        // Unfold the purity corrected data
        TFile* f = new TFile((output_folder + namef_ntuple_eec_paircorrections).c_str());
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

        TH2D* hunfolded_ratio      = new TH2D("hunfolded_ratio"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hpuritycorrected     = new TH2D("hpuritycorrected" ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hpuritycorrected_ref = new TH2D("hpuritycorrected_ref","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        ntuple_data->Project("hpuritycorrected"    ,"jet_pt:R_L","purity");
        ntuple_data->Project("hpuritycorrected_ref","jet_pt:R_L","purity");
        
        RooUnfoldBayes unfold(response, hpuritycorrected, niter);

        TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
        
        hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected_ref,1,1);
        hunfolded_ratio->Smooth();

        TH1F* hcorr_jet[nbin_jet_pt];
        TH1F* hcorr_jet_centroid[nbin_jet_pt];
        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_eec_syst[nbin_jet_pt]; 
        TH1F* hcorr_tau[nbin_jet_pt]; 
        TH1F* hcorr_tau_syst[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge_syst[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge_syst[nbin_jet_pt]; 
        TH1F* hcorr_eec_total[nbin_jet_pt]; 
        
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
        
        if (do_print) 
                c->Print(Form("./plots/unfolded2d_unf-niter%i_ratio_incsyst.pdf",niter));

        // Fill the NOIMNAL histograms
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)         ,"", 1  ,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                hcorr_jet_centroid[bin] = new TH1F(Form("hcorr_jet_centroid%i" ,bin),"", 200,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 

                hcorr_eec[bin]          = new TH1F(Form("hcorr_eec%i",bin)         ,"", nbin_rl_nominal,rl_nominal_binning );
                hcorr_eec_syst[bin]     = new TH1F(Form("hcorr_eec_syst%i",bin)    ,"", nbin_rl_nominal,rl_nominal_binning );
                hcorr_tau[bin]          = new TH1F(Form("hcorr_tau%i",bin)         ,"", nbin_rl_nominal,tau_nominal_binning);
                hcorr_tau_syst[bin]     = new TH1F(Form("hcorr_tau_syst%i",bin)    ,"", nbin_rl_nominal,tau_nominal_binning );
                hcorr_eec_eqcharge[bin]       = new TH1F(Form("hcorr_eec_eqcharge%i",bin)      ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                hcorr_eec_eqcharge_syst[bin]  = new TH1F(Form("hcorr_eec_eqcharge_syst%i",bin) ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                hcorr_eec_neqcharge[bin]      = new TH1F(Form("hcorr_eec_neqcharge%i",bin)     ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                hcorr_eec_neqcharge_syst[bin] = new TH1F(Form("hcorr_eec_neqcharge_syst%i",bin),"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                hcorr_eec_total[bin]          = new TH1F(Form("hcorr_eec_total%i",bin)    ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                
                set_histogram_style(hcorr_eec[bin]               , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_eec_syst[bin]          , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size  );
                set_histogram_style(hcorr_tau[bin]               , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_tau_syst[bin]          , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size  );
                set_histogram_style(hcorr_eec_eqcharge[bin]      , corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size+1);
                set_histogram_style(hcorr_eec_eqcharge_syst[bin] , corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size+1);
                set_histogram_style(hcorr_eec_neqcharge[bin]     , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_eec_neqcharge_syst[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size  );
                
                hcorr_eec[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_tau[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_tau_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_eqcharge[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_eqcharge_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_neqcharge[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_neqcharge_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                
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

                        hcorr_eec[bin]->Fill(R_L,event_weight*purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_eec_syst[bin]->Fill(R_L,event_weight*purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_tau[bin]->Fill(R_L*jet_pt_centroid,event_weight*purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_tau_syst[bin]->Fill(R_L*jet_pt_centroid,event_weight*purity*unfolding_weight*weight_pt/efficiency);

                        if (eq_charge > 0) {
                                hcorr_eec_eqcharge[bin]->Fill(R_L,event_weight*purity*unfolding_weight*weight_pt/efficiency);
                                hcorr_eec_eqcharge_syst[bin]->Fill(R_L,event_weight*purity*unfolding_weight*weight_pt/efficiency);
                        }
                        else if (eq_charge < 0) {
                                hcorr_eec_neqcharge[bin]->Fill(R_L,event_weight*purity*unfolding_weight*weight_pt/efficiency);
                                hcorr_eec_neqcharge_syst[bin]->Fill(R_L,event_weight*purity*unfolding_weight*weight_pt/efficiency);
                        }
                }

                hcorr_eec_total[bin]->Add(hcorr_eec_eqcharge[bin],hcorr_eec_neqcharge[bin],1,1);
                hcorr_eec_eqcharge[bin]->Divide(hcorr_eec_total[bin]);
                hcorr_eec_neqcharge[bin]->Divide(hcorr_eec_total[bin]);
                hcorr_eec_eqcharge_syst[bin]->Divide(hcorr_eec_total[bin]);
                hcorr_eec_neqcharge_syst[bin]->Divide(hcorr_eec_total[bin]);
                
                hcorr_eec[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_eec_syst[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                
                hcorr_tau[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_tau_syst[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                
        }

        // Include the systematics in the whole deal
        const int nsyst = sizeof(available_systematics)/sizeof(available_systematics[0]);
        TFile* fsyst[nsyst];
        TH1F* hdev_eec[nbin_jet_pt];
        TH1F* hdev_tau[nbin_jet_pt];
        TH1F* hdev_eec_eqcharge[nbin_jet_pt];
        TH1F* hdev_eec_neqcharge[nbin_jet_pt];

        std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                fsyst[syst_index] = new TFile((output_folder + devfromnom_namef[available_systematics[syst_index]]).c_str());
                
                if (fsyst[syst_index]->IsZombie()) 
                        continue;
                
                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hdev_eec[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_deviations_eec%i",bin));

                        set_histo_with_systematics(hdev_eec[bin], hcorr_eec[bin], hcorr_eec_syst[bin], systematic_errtype[available_systematics[syst_index]]);

                        if (bin!=nbin_jet_pt-1) 
                                std::cout<<" & ";
                        else
                                std::cout<<" \\\\ ";
                }

                std::cout<<std::endl;

                delete fsyst[syst_index];
        }

        std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                fsyst[syst_index] = new TFile((output_folder + devfromnom_namef[available_systematics[syst_index]]).c_str());
                
                if (fsyst[syst_index]->IsZombie()) 
                        continue;

                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hdev_tau[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_deviations_tau%i",bin));

                        set_histo_with_systematics(hdev_tau[bin], hcorr_tau[bin], hcorr_tau_syst[bin], systematic_errtype[available_systematics[syst_index]]);

                        if (bin!=nbin_jet_pt-1) 
                                std::cout<<" & ";
                        else
                                std::cout<<" \\\\ ";
                }

                std::cout<<std::endl;

                delete fsyst[syst_index];
        }

        std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                fsyst[syst_index] = new TFile((output_folder + devfromnom_namef[available_systematics[syst_index]]).c_str());
                
                if (fsyst[syst_index]->IsZombie()) 
                        continue;
                
                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hdev_eec_eqcharge[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_deviations_eec_eqcharge%i",bin));

                        set_histo_with_systematics(hdev_eec_eqcharge[bin], hcorr_eec_eqcharge[bin], hcorr_eec_eqcharge_syst[bin], systematic_errtype[available_systematics[syst_index]]);

                        if (bin!=nbin_jet_pt-1) 
                                std::cout<<" & ";
                        else
                                std::cout<<" \\\\ ";
                }

                std::cout<<std::endl;

                delete fsyst[syst_index];
        }

        std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                fsyst[syst_index] = new TFile((output_folder + devfromnom_namef[available_systematics[syst_index]]).c_str());
                
                if (fsyst[syst_index]->IsZombie()) 
                        continue;
                
                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hdev_eec_neqcharge[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_deviations_eec_neqcharge%i",bin));

                        set_histo_with_systematics(hdev_eec_neqcharge[bin], hcorr_eec_neqcharge[bin], hcorr_eec_neqcharge_syst[bin], systematic_errtype[available_systematics[syst_index]]);

                        if (bin!=nbin_jet_pt-1) 
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
        TLegend* l_data_chargedeec;

        TLatex* lhcbprint = new TLatex();
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_tau[bin],"E X0");
                s_data->Add(hcorr_tau_syst[bin],"E2");
                l_data_tau->AddEntry(hcorr_tau[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lpf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L} #LT p_{T,jet} #GT(GeV);#Sigma_{EEC}(R_{L})");
        l_data_tau->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        draw_lhcb_tag(lhcbprint);

        if (do_print) 
                c->Print(Form("./plots/corrtau_unf-niter%i_2dunf_incsyst.pdf",niter));

        s_data = new THStack();
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_eec[bin],"E X0");
                s_data->Add(hcorr_eec_syst[bin],"E2");
                l_data->AddEntry(hcorr_eec[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lpf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        s_data->SetMaximum(1.05);
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);

        draw_lhcb_tag(lhcbprint);

        if (do_print) 
                c->Print(Form("./plots/correec_unf-niter%i_2dunf_incsyst.pdf",niter));

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data = new THStack();
                s_data->Add(hcorr_eec_eqcharge[bin] ,"E X0");
                s_data->Add(hcorr_eec_eqcharge_syst[bin],"E2");
                s_data->Add(hcorr_eec_neqcharge[bin] ,"E X0");
                s_data->Add(hcorr_eec_neqcharge_syst[bin],"E2");

                l_data_chargedeec = new TLegend(1-gPad->GetRightMargin()-0.21,gPad->GetBottomMargin()+0.01,1-gPad->GetRightMargin()-0.01,gPad->GetBottomMargin()+0.20);
                l_data_chargedeec->SetHeader(Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l_data_chargedeec->AddEntry(hcorr_eec_eqcharge[bin], "eq. charge correlations","lfp");
                l_data_chargedeec->AddEntry(hcorr_eec_neqcharge[bin],"op. charge correlations","lfp");

                s_data->Draw("NOSTACK");
                s_data->SetMaximum(1);
                s_data->SetMinimum(0);
                s_data->SetTitle(";R_{L};Charged #Sigma_{EEC}(R_{L})");
                l_data_chargedeec->Draw("SAME");
                gPad->SetLogx(1);
                gPad->SetLogy(0);
                
                draw_lhcb_tag(lhcbprint);

                if (do_print)
                        c->Print(Form("./plots/corrchargedeec_jetptbin%i_unf-niter%i_2dunf_incsyst.pdf",bin,niter));
        }
}
