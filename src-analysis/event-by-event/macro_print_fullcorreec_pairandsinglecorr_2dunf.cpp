#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorreec_pairandsinglecorr_2dunf(int niter = nominal_niter, bool do_print = true, bool do_jet_unfolding = false, bool apply_alice_factor = false)
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + namef_histos_pairandsinglecorr_eec).c_str(),"RECREATE");
        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_ntuple_eec_pairandsinglecorr).c_str()); 
        if (fcorr->IsZombie()) 
                return;
        
        TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        // Set the branches of data
        float R_L, jet_pt, weight_pt, event_weight, efficiency, purity, efficiency_relerror, purity_relerror, eq_charge;
        float h1_efficiency, h2_efficiency, h1_purity, h2_purity;
        set_data_ntuple_branches(ntuple_data, &event_weight, &R_L, &jet_pt, &weight_pt, &efficiency, &purity, &efficiency_relerror, &purity_relerror, &eq_charge);
        ntuple_data->SetBranchAddress("h1_efficiency",&h1_efficiency);
        ntuple_data->SetBranchAddress("h2_efficiency",&h2_efficiency);
        ntuple_data->SetBranchAddress("h1_purity",&h1_purity);
        ntuple_data->SetBranchAddress("h2_purity",&h2_purity);
        
        // Unfold the purity corrected pairs
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

        TH2D* hunfolded_ratio      = new TH2D("hunfolded_ratio"     ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hpuritycorrected     = new TH2D("hpuritycorrected"    ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hpuritycorrected_ref = new TH2D("hpuritycorrected_ref","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        ntuple_data->Project("hpuritycorrected" , "jet_pt:R_L","purity");
        ntuple_data->Project("hpuritycorrected_ref", "jet_pt:R_L","purity");
        
        RooUnfoldBayes unfold(response, hpuritycorrected, niter);

        TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
        
        hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected_ref,1,1);
        hunfolded_ratio->Smooth();

        // Unfold the purity corrected jets
        TNtuple* ntuple_jet_unfolding = (TNtuple*) f->Get(name_ntuple_mcreco_jet.c_str());
        
        float jet_pt_unfolding_reco, jet_pt_unfolding_truth;
        set_unfolding_jet_ntuple_branches(ntuple_jet_unfolding, &jet_pt_unfolding_reco, &jet_pt_unfolding_truth);
        
        TH1D* hpurcorr_jet = new TH1D("hpurcorr_jet","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH1D* hmeas_jet    = new TH1D("hmeas_jet"   ,"",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH1D* htrue_jet    = new TH1D("htrue_jet"   ,"",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);

        TH2D* hresponse_jet = new TH2D("hresponse_jet","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        for (int evt = 0 ; evt < ntuple_jet_unfolding->GetEntries() ; evt++) {
                ntuple_jet_unfolding->GetEntry(evt);

                hresponse_jet->Fill(jet_pt_unfolding_reco, jet_pt_unfolding_truth);
        }

        RooUnfoldResponse* response_jet = new RooUnfoldResponse(hmeas_jet, htrue_jet, hresponse_jet, "response_jet");
        
        TH1D* hunfolded_ratio_jet      = new TH1D("hunfolded_ratio_jet"  ,"",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH1D* hpuritycorrected_jet     = new TH1D("hpuritycorrected_jet" ,"",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH1D* hpuritycorrected_ref_jet = new TH1D("hpuritycorrected_ref_jet","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        ntuple_jet->Project("hpuritycorrected_jet" , "jet_pt", "jet_purity");
        ntuple_jet->Project("hpuritycorrected_ref_jet", "jet_pt", "jet_purity");
        
        RooUnfoldBayes unfold_jet(response_jet, hpuritycorrected_jet, 2);

        TH1D* hunfolded_bayes_jet = (TH1D*) unfold_jet.Hunfold();
        
        hunfolded_ratio_jet->Divide(hunfolded_bayes_jet,hpuritycorrected_ref_jet,1,1);
        hunfolded_ratio_jet->Smooth();

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
                c->Print(Form("./plots/unfolded2d_unf-niter%i_ratio_allcorr.pdf",niter));

        hunfolded_ratio_jet->Draw("col");
        hunfolded_ratio_jet->SetTitle("Jet : Purity Corrected Unfolded/Purity Corrected;p_{T,jet} (GeV);");
        hunfolded_ratio_jet->GetXaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        if (do_print) 
                c->Print(Form("./plots/unfolded1d_jetpt_unf-niter%i_ratio_allcorr.pdf",niter));

        hresponse_jet->Draw("col text");
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        if (do_print) 
                c->Print("./plots/jet_pt_response_matrix_allcorr.pdf");

        // Fill the histograms
        TH1F* hcorr_jet[nbin_jet_pt];
        TH1F* hcorr_jet_centroid[nbin_jet_pt];
        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_npair[nbin_jet_pt]; 
        TH1F* hcorr_tau[nbin_jet_pt]; 
        
        TH1F* hcorr_eec_eqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_total[nbin_jet_pt]; // necessary due to the difference in the type of binning
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)         ,"", 1  ,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                hcorr_jet_centroid[bin] = new TH1F(Form("hcorr_jet_centroid%i" ,bin),"", 200,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 

                hcorr_eec[bin]          = new TH1F(Form("hcorr_eec%i",bin)         ,"", nbin_rl_nominal,rl_nominal_binning );
                hcorr_tau[bin]          = new TH1F(Form("hcorr_tau%i",bin)         ,"", nbin_rl_nominal,tau_nominal_binning);
                hcorr_npair[bin]        = new TH1F(Form("hcorr_npair%i",bin)       ,"", nbin_rl_nominal,rl_nominal_binning );

                hcorr_eec_eqcharge[bin]  = new TH1F(Form("hcorr_eec_eqcharge%i",bin) ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                hcorr_eec_neqcharge[bin] = new TH1F(Form("hcorr_eec_neqcharge%i",bin),"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                hcorr_eec_total[bin]     = new TH1F(Form("hcorr_eec_total%i",bin)    ,"",nbin_chargedeec_nominal,rl_chargedeec_binning);
                
                set_histogram_style(hcorr_eec[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_tau[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_npair[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                set_histogram_style(hcorr_eec_eqcharge[bin] , std_marker_color_jet_pt[bin] , std_line_width, std_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_eec_neqcharge[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        
                ntuple_jet->Project(Form("hcorr_jet%i" ,bin)         , "jet_pt",jet_full_corr[bin]);
                ntuple_jet->Project(Form("hcorr_jet_centroid%i" ,bin), "jet_pt",jet_full_corr[bin]);

                // Apply the unfolding factor
                if (do_jet_unfolding){
                        hcorr_jet[bin]->Scale(hunfolded_ratio_jet->GetBinContent(bin+2));

                        std::cout<<"Scaling by "<<hunfolded_ratio_jet->GetBinContent(bin+2)<<std::endl;
                }

                double jet_pt_centroid = hcorr_jet_centroid[bin]->GetMean();
                for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++) {
                        ntuple_data->GetEntry(entry);

                        if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) 
                                continue;

                        if(h1_efficiency <= 0 ||h2_efficiency <= 0)
                                continue;
                        
                        if(h1_purity <= 0 ||h2_purity <= 0)
                                continue;
                        
                        // double unfolding_weight = 1.;//hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt));
                        // // if (unfolding_weight <= 0) 
                        // //         unfolding_weight = 1;

                        double correction = h1_purity*h2_purity/efficiency;

                        hcorr_eec[bin]->Fill(R_L,event_weight*weight_pt*correction);
                        hcorr_tau[bin]->Fill(R_L*jet_pt_centroid,event_weight*weight_pt*correction);
                        hcorr_npair[bin]->Fill(R_L,event_weight*correction);
                        
                        if (eq_charge > 0)
                                hcorr_eec_eqcharge[bin]->Fill(R_L,event_weight*weight_pt*correction);
                        else if (eq_charge < 0)
                                hcorr_eec_neqcharge[bin]->Fill(R_L,event_weight*weight_pt*correction);
                }

                hcorr_tau[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_eec[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_npair[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");

                hcorr_eec_total[bin]->Add(hcorr_eec_eqcharge[bin],hcorr_eec_neqcharge[bin],1,1);
                hcorr_eec_eqcharge[bin]->Divide(hcorr_eec_total[bin]);
                hcorr_eec_neqcharge[bin]->Divide(hcorr_eec_total[bin]);

                fout->cd();
                hcorr_eec[bin]->Write();
                hcorr_tau[bin]->Write();
                hcorr_eec_eqcharge[bin]->Write();
                hcorr_eec_neqcharge[bin]->Write();
                gROOT->cd();
        }

        THStack* s_data     = new THStack();
        TLegend* l_data     = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);
        TLegend* l_data_tau = new TLegend(gPad->GetLeftMargin()+0.01,1-gPad->GetTopMargin()-0.01,gPad->GetLeftMargin()+0.21,1-gPad->GetTopMargin()-0.21);

        TLegend* l_data_chargedeec;

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_tau[bin],"E");
                l_data_tau->AddEntry(hcorr_tau[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L} #LT p_{T,jet} #GT(GeV);#Sigma_{EEC}(R_{L})");
        l_data_tau->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

        if (do_print) 
                c->Print(Form("./plots/corrtau_unf-niter%i_jetptunf-%s_2dunf_allcorr.pdf",niter,(do_jet_unfolding)?"yes":"no"));

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data = new THStack();
                s_data->Add(hcorr_eec_eqcharge[bin] ,"E1 ");
                s_data->Add(hcorr_eec_neqcharge[bin],"E1 ");

                
                l_data_chargedeec = new TLegend(1-gPad->GetRightMargin()-0.21,gPad->GetBottomMargin()+0.01,1-gPad->GetRightMargin()-0.01,gPad->GetBottomMargin()+0.20);
                l_data_chargedeec->SetHeader(Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l_data_chargedeec->AddEntry(hcorr_eec_eqcharge[bin], "eq. charge correlations","lfp");
                l_data_chargedeec->AddEntry(hcorr_eec_neqcharge[bin],"op. charge correlations","lfp");

                s_data->Draw("NOSTACK");
                s_data->SetMaximum(0.8);
                s_data->SetMinimum(0.2);
                s_data->SetTitle(";R_{L};Charged #Sigma_{EEC}(R_{L})");
                l_data_chargedeec->Draw("SAME");
                gPad->SetLogx(1);
                gPad->SetLogy(0);
                
                tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

                if (do_print)
                        c->Print(Form("./plots/corrchargedeec_jetptbin%i_unf-niter%i_jetptunf-%s_2dunf_allcorr.pdf", bin, niter,(do_jet_unfolding)?"yes":"no"));
        }

        s_data = new THStack();
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_eec[bin],"E");
                l_data->AddEntry(hcorr_eec[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        // s_data->SetMaximum(1.2);
        // s_data->SetMaximum(0);
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

        if (do_print) 
                c->Print(Form("./plots/correec_unf-niter%i_jetptunf-%s_2dunf_allcorr.pdf",niter,(do_jet_unfolding)?"yes":"no"));
}
