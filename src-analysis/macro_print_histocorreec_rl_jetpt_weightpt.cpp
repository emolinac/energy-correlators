#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_histocorreec_rl_jetpt_weightpt(int niter = 4, int niter_jet = 4)
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niter%i_niterjet%i.root",niter,niter_jet)).c_str(),"RECREATE");
        TFile* f    = new TFile((output_folder + namef_ntuple_eec_paircorrections).c_str());
        TFile* fjet = new TFile((output_folder + namef_ntuple_jet_purity).c_str());
        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_3dpaircorr_rl_jetpt_weightpt_histos).c_str()); 
        if (fcorr->IsZombie()) 
                return;

        TH3F* h_npair            = (TH3F*) fcorr->Get("h_npair");
        TH3F* h_npair_wmuon      = (TH3F*) fcorr->Get("h_npair_wmuon");
        TH3F* h_eqchnpair        = (TH3F*) fcorr->Get("h_eqchnpair");
        TH3F* h_eqchnpair_wmuon  = (TH3F*) fcorr->Get("h_eqchnpair_wmuon");
        TH3F* h_neqchnpair       = (TH3F*) fcorr->Get("h_neqchnpair");
        TH3F* h_neqchnpair_wmuon = (TH3F*) fcorr->Get("h_neqchnpair_wmuon");
        
        TH3F* h_efficiency  = (TH3F*) fcorr->Get("hefficiency");
        TH3F* h_purity      = (TH3F*) fcorr->Get("hpurity");

        TH1F* h_njet           = (TH1F*) fcorr->Get("h_njet");
        TH1F* h_njet_wmuoneff  = (TH1F*) fcorr->Get("h_njet_wmuoneff");
        TH1F* h_efficiency_jet = (TH1F*) fcorr->Get("hefficiency_jet");
        TH1F* h_purity_jet     = (TH1F*) fcorr->Get("hpurity_jet");

        // Correct the jets
        TNtuple* ntuple_jet_unfolding = (TNtuple*) fjet->Get(name_ntuple_jetpurity.c_str());
        
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
        h_njet_purity_corrected->Multiply(h_njet_wmuoneff,h_purity_jet,1,1);

        RooUnfoldBayes unfold_jet(response_jet, h_njet_purity_corrected, niter_jet);

        TH1D* h_njet_unfolded = (TH1D*) unfold_jet.Hreco();

        h_njet_unfolded->Divide(h_efficiency_jet);
        
        // Correct the npairs
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());
        
        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth, h1_pt_reco, h1_pt_truth, h2_pt_reco, h2_pt_truth, eq_charge_reco;
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

        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                ntuple->GetEntry(evt);

                if (R_L_truth != -999)
                        response_npair->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
                if (R_L_truth != -999 && eq_charge_reco > 0)
                        response_eqchnpair->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
                if (R_L_truth != -999 && eq_charge_reco < 0)
                        response_neqchnpair->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
        }

        TH3F* h_npair_purity_corrected      = new TH3F("h_npair_purity_corrected",     "",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_eqchnpair_purity_corrected  = new TH3F("h_eqchnpair_purity_corrected", "",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_neqchnpair_purity_corrected = new TH3F("h_neqchnpair_purity_corrected","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        
        h_npair_purity_corrected->Multiply(h_npair_wmuon,h_purity,1,1);
        h_eqchnpair_purity_corrected->Multiply(h_eqchnpair_wmuon,h_purity,1,1);
        h_neqchnpair_purity_corrected->Multiply(h_neqchnpair_wmuon,h_purity,1,1);
        
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

        // Plot the EECs
        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_eqcheec[nbin_jet_pt]; 
        TH1F* hcorr_neqcheec[nbin_jet_pt]; 
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                int nominal_jet_pt_bin = bin + 3;

                hcorr_eec[bin]      = new TH1F(Form("hcorr_eec%i",bin)     , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hcorr_eqcheec[bin]  = new TH1F(Form("hcorr_eqcheec%i",bin) , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hcorr_neqcheec[bin] = new TH1F(Form("hcorr_neqcheec%i",bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);

                set_histogram_style(hcorr_eec[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_eqcheec[bin], corr_marker_color_jet_pt[bin], std_line_width-1, std_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_neqcheec[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+2);
                
                project_nominal_phase_space(h_eec_unfolded_2d     , hcorr_eec[bin]     , nominal_jet_pt_bin);
                project_nominal_phase_space(h_eqcheec_unfolded_2d , hcorr_eqcheec[bin] , nominal_jet_pt_bin);
                project_nominal_phase_space(h_neqcheec_unfolded_2d, hcorr_neqcheec[bin], nominal_jet_pt_bin);
                
                hcorr_eqcheec[bin]->Divide(hcorr_eec[bin]);
                hcorr_neqcheec[bin]->Divide(hcorr_eec[bin]);

                hcorr_eec[bin]->Scale(1./h_njet_unfolded->GetBinContent(bin + 3),"width");
        
                fout->cd();
                hcorr_eec[bin]->Write();
                hcorr_eqcheec[bin]->Write();
                hcorr_neqcheec[bin]->Write();
                gROOT->cd();
        }

        // Plot the EECs AS A FUNCTION OF R_l * <p^2_{T,jet}>TH1F* hcorr_tau[nbin_jet_pt];
        TH1F* hcorr_tau[nbin_jet_pt];        
        double tau_binning[nbin_jet_pt][nbin_rl_nominal + 1];
        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                double avge_pt2_jet = (jet_pt_binning[i + 1] + jet_pt_binning[i])/2.;
                
                get_tau_binning_from_eec_binning(tau_binning[i], rl_nominal_binning, avge_pt2_jet);
                
                hcorr_tau[i] = new TH1F(Form("hcorr_tau%i",i),"",nbin_rl_nominal,tau_binning[i]);
                get_tau_from_uoflow_eec(hcorr_eec[i], hcorr_tau[i], avge_pt2_jet);
                
                set_histogram_style(hcorr_tau[i], corr_marker_color_jet_pt[i], std_line_width, corr_marker_style_jet_pt[i], std_marker_size + 1);

                fout->cd();
                hcorr_tau[i]->Write();
                gROOT->cd();
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1480);
        c->Draw();

        TLatex* tex = new TLatex();
        TLatex* lhcbprint = new TLatex();
        set_lhcb_watermark_properties(tex);

        // Print EEC
        THStack* s_data = new THStack();
        TLegend* l_data = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_eec[bin],"EP");
                l_data->AddEntry(hcorr_eec[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"p");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        s_data->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
        s_data->SetMaximum(1.3);
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");
        draw_lhcb_tag(lhcbprint);

        c->Print(Form("./plots/correec_rl_jetpt_weightpt_unf-niter%i.pdf",niter));

        // Print EEC as a function of R_l * <p^2_{T,jet}>
        s_data = new THStack();
        l_data = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_tau[bin],"EP");
                l_data->AddEntry(hcorr_tau[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"p");
        }

        TH1F* frame = gPad->DrawFrame(tau_binning[0][0], 0.0, tau_binning[2][nbin_rl_nominal], 0.08);
        
        s_data->Draw("NOSTACK SAME");
        frame->SetTitle(";R_{L}#LT p_{T,jet} #GT;#Sigma_{EEC}(R_{L})#times ln(#LT p_{T,jet} #GT)/#LT p_{T,jet} #GT");
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

        draw_lhcb_tag(lhcbprint);

        c->Print(Form("./plots/corrtau_rl_jetpt_weightpt_unf-niter%i.pdf",niter));

        // Print the charged EECs
        s_data = new THStack();
        l_data = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_eqcheec[bin],"EP X0");
                s_data->Add(hcorr_eqcheec[bin],"HIST L");
                s_data->Add(hcorr_neqcheec[bin],"EP X0 C");
                s_data->Add(hcorr_neqcheec[bin],"HIST L");
                l_data->AddEntry(hcorr_eqcheec[bin],Form("eq. charge: %.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"p");
                l_data->AddEntry(hcorr_neqcheec[bin],Form("neq. charge: %.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"p");
        }

        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L};");
        s_data->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
        s_data->SetMaximum(1.05);
        s_data->SetMinimum(0.21);
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

        draw_lhcb_tag(lhcbprint);

        c->Print(Form("./plots/corrchargedeec_rl_jetpt_weightpt_unf-niter%i.pdf",niter));
}
