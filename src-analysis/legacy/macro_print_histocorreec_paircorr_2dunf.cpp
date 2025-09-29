#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_histocorreec_paircorr_2dunf(int niter = nominal_niter, bool do_print = true, bool do_jet_unfolding = false)
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + namef_histopaircorr_eec).c_str(),"RECREATE");
        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_paircorr_histos).c_str()); 
        if (fcorr->IsZombie()) 
                return;

        TH2F* h_eec        = (TH2F*) fcorr->Get("h_eec");
        TH2F* h_efficiency = (TH2F*) fcorr->Get("hefficiency");
        TH2F* h_purity     = (TH2F*) fcorr->Get("hpurity");

        TH1F* h_njet           = (TH1F*) fcorr->Get("h_njet");
        TH1F* h_njet_wmuoneff  = (TH1F*) fcorr->Get("h_njet_wmuoneff");
        TH1F* h_efficiency_jet = (TH1F*) fcorr->Get("hefficiency_jet");
        TH1F* h_purity_jet     = (TH1F*) fcorr->Get("hpurity_jet");

        TH1F* h_muon_eff = new TH1F("h_muon_eff","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        h_muon_eff->Divide(h_njet, h_njet_wmuoneff);
        
        // Correct the jets
        TFile* f = new TFile((output_folder + namef_ntuple_eec_paircorrections).c_str());

        TNtuple* ntuple_jet_unfolding = (TNtuple*) f->Get(name_ntuple_mcreco_jet.c_str());
        
        float jet_pt_unfolding_reco, jet_pt_unfolding_truth;
        set_unfolding_jet_ntuple_branches(ntuple_jet_unfolding, &jet_pt_unfolding_reco, &jet_pt_unfolding_truth);
        
        TH1D* hmeas_jet = new TH1D("hmeas_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1D* htrue_jet = new TH1D("htrue_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2D* hresp_jet = new TH2D("hresp_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        
        for (int evt = 0 ; evt < ntuple_jet_unfolding->GetEntries() ; evt++) {
                ntuple_jet_unfolding->GetEntry(evt);

                hresp_jet->Fill(jet_pt_unfolding_reco, jet_pt_unfolding_truth);
        }

        RooUnfoldResponse* response_jet = new RooUnfoldResponse(hmeas_jet, htrue_jet, hresp_jet, "response_jet");
        
        TH1F* h_njet_purity_corrected = new TH1F("h_njet_purity_corrected","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        h_njet_purity_corrected->Multiply(h_njet,h_purity_jet,1,1);

        RooUnfoldBayes unfold_jet(response_jet, h_njet_purity_corrected, 4);

        TH1D* h_njet_unfolded = (TH1D*) unfold_jet.Hunfold();

        h_njet_unfolded->Divide(h_efficiency_jet);

        h_njet_unfolded->Divide(h_muon_eff);
        
        // Correct the eecs
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());
        
        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        
        TH2D* hmeas = new TH2D("hmeas" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH2D* htrue = new TH2D("htrue" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);

        RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");
        
        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                ntuple->GetEntry(evt);

                if (abs(R_L_truth - R_L_reco) < rl_resolution) 
                        response->Fill(R_L_reco, jet_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
        }

        TH2F* h_eec_purity_corrected = new TH2F("h_eec_purity_corrected","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        h_eec_purity_corrected->Multiply(h_eec,h_purity,1,1);
        
        RooUnfoldBayes unfold(response, h_eec_purity_corrected, niter);

        TH2D* h_eec_unfolded = (TH2D*) unfold.Hunfold();

        h_eec_unfolded->Divide(h_efficiency);

        for (int i = 1 ; i <= h_eec_unfolded->GetNbinsX(); i++) {
                for (int j = 1 ; j <= h_eec_unfolded->GetNbinsY(); j++) {
                        double reweight = h_purity_jet->GetBinContent(j)/h_efficiency_jet->GetBinContent(j)/h_muon_eff->GetBinContent(j);

                        h_eec_unfolded->SetBinContent(i, j, h_eec_unfolded->GetBinContent(i, j) * reweight);
                        h_eec_unfolded->SetBinError(i, j, h_eec_unfolded->GetBinError(i, j) * reweight);
                }
        }

        // Fill the histograms
        TH1F* hcorr_jet_centroid[nbin_jet_pt];
        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_npair[nbin_jet_pt]; 
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_eec[bin]   = new TH1F(Form("hcorr_eec%i",bin)   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hcorr_npair[bin] = new TH1F(Form("hcorr_npair%i",bin) , "", nbin_rl_nominal, rl_nominal_binning );
                
                set_histogram_style(hcorr_eec[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_npair[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                for (int i = 1 ; i <= hcorr_eec[bin]->GetNbinsX(); i++) {
                        if (i == 1 || i == hcorr_eec[bin]->GetNbinsX()) {
                                hcorr_eec[bin]->SetBinContent(i, 0);
                                hcorr_eec[bin]->SetBinError(i, 0);        

                                continue;
                        }

                        hcorr_eec[bin]->SetBinContent(i, h_eec_unfolded->GetBinContent(i, bin + 2));
                        hcorr_eec[bin]->SetBinError(i, h_eec_unfolded->GetBinError(i, bin + 2));
                }        

                hcorr_eec[bin]->Scale(1./h_njet_unfolded->GetBinContent(bin + 2),"width");

                fout->cd();
                hcorr_eec[bin]->Write();
                gROOT->cd();
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        THStack* s_data     = new THStack();
        TLegend* l_data     = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_eec[bin],"E");
                l_data->AddEntry(hcorr_eec[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        s_data->SetMaximum(1.2);
        s_data->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        
        tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

        if (do_print) 
                c->Print(Form("./plots/correec_unf-niter%i_jetptunf-%s_2dunf_avgeeventweight.pdf",niter,(do_jet_unfolding)?"yes":"no"));
}
