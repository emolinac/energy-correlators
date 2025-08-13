#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorreec(int niter = nominal_niter)
{
        // Open the necessary files
        TFile* fout = new TFile((output_folder + namef_histos_corr_eec).c_str(),"RECREATE");
        gROOT->cd();

        TFile* fcorr = new TFile((output_folder + namef_ntuple_eec_corr).c_str()); 
        if (fcorr->IsZombie()) 
                return;
        
        TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        // Set the branches of data
        float R_L, jet_pt, weight_pt;
        float efficiency, purity, efficiency_relerror, purity_relerror;
        set_data_ntuple_branches(ntuple_data, &R_L, &jet_pt, &weight_pt, &efficiency, &purity, &efficiency_relerror, &purity_relerror);
        
        // Unfold the purity corrected data
        TFile* f = new TFile((output_folder + namef_ntuple_eec_paircorrections).c_str());
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        
        TH3D* hpurcorr = new TH3D("hpurcorr","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight,weight_binning);
        TH3D* hmeas    = new TH3D("hmeas"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight,weight_binning);
        TH3D* htrue    = new TH3D("htrue"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight,weight_binning);
        RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");
        
        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                // Access entry of ntuple
                ntuple->GetEntry(evt);

                if (abs(R_L_truth - R_L_reco) < rl_resolution) 
                        response->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
        }

        TH2D* hunfolded_ratio     = new TH2D("hunfolded_ratio"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH3D* hpuritycorrected    = new TH3D("hpuritycorrected" ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight,weight_binning);
        TH3D* hpuritycorrected2   = new TH3D("hpuritycorrected2","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight,weight_binning);
        
        ntuple_data->Project("hpuritycorrected" ,"weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
        ntuple_data->Project("hpuritycorrected2","weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
        
        RooUnfoldBayes unfold(response, hpuritycorrected, niter);
        TH3D* hunfolded_bayes            = (TH3D*) unfold.Hunfold();
        TH2D* hunfolded_bayes_rl_jet_pt   = (TH2D*) hunfolded_bayes->Project3D("yx");
        TH2D* hpuritycorrected2_rl_jet_pt = (TH2D*) hpuritycorrected2->Project3D("yx");
        hunfolded_ratio->Divide(hunfolded_bayes_rl_jet_pt,hpuritycorrected2_rl_jet_pt,1,1);
        
        TH1F* hcorr_jet[nbin_jet_pt];
        TH1F* hcorr_jet_centroid[nbin_jet_pt];
        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_eec_nounf[nbin_jet_pt]; 
        TH1F* hcorr_tau[nbin_jet_pt]; 
        TH1F* hcorr_tau_nounf[nbin_jet_pt]; 
        
        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        // gStyle->SetPaintTextFormat("4.2f");
        // hunfolded_ratio->Draw("col text");
        // hunfolded_ratio->SetTitle("Purity Corrected Unfolded/Purity Corrected;R_{L};p_{T,jet} (GeV)");
        // hunfolded_ratio->GetXaxis()->SetRangeUser(rl_min, rl_max);
        // hunfolded_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
        // gPad->SetLogx(0);
        // gPad->SetLogy(1);
        // c->Print(Form("./plots/unfolded3d_niter%_ratio_sepyears_linearbinning_unfbinvarv2.pdf",niter));

        // Fill the histograms
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)  ,"", 1,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                hcorr_jet_centroid[bin] = new TH1F(Form("hcorr_jet_centroid%i" ,bin),"", 200,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 

                hcorr_eec[bin]          = new TH1F(Form("hcorr_eec%i",bin)         ,"", nbin_rl_nominal,rl_nominal_binning);
                hcorr_eec_nounf[bin]    = new TH1F(Form("hcorr_eec_nounf%i",bin)   ,"", nbin_rl_nominal,rl_nominal_binning);
                hcorr_tau[bin]          = new TH1F(Form("hcorr_tau%i",bin)         ,"", nbin_rl_nominal,tau_nominal_binning);
                hcorr_tau_nounf[bin]    = new TH1F(Form("hcorr_tau_nounf%i",bin)   ,"", nbin_rl_nominal,tau_nominal_binning);
                
                set_histogram_style(hcorr_eec[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_tau[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                ntuple_jet->Project(Form("hcorr_jet%i" ,bin)         , "jet_pt",jet_full_corr[bin]);
                ntuple_jet->Project(Form("hcorr_jet_centroid%i" ,bin), "jet_pt",jet_full_corr[bin]);
                
                double jet_pt_centroid = hcorr_jet_centroid[bin]->GetMean();    
                for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++) {
                        ntuple_data->GetEntry(entry);

                        if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) 
                                continue;
                        
                        // if (efficiency_relerror>corr_rel_error) 
                        //         continue;
                        
                        // if (purity_relerror>corr_rel_error) 
                        //         continue;
                        
                        if (efficiency <= 0 || efficiency > 1) 
                                efficiency = 1;
                        
                        if (purity <= 0 || purity > 1) 
                                purity = 1;
                        
                        double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt,weight_pt));
                        if (unfolding_weight <= 0) 
                                unfolding_weight = 1;

                        hcorr_eec[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_eec_nounf[bin]->Fill(R_L,purity*weight_pt/efficiency);
                        hcorr_tau[bin]->Fill(R_L*jet_pt_centroid,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_tau_nounf[bin]->Fill(R_L*jet_pt_centroid,purity*weight_pt/efficiency);
                }

                hcorr_eec[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_eec_nounf[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");

                fout->cd();
                hcorr_eec[bin]->Write();
                hcorr_eec_nounf[bin]->Write();
                hcorr_tau[bin]->Write();
                hcorr_tau_nounf[bin]->Write();
                gROOT->cd();
        }

        THStack* s_data = new THStack();
        TLegend* l_data = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_eec[bin],"E");
                l_data->AddEntry(hcorr_eec[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lf");
        }
        
        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        
        tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

        c->Print(Form("./plots/fullcorreec_niter%i.pdf",niter));
}
