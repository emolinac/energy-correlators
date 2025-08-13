#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorreec_data2mc(int niter = 4, bool apply_unfolding = true, bool apply_jet_unfolding = true)
{
        gStyle->SetPadTopMargin(0.08);

        TFile* fcorr = new TFile((output_folder + namef_ntuple_eec_paircorr).c_str()); 
        if (fcorr->IsZombie()) 
                return;
        
        TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        float R_L, efficiency, purity, jet_pt, weight, weight_pt;
        
        ntuple_data->SetBranchAddress("R_L"       , &R_L);
        ntuple_data->SetBranchAddress("efficiency", &efficiency);
        ntuple_data->SetBranchAddress("purity"    , &purity);
        ntuple_data->SetBranchAddress("jet_pt"    , &jet_pt);
        ntuple_data->SetBranchAddress("weight"    , &weight);
        ntuple_data->SetBranchAddress("weight_pt" , &weight_pt);
        
        // Unfold the purity corrected npairs
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

        TH2D* hunfolded_ratio     = new TH2D("hunfolded_ratio"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hpuritycorrected    = new TH2D("hpuritycorrected" ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hpuritycorrected2   = new TH2D("hpuritycorrected2","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        ntuple_data->Project("hpuritycorrected" ,"jet_pt:R_L","purity*weight");
        ntuple_data->Project("hpuritycorrected2","jet_pt:R_L","purity*weight");
        
        RooUnfoldBayes unfold(response, hpuritycorrected, niter);
        TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
        hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected2,1,1);
        hunfolded_ratio->Smooth();

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        gStyle->SetOptStat("");
        gStyle->SetPaintTextFormat("4.2f");
        hunfolded_ratio->Draw("col text");
        hunfolded_ratio->SetTitle("Purity Correction;R_{L};p_{T,jet}(GeV)");
        hunfolded_ratio->Smooth();
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/unfolding_effect_smooth.pdf");

        
        // Unfold the purity corrected jets
        TNtuple* ntuple_jet_unf = (TNtuple*) f->Get(name_ntuple_mcreco_jet.c_str());

        float unf_jet_pt_reco, unf_jet_pt_truth;
        ntuple_jet_unf->SetBranchAddress("jet_pt", &unf_jet_pt_reco);
        ntuple_jet_unf->SetBranchAddress("jet_pt_truth", &unf_jet_pt_truth);
        
        TH1D* hpurcorr_jet = new TH1D("hpurcorr_jet","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH1D* hmeas_jet    = new TH1D("hmeas_jet"   ,"",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH1D* htrue_jet    = new TH1D("htrue_jet"   ,"",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        TH2D* hresponse_jet = new TH2D("hresponse_jet","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        for (int evt = 0 ; evt < ntuple_jet_unf->GetEntries() ; evt++) {
                ntuple_jet_unf->GetEntry(evt);

                if (unf_jet_pt_truth != -999) 
                        hresponse_jet->Fill(unf_jet_pt_reco, unf_jet_pt_truth);
        }

        RooUnfoldResponse* response_jet = new RooUnfoldResponse(hmeas_jet, htrue_jet, hresponse_jet, "response_jet");

        TH1D* hunfolded_ratio_jet   = new TH1D("hunfolded_ratio_jet"  ,"",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH1D* hpuritycorrected_jet  = new TH1D("hpuritycorrected_jet" ,"",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH1D* hpuritycorrected2_jet = new TH1D("hpuritycorrected2_jet","",nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        
        ntuple_jet->Project("hpuritycorrected_jet" ,"jet_pt","jet_purity");
        ntuple_jet->Project("hpuritycorrected2_jet","jet_pt","jet_purity");
        
        RooUnfoldBayes unfold_jet(response_jet, hpuritycorrected_jet, 3);
        TH1D* hunfolded_bayes_jet = (TH1D*) unfold_jet.Hunfold();
        hunfolded_ratio_jet->Divide(hunfolded_bayes_jet,hpuritycorrected2_jet,1,1);
        hunfolded_ratio_jet->Smooth();

        c->Draw();

        gStyle->SetOptStat("");
        gStyle->SetPaintTextFormat("4.2f");
        hunfolded_ratio_jet->Draw();
        hunfolded_ratio_jet->SetTitle("Purity Correction;R_{L};p_{T,jet}(GeV)");
        hunfolded_ratio_jet->Smooth();
        // gPad->SetLogx(1);
        // gPad->SetLogy(1);
        c->Print("./plots/unfolding_effect_jetpt_smooth.pdf");


        TH1F* hcorr_jet[nbin_jet_pt];
        TH1F* hcorr_eec[nbin_jet_pt]; 
        
        // Fill the histograms
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_jet[bin] = new TH1F(Form("hcorr_jet%i" ,bin),"", 1,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                hcorr_eec[bin] = new TH1F(Form("hcorr_eec%i",bin),"", nbin_rl_nominal, rl_nominal_binning);
                
                set_histogram_style(hcorr_eec[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size);
                
                ntuple_jet->Project(Form("hcorr_jet%i" ,bin), "jet_pt", jet_full_corr[bin]);

                if (apply_jet_unfolding)
                        hcorr_jet[bin]->Scale(hunfolded_ratio_jet->GetBinContent(bin + 2)); // hunfolded_ratio_jet has two more bins
        
                for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++) {
                        ntuple_data->GetEntry(entry);

                        if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) 
                                continue;
                        
                        if (efficiency <= 0 || purity <= 0) 
                                continue;
                        
                        if (efficiency > 1 || purity > 1) 
                                continue;
                        
                        double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt));
                        
                        if (apply_unfolding)
                                hcorr_eec[bin]->Fill(R_L,unfolding_weight*weight*weight_pt*purity/efficiency);
                        else
                                hcorr_eec[bin]->Fill(R_L,weight*weight_pt*purity/efficiency);
                }

                hcorr_eec[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
        }

        // Simulations Section
        TFile* fmc   = new TFile((output_folder+namef_ntuple_mc_eec).c_str());
        
        TNtuple* ntuple_mc     = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
        TNtuple* ntuple_mc_jet = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());
        
        float R_L_mc, jet_pt_mc, weight_pt_mc;
        ntuple_mc->SetBranchAddress("R_L",&R_L_mc);
        ntuple_mc->SetBranchAddress("jet_pt",&jet_pt_mc);
        ntuple_mc->SetBranchAddress("weight_pt",&weight_pt_mc);
        
        TH1F* hmc[nbin_jet_pt]; 
        TH1F* hmc_jet[nbin_jet_pt];

        c = new TCanvas("c", "", 1920, 600);
        c->Draw();
        c->Divide(3,1);

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        THStack* s_data[3];
        TLegend* l_data[3];

        for(int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hmc[bin]     = new TH1F(Form("hmc[%i]", bin)    ,"",nbin_rl_nominal,rl_nominal_binning);
                hmc_jet[bin] = new TH1F(Form("hmc_jet[%i]", bin),"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]);
                
                set_histogram_style(hmc[bin], corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin], std_marker_size);
                
                // Fill and normalize MC        
                for(int entry = 0 ; entry < ntuple_mc->GetEntries() ; entry++) {
                        ntuple_mc->GetEntry(entry);

                        if(jet_pt_mc < jet_pt_binning[bin] || jet_pt_mc > jet_pt_binning[bin+1]) 
                                continue;
                        
                        hmc[bin]->Fill(R_L_mc, weight_pt_mc);
                }
                
                ntuple_mc_jet->Project(Form("hmc_jet[%i]" ,bin),"jet_pt",pair_jet_pt_cut[bin]);
                hmc[bin]->Scale(1./hmc_jet[bin]->Integral(),"width");
        }

        // Draw the log binning histos
        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                c->cd(bin+1);

                s_data[bin] = new THStack();
                l_data[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,1.-gPad->GetTopMargin(),1-gPad->GetRightMargin()-0.01,1.-gPad->GetTopMargin()-0.21);

                s_data[bin]->Add(hmc[bin],"E");
                s_data[bin]->Add(hcorr_eec[bin],"E");
                
                s_data[bin]->Draw("NOSTACK");
                
                s_data[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};N_{pairs}", jet_pt_binning[bin], jet_pt_binning[bin+1]));

                s_data[bin]->SetMaximum(1.5);
                s_data[bin]->SetMinimum(0.);
                
                l_data[bin]->AddEntry(hmc[bin]      , "mc"  , "p");
                l_data[bin]->AddEntry(hcorr_eec[bin], "data", "p");

                l_data[bin]->Draw("SAME");    

                gPad->SetLogx(1);
                // gPad->SetLogy(1);
        }

        c->Print(Form("./plots/data2mc_eec_paircorr_fullsim_corr-smoothed_unf-%s-smoothed_unfjet-%s-smoothed_altcorrbinning.pdf",(apply_unfolding)?"yes":"no",(apply_jet_unfolding)?"yes":"no"));
}