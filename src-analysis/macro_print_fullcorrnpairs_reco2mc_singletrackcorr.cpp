#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorrnpairs_reco2mc_singletrackcorr(int niter = 10, bool do_print = true)
{
        gStyle->SetPadTopMargin(0.08);

        TFile* fcorr = new TFile((output_folder + "ntuple_corrmcrecoe2c.root").c_str()); 
        if (fcorr->IsZombie()) 
                return;
        
        TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        float R_L, efficiency, purity, jet_pt, weight;
        
        ntuple_data->SetBranchAddress("R_L"       , &R_L);
        ntuple_data->SetBranchAddress("efficiency", &efficiency);
        ntuple_data->SetBranchAddress("purity"    , &purity);
        ntuple_data->SetBranchAddress("jet_pt"    , &jet_pt);
        ntuple_data->SetBranchAddress("weight"    , &weight);

        // Unfold the purity corrected data
        TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections).c_str());
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
        
        ntuple_data->Project("hpuritycorrected" ,"jet_pt:R_L","purity");
        ntuple_data->Project("hpuritycorrected2","jet_pt:R_L","purity");
        
        RooUnfoldBayes unfold(response, hpuritycorrected, niter);
        TH2D* hunfolded_bayes   = (TH2D*) unfold.Hunfold();
        hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected2,1,1);

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        gStyle->SetOptStat("");
        gStyle->SetPaintTextFormat("4.2f");
        hunfolded_ratio->Draw("col text");
        hunfolded_ratio->SetTitle("Purity Correction;R_{L};p_{T,jet}(GeV)");
        hunfolded_ratio->Smooth();
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        c->Print("./plots/unfolding_effect.pdf");

        
        TH1F* hcorr_jet[nbin_jet_pt];
        TH1F* hcorr_npairs[nbin_jet_pt]; 
        
        // Fill the histograms
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_jet[bin] = new TH1F(Form("hcorr_jet%i" ,bin),"", 1,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                
                hcorr_npairs[bin] = new TH1F(Form("hcorr_npairs%i",bin),"", nbin_rl_nominal, rl_nominal_binning);
                
                set_histogram_style(hcorr_npairs[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                
                ntuple_jet->Project(Form("hcorr_jet%i" ,bin), "jet_pt", jet_full_corr[bin]);
        
                for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++) {
                        ntuple_data->GetEntry(entry);

                        if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) 
                                continue;
                        
                        if (efficiency <= 0 || purity <= 0) 
                                continue;
                        
                        if (efficiency > 1 || purity > 1) 
                                continue;

                        double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt));
                        if (unfolding_weight <= 0) 
                                unfolding_weight = 1;
                        
                        hcorr_npairs[bin]->Fill(R_L,weight*unfolding_weight*purity/efficiency);
                }

                hcorr_npairs[bin]->Scale(1./hcorr_jet[bin]->Integral());
        }

        // Simulations Section
        TFile* fmc   = new TFile((output_folder+namef_ntuple_mc_e2c).c_str());
        
        TNtuple* ntuple_mc     = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
        TNtuple* ntuple_mc_jet = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());
        
        float R_L_mc, jet_pt_mc;
        ntuple_mc->SetBranchAddress("R_L",&R_L_mc);
        ntuple_mc->SetBranchAddress("jet_pt",&jet_pt_mc);
        
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
                hmc[bin]     = new TH1F(Form("hmc[%i]", bin),"",nbin_rl_nominal,rl_nominal_binning);
                hmc_jet[bin] = new TH1F(Form("hmc_jet[%i]", bin),"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]);
                
                set_histogram_style(hmc[bin], corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin], std_marker_size);
                
                // Fill and normalize MC        
                for(int entry = 0 ; entry < ntuple_mc->GetEntries() ; entry++) {
                        ntuple_mc->GetEntry(entry);

                        if(jet_pt_mc < jet_pt_binning[bin] || jet_pt_mc > jet_pt_binning[bin+1]) 
                                continue;
                        
                        hmc[bin]->Fill(R_L_mc);
                }
                
                ntuple_mc_jet->Project(Form("hmc_jet[%i]" ,bin),"jet_pt",pair_jet_pt_cut[bin]);
                hmc[bin]->Scale(1./hmc_jet[bin]->Integral());
        }

        // Draw the log binning histos
        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                c->cd(bin+1);

                s_data[bin] = new THStack();
                l_data[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,gPad->GetBottomMargin()+0.01,1-gPad->GetRightMargin()-0.01,gPad->GetBottomMargin()+0.16);

                s_data[bin]->Add(hmc[bin],"E");
                s_data[bin]->Add(hcorr_npairs[bin],"E");
                
                s_data[bin]->Draw("NOSTACK");
                
                s_data[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};N_{pairs}", jet_pt_binning[bin], jet_pt_binning[bin+1]));
                
                l_data[bin]->AddEntry(hmc[bin]         , "mc"  , "p");
                l_data[bin]->AddEntry(hcorr_npairs[bin], "data", "p");

                l_data[bin]->Draw("SAME");    

                gPad->SetLogx(1);
        }

        c->Print("./plots/reco2mc_npairs_multiplicity_fullconditions_wredundantjetid_eventweightapplied_pseudatajetmatchedjet_unfolded.pdf");
}