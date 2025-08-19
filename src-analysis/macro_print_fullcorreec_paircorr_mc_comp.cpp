#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorreec_paircorr_mc_comp(int niter = 4, bool do_print = true)
{
        gStyle->SetPadTopMargin(0.08);

        TFile* fcorr = new TFile((output_folder + namef_ntuple_eec_paircorr).c_str()); 
        if (fcorr->IsZombie()) 
                return;
        
        TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        // Set the branches of data
        float R_L, jet_pt, weight_pt;
        float efficiency, purity, efficiency_factorized, purity_factorized;
        ntuple_data->SetBranchAddress("R_L", &R_L);
        ntuple_data->SetBranchAddress("jet_pt", &jet_pt);
        ntuple_data->SetBranchAddress("weight_pt", &weight_pt);
        ntuple_data->SetBranchAddress("efficiency", &efficiency);
        ntuple_data->SetBranchAddress("purity", &purity);
        
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
        
        ntuple_data->Project("hpuritycorrected"    ,"jet_pt:R_L");
        ntuple_data->Project("hpuritycorrected_ref","jet_pt:R_L");
        
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
        
        gStyle->SetPaintTextFormat("4.2f");
        
        // Fill the NOIMNAL histograms
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)         ,"", 1  ,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                hcorr_jet_centroid[bin] = new TH1F(Form("hcorr_jet_centroid%i" ,bin),"", 200,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 

                hcorr_eec[bin]          = new TH1F(Form("hcorr_eec%i",bin)         ,"", nbin_rl_nominal,rl_nominal_binning );
                hcorr_eec_syst[bin]     = new TH1F(Form("hcorr_eec_syst%i",bin)    ,"", nbin_rl_nominal,rl_nominal_binning );
                hcorr_tau[bin]          = new TH1F(Form("hcorr_tau%i",bin)         ,"", nbin_rl_nominal,tau_nominal_binning);
                hcorr_tau_syst[bin]     = new TH1F(Form("hcorr_tau_syst%i",bin)    ,"", nbin_rl_nominal,tau_nominal_binning );
        
                set_histogram_style(hcorr_eec[bin]     , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_eec_syst[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size  );
                set_histogram_style(hcorr_tau[bin]     , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hcorr_tau_syst[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size  );
        
                hcorr_eec[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_tau[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_tau_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                
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

                        hcorr_eec[bin]->Fill(R_L,unfolding_weight*weight_pt/efficiency);
                        hcorr_eec_syst[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_tau[bin]->Fill(R_L*jet_pt_centroid,purity*unfolding_weight*weight_pt/efficiency);
                        hcorr_tau_syst[bin]->Fill(R_L*jet_pt_centroid,purity*unfolding_weight*weight_pt/efficiency);
                }

                hcorr_eec[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_eec_syst[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                
                hcorr_tau[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                hcorr_tau_syst[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
                
        }

        // // Include the systematics in the whole deal
        // const int nsyst = sizeof(available_systematics)/sizeof(available_systematics[0]);
        // TFile* fsyst[nsyst];
        // TH1F* hdev[nbin_jet_pt];
        // TH1F* hdev_tau[nbin_jet_pt];

        // std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        // std::cout<<"\\hline"<<std::endl;
        // for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
        //         fsyst[syst_index] = new TFile((output_folder + devfromnom_namef[available_systematics[syst_index]]).c_str());
                
        //         if (fsyst[syst_index]->IsZombie()) 
        //                 continue;
                
        //         std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                
        //         for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
        //                 hdev[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_deviations_eec%i",bin));

        //                 set_histo_with_systematics(hdev[bin], hcorr_eec[bin], hcorr_eec_syst[bin], systematic_errtype[available_systematics[syst_index]]);

        //                 if (bin!=nbin_jet_pt-1) 
        //                         std::cout<<" & ";
        //                 else
        //                         std::cout<<" \\\\ ";
        //         }

        //         std::cout<<std::endl;

        //         delete fsyst[syst_index];
        // }

        // std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        // std::cout<<"\\hline"<<std::endl;
        // for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
        //         fsyst[syst_index] = new TFile((output_folder + devfromnom_namef[available_systematics[syst_index]]).c_str());
                
        //         if (fsyst[syst_index]->IsZombie()) 
        //                 continue;

        //         std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
        //         for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
        //                 hdev_tau[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_deviations_tau%i",bin));

        //                 set_histo_with_systematics(hdev_tau[bin], hcorr_tau[bin], hcorr_tau_syst[bin], systematic_errtype[available_systematics[syst_index]]);

        //                 if (bin!=nbin_jet_pt-1) 
        //                         std::cout<<" & ";
        //                 else
        //                         std::cout<<" \\\\ ";
        //         }

        //         std::cout<<std::endl;

        //         delete fsyst[syst_index];
        // }

        // Simulations Section
        TFile* fmc   = new TFile((output_folder+namef_ntuple_mc_eec).c_str());
        
        TNtuple* ntuple_mc         = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
        TNtuple* ntuple_mc_jet     = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());
        TNtuple* ntuple_mcreco     = (TNtuple*) fmc->Get((name_ntuple_mcreco).c_str());
        TNtuple* ntuple_mcreco_jet = (TNtuple*) fmc->Get((name_ntuple_mcreco_jet).c_str());

        // Set the branches of data
        float R_L_mc, jet_pt_mc, weight_pt_mc;
        ntuple_mc->SetBranchAddress("R_L",&R_L_mc);
        ntuple_mc->SetBranchAddress("jet_pt",&jet_pt_mc);
        ntuple_mc->SetBranchAddress("weight_pt",&weight_pt_mc);
        
        // Set the branches of data
        float R_L_mcreco, jet_pt_mcreco, weight_pt_mcreco;
        ntuple_mcreco->SetBranchAddress("R_L",&R_L_mcreco);
        ntuple_mcreco->SetBranchAddress("jet_pt",&jet_pt_mcreco);
        ntuple_mcreco->SetBranchAddress("weight_pt",&weight_pt_mcreco);
        
        TH1F* hmc[nbin_jet_pt]; 
        TH1F* hmcreco[nbin_jet_pt]; 
        TH1F* hmc_jet[nbin_jet_pt]; 
        TH1F* hmcreco_jet[nbin_jet_pt]; 
        
        TCanvas* c = new TCanvas("c","",1800,600);
        c->Draw();
        c->Divide(3,1);
        
        TLatex* tex = new TLatex();
        tex->SetTextColorAlpha(16,0.3);
        tex->SetTextSize(0.1991525);
        tex->SetTextAngle(26.15998);
        tex->SetLineWidth(2);

        THStack* s_data[3];
        TLegend* l_data[3];

        for(int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hmc[bin]         = new TH1F(Form("hmc[%i]" ,bin)      ,"",nbin_rl_nominal,rl_nominal_binning);
                hmcreco[bin]     = new TH1F(Form("hmcreco[%i]" ,bin)  ,"",nbin_rl_nominal,rl_nominal_binning);
                
                hmc_jet[bin]     = new TH1F(Form("hmc_jet[%i]" ,bin)     ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]);
                hmcreco_jet[bin] = new TH1F(Form("hmcreco_jet[%i]" ,bin) ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]);
                
                set_histogram_style(hmc[bin]        , corr_marker_color_jet_pt[1], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
                set_histogram_style(hmcreco[bin]    , corr_marker_color_jet_pt[2], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
                
                // Fill and normalize MC        
                for(int entry = 0 ; entry < ntuple_mc->GetEntries() ; entry++) {
                        ntuple_mc->GetEntry(entry);

                        if(jet_pt_mc<jet_pt_binning[bin]||jet_pt_mc>jet_pt_binning[bin+1]) 
                                continue;
                        
                        hmc[bin]->Fill(R_L_mc,weight_pt_mc);
                }
                
                ntuple_mc_jet->Project(Form("hmc_jet[%i]" ,bin),"jet_pt",pair_jet_pt_cut[bin]);
                hmc[bin]->Scale(1./hmc_jet[bin]->Integral(),"width");
                
                // Fill and normalize MCReco
                for(int entry = 0 ; entry < ntuple_mcreco->GetEntries() ; entry++) {
                        ntuple_mcreco->GetEntry(entry);

                        if(jet_pt_mcreco<jet_pt_binning[bin]||jet_pt_mcreco>jet_pt_binning[bin+1]) 
                                continue;
                        
                        hmcreco[bin]->Fill(R_L_mcreco,weight_pt_mcreco);
                }

                ntuple_mcreco_jet->Project(Form("hmcreco_jet[%i]" ,bin),"jet_pt",pair_jet_pt_cut[bin]);
                hmcreco[bin]->Scale(1./hmcreco_jet[bin]->Integral(),"width");
        }

        // Draw the log binning histos
        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                c->cd(bin+1);
                s_data[bin] = new THStack();
                l_data[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,1-gPad->GetTopMargin()-0.15,1-gPad->GetRightMargin()-0.01,1-gPad->GetTopMargin()-0.01);

                s_data[bin]->Add(hmc[bin],"E");
                s_data[bin]->Add(hmcreco[bin],"E");
                s_data[bin]->Add(hcorr_eec[bin],"E");
                s_data[bin]->Draw("NOSTACK");
                s_data[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};#Sigma_{EEC}(R_{L})",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s_data[bin]->SetMaximum(1.);
                s_data[bin]->SetMinimum(40E-03);

                l_data[bin]->AddEntry(hmc[bin]      ,"mc"    ,"p");
                l_data[bin]->AddEntry(hmcreco[bin]  ,"mcreco","p");
                l_data[bin]->AddEntry(hcorr_eec[bin],"data"  ,"p");
                gPad->SetLogx(1);
                l_data[bin]->Draw("SAME");    
        }

        c->Print(Form("./plots/data2mc_correec_unf-niter%i.pdf",niter));
}