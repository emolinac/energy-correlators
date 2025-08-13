#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorrnpairs_data2mc(int niter = 15, bool do_print = true)
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
                        
                        hcorr_npairs[bin]->Fill(R_L,weight*purity/efficiency);
                }

                hcorr_npairs[bin]->Scale(1./hcorr_jet[bin]->Integral());
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

        TCanvas* c = new TCanvas("c", "", 1920, 600);
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

                // s_data[bin]->SetMaximum(6.);
                s_data[bin]->SetMinimum(0.);
                
                l_data[bin]->AddEntry(hmc[bin]      , "mc"  , "p");
                l_data[bin]->AddEntry(hcorr_npairs[bin], "data", "p");

                l_data[bin]->Draw("SAME");    

                gPad->SetLogx(1);
                // gPad->SetLogy(1);
        }

        c->Print("./plots/data2mc_npairs_paircorr_fullsim_smooth.pdf");
}