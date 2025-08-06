#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorrhadronmultiplicity_reco2mc(int niter = 15, bool do_print = true)
{
        gStyle->SetPadTopMargin(0.08);

        TFile* fcorr = new TFile((output_folder + "ntuple_mcreco_hadronmultiplicity.root").c_str()); 
        if (fcorr->IsZombie()) 
                return;
        
        TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
        TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
        
        float h1_eta, efficiency, purity, jet_pt, weight;
        
        ntuple_data->SetBranchAddress("h1_eta"    , &h1_eta);
        ntuple_data->SetBranchAddress("efficiency", &efficiency);
        ntuple_data->SetBranchAddress("purity"    , &purity);
        ntuple_data->SetBranchAddress("jet_pt"    , &jet_pt);
        ntuple_data->SetBranchAddress("weight"    , &weight);
        
        TH1F* hcorr_jet[nbin_jet_pt];
        TH1F* hcorr_hadrons[nbin_jet_pt]; 
        
        // Fill the histograms
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_jet[bin] = new TH1F(Form("hcorr_jet%i" ,bin),"", 1,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
                
                hcorr_hadrons[bin] = new TH1F(Form("hcorr_hadrons%i",bin),"", sl_eta_nbins,sl_eta_binning);
                
                set_histogram_style(hcorr_hadrons[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                
                ntuple_jet->Project(Form("hcorr_jet%i" ,bin), "jet_pt", jet_full_corr[bin]);
        
                for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++) {
                        ntuple_data->GetEntry(entry);

                        if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) 
                                continue;
                        
                        if (efficiency <= 0 || purity <= 0) 
                                continue;
                        
                        hcorr_hadrons[bin]->Fill(h1_eta,weight*purity/efficiency);
                }

                hcorr_hadrons[bin]->Scale(1./hcorr_jet[bin]->Integral());
        }

        // Simulations Section
        TFile* fmc   = new TFile((output_folder+namef_ntuple_hadron).c_str());
        
        TNtuple* ntuple_mc         = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
        TNtuple* ntuple_mc_jet     = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());
        TNtuple* ntuple_mcreco     = (TNtuple*) fmc->Get((name_ntuple_mcreco).c_str());
        TNtuple* ntuple_mcreco_jet = (TNtuple*) fmc->Get((name_ntuple_mcreco_jet).c_str());

        float h_eta_mc, jet_pt_mc;
        ntuple_mc->SetBranchAddress("h_eta",&h_eta_mc);
        ntuple_mc->SetBranchAddress("jet_pt",&jet_pt_mc);
        
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
                hmc[bin]     = new TH1F(Form("hmc[%i]", bin),"",sl_eta_nbins,sl_eta_binning);
                hmc_jet[bin] = new TH1F(Form("hmc_jet[%i]", bin),"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]);
                
                set_histogram_style(hmc[bin], corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin], std_marker_size);
                
                // Fill and normalize MC        
                for(int entry = 0 ; entry < ntuple_mc->GetEntries() ; entry++) {
                        ntuple_mc->GetEntry(entry);

                        if(jet_pt_mc < jet_pt_binning[bin] || jet_pt_mc > jet_pt_binning[bin+1]) 
                                continue;
                        
                        hmc[bin]->Fill(h_eta_mc);
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
                s_data[bin]->Add(hcorr_hadrons[bin],"E");
                
                s_data[bin]->Draw("NOSTACK");
                
                s_data[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;P(GeV);Hadron Multiplicity", jet_pt_binning[bin], jet_pt_binning[bin+1]));
                
                l_data[bin]->AddEntry(hmc[bin]          , "mc"  , "p");
                l_data[bin]->AddEntry(hcorr_hadrons[bin], "data", "p");

                // l_data[bin]->Draw("SAME");    
        }

        c->Print("./plots/reco2mc_hadronmultiplicity_fullconditions_wredundantjetid.pdf");
}