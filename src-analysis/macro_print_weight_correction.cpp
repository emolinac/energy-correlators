#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_weight_correction(bool do_print = true)
{
        gStyle->SetPadTopMargin(0.08);

        // Simulations Section
        TFile* fmc = new TFile((output_folder+namef_ntuple_eec_paircorrections).c_str());
        
        TNtuple* ntuple_mc = (TNtuple*) fmc->Get((name_ntuple_correction_reco).c_str());
        
        float R_L_mc, jet_pt_mc, weight_pt_mc, weight_pt_reco;
        ntuple_mc->SetBranchAddress("R_L",&R_L_mc);
        ntuple_mc->SetBranchAddress("jet_pt",&jet_pt_mc);
        ntuple_mc->SetBranchAddress("weight_pt",&weight_pt_mc);
        ntuple_mc->SetBranchAddress("weight_pt_truth",&weight_pt_reco);
        
        TH1F* hmc_truthweights[nbin_jet_pt]; 
        TH1F* hmc[nbin_jet_pt]; 
        TH1F* hratio[nbin_jet_pt]; 
        
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
                hratio[bin]               = new TH1F(Form("hratio[%i]" ,bin)      ,"",nbin_rl_nominal,rl_nominal_binning);
                hmc_truthweights[bin] = new TH1F(Form("hmc_truthweights[%i]" ,bin),"",nbin_rl_nominal,rl_nominal_binning);
                hmc[bin]              = new TH1F(Form("hmc[%i]" ,bin)             ,"",nbin_rl_nominal,rl_nominal_binning);
                
                set_histogram_style(hmc_truthweights[bin], corr_marker_color_jet_pt[1], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
                set_histogram_style(hmc[bin]             , corr_marker_color_jet_pt[2], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
                set_histogram_style(hratio[bin]          , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin] , std_marker_size);
                
                // Fill and normalize MCReco
                for(int entry = 0 ; entry < ntuple_mc->GetEntries() ; entry++) {
                        ntuple_mc->GetEntry(entry);

                        if(jet_pt_mc<jet_pt_binning[bin]||jet_pt_mc>jet_pt_binning[bin+1]) 
                                continue;

                        if(weight_pt_reco == -999)
                                continue;
                        
                        hmc[bin]->Fill(R_L_mc,weight_pt_mc);
                        hmc_truthweights[bin]->Fill(R_L_mc,weight_pt_reco);
                }
        }

        // Draw the log binning histos
        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                c->cd(bin+1);
                s_data[bin] = new THStack();
                l_data[bin] = new TLegend(1-gPad->GetRightMargin()-0.31,1-gPad->GetTopMargin()-0.15,1-gPad->GetRightMargin()-0.01,1-gPad->GetTopMargin()-0.01);

                // s_data[bin]->Add(hmc_truthweights[bin],"E");
                // s_data[bin]->Add(hmc[bin],"E");

                hratio[bin]->Divide(hmc_truthweights[bin],hmc[bin],1,1);

                s_data[bin]->Add(hratio[bin],"E");
                s_data[bin]->Draw("NOSTACK");
                s_data[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};#Sigma_{EEC}(R_{L})",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s_data[bin]->SetMaximum(1.5);
                s_data[bin]->SetMinimum(0.5);

                l_data[bin]->AddEntry(hratio[bin],"mc2mc-w-recoweight","p");

                gPad->SetLogx(1);
                
                l_data[bin]->Draw("SAME");    
        }

        c->Print("./plots/weightcorrection.pdf");

        TFile* fweightcorr = new TFile((output_folder+namef_histos_weight_corr).c_str(),"RECREATE");
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++)
                hratio[bin]->Write(Form("h_weightcorr_%i",bin));
}
