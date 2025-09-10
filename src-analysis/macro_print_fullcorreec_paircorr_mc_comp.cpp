#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorreec_paircorr_mc_comp()
{
        gStyle->SetPadTopMargin(0.08);

        TString variation = "3dunf-7niter-ct";

        TFile* fdata = new TFile((output_folder + "histos_histopaircorr_eec_3dunf-niter7_ct.root").c_str());

        TH1F* hcorr_eec[nbin_jet_pt];
        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                hcorr_eec[i] = (TH1F*) fdata->Get(Form("hcorr_eec%i",i));

                set_histogram_style(hcorr_eec[i], corr_marker_color_jet_pt[i], std_line_width, corr_marker_style_jet_pt[i] , std_marker_size);
        }

        // Simulations Section
        TFile* fmc   = new TFile((output_folder+namef_ntuple_mc_eec).c_str());
        
        TNtuple* ntuple_mc     = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
        TNtuple* ntuple_mc_jet = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());
        
        TH1F* htruth_eec[nbin_jet_pt]; 
        TH1F* htruth_jet[nbin_jet_pt]; 
        
        for(int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                htruth_jet[bin] = new TH1F(Form("htruth_jet%i" ,bin),"",200,jet_pt_binning[bin],jet_pt_binning[bin + 1]);
                ntuple_mc_jet->Project(Form("htruth_jet%i" ,bin),"jet_pt");

                htruth_eec[bin] = new TH1F(Form("htruth_eec%i",bin),"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                ntuple_mc->Project(Form("htruth_eec%i",bin),"R_L",eec_jet_pt_cut[bin]);

                for (int i = 1 ; i <= htruth_eec[bin]->GetNbinsX(); i++) {
                        if (i == 1 || i == htruth_eec[bin]->GetNbinsX()) {
                                htruth_eec[bin]->SetBinContent(i, 0);
                                htruth_eec[bin]->SetBinError(i, 0);        

                                continue;
                        }
                }

                htruth_eec[bin]->Scale(1./htruth_jet[bin]->Integral(),"width");

                set_histogram_style(htruth_eec[bin], corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
        }

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

        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                c->cd(bin+1);
                s_data[bin] = new THStack();
                l_data[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,1-gPad->GetTopMargin()-0.15,1-gPad->GetRightMargin()-0.01,1-gPad->GetTopMargin()-0.01);

                s_data[bin]->Add(htruth_eec[bin],"E");
                s_data[bin]->Add(hcorr_eec[bin],"E");
                s_data[bin]->Draw("NOSTACK");
                s_data[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};#Sigma_{EEC}(R_{L})",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s_data[bin]->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                
                s_data[bin]->SetMaximum(1.3);
                // s_data[bin]->SetMinimum(0.5);
                
                l_data[bin]->AddEntry(htruth_eec[bin],"mc"    ,"p");
                l_data[bin]->AddEntry(hcorr_eec[bin],"data"  ,"p");
                gPad->SetLogx(1);
                l_data[bin]->Draw("SAME");    
        }

        c->Print("./plots/data2mc_correec_"+variation+".pdf");

        TCanvas* c_ratio = new TCanvas("c_ratio","",1800,400);
        c_ratio->Draw();
        c_ratio->Divide(3,1);
        
        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                c_ratio->cd(bin+1);
                s_data[bin] = new THStack();
                l_data[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,1-gPad->GetTopMargin()-0.15,1-gPad->GetRightMargin()-0.01,1-gPad->GetTopMargin()-0.01);

                hcorr_eec[bin]->Divide(htruth_eec[bin]);
                
                s_data[bin]->Add(hcorr_eec[bin],"E");
                s_data[bin]->Draw("NOSTACK");
                s_data[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};#Sigma_{EEC}(R_{L})",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s_data[bin]->SetMaximum(1.5);
                s_data[bin]->SetMinimum(0.5);
                s_data[bin]->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);

                l_data[bin]->AddEntry(htruth_eec[bin],"data/mc","p");
                gPad->SetLogx(1);
                l_data[bin]->Draw("SAME");    
        }

        c_ratio->Print("./plots/data2mcratio_correec_"+variation+".pdf");
}