#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorreec_paircorr_mc_comp(int niter = 4, bool do_print = true, bool apply_weight_corr = false)
{
        gStyle->SetPadTopMargin(0.08);

        TString variation = "3dunf";

        TFile* fdata = new TFile((output_folder + namef_histos_paircorr_eec_3dunf).c_str());

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

                htruth_eec[bin] = new TH1F(Form("htruth_eec%i",bin),"",nbin_rl_nominal,rl_nominal_binning);
                ntuple_mc->Project(Form("htruth_eec%i",bin),"R_L",eec_jet_pt_cut[bin]);

                htruth_eec[bin]->Scale(1./htruth_jet[bin]->Integral(),"width");
                // htruth_eec[bin]->Scale(2./htruth_jet[bin]->Integral(),"width"); // For 4d and 5d 

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
                // s_data[bin]->SetMaximum(1.7);
                // s_data[bin]->SetMinimum(40E-03);

                l_data[bin]->AddEntry(htruth_eec[bin],"mc"    ,"p");
                l_data[bin]->AddEntry(hcorr_eec[bin],"data"  ,"p");
                gPad->SetLogx(1);
                l_data[bin]->Draw("SAME");    
        }

        c->Print(Form("./plots/data2mc_correec_unf-niter%i_"+variation+".pdf",niter));

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

                l_data[bin]->AddEntry(htruth_eec[bin],"data/mc","p");
                gPad->SetLogx(1);
                l_data[bin]->Draw("SAME");    
        }

        c_ratio->Print(Form("./plots/data2mcratio_correec_unf-niter%i_"+variation+".pdf",niter));
}