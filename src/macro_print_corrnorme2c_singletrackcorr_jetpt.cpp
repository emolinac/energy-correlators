#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_corrnorme2c_singletrackcorr_jetpt()
{
    // Open the necessary files
    TFile* fcorr = new TFile((output_folder+namef_ntuple_e2c_corr).c_str()); 
    if(fcorr->IsZombie()) return;
    
    TFile* ftest = new TFile("./test.root","RECREATE");
    gROOT->cd();

    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
    
    TH1F* hcorr_data[Nbin_jet_pt]; 
    TH1F* hall_data[Nbin_jet_pt];  
    TH1F* hcorr_jet[Nbin_jet_pt]; 
    TH1F* hall_jet[Nbin_jet_pt];
    TH1F* hallref_jet[Nbin_jet_pt];
    TH1F* hcorrref_jet[Nbin_jet_pt];

    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend();
    
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        hallref_jet[bin]  = new TH1F(Form("hallref_jet[%i]" ,bin) ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]); 
        hcorrref_jet[bin] = new TH1F(Form("hcorrref_jet[%i]" ,bin) ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]); 

        hall_jet[bin]     = new TH1F(Form("hall_jet[%i]" ,bin) ,"",Nbin_R_L,R_L_min,R_L_max);
        hcorr_jet[bin]    = new TH1F(Form("hcorr_jet[%i]" ,bin),"",Nbin_R_L,R_L_min,R_L_max);
        hcorr_data[bin]   = new TH1F(Form("hcorr_data[%i]",bin),"",Nbin_R_L,R_L_min,R_L_max);
        hall_data[bin]    = new TH1F(Form("hall_data[%i]" ,bin),"",Nbin_R_L,R_L_min,R_L_max);
        
        set_histogram_style(hcorr_data[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        set_histogram_style(hall_data[bin] , std_marker_color_jet_pt[bin] , std_line_width, std_marker_style_jet_pt[bin] , std_marker_size+1);

        // Project into the histograms
        ntuple_data->Project(Form("hcorr_data[%i]",bin),"R_L",e2c_jetpt_full_corr_singletrack[bin]);
        ntuple_data->Project(Form("hall_data[%i]" ,bin),"R_L",e2c_jetpt_cut[bin]);
        ntuple_jet->Project(Form("hcorrref_jet[%i]" ,bin),"jet_pt",jet_full_corr[bin]);
        ntuple_jet->Project(Form("hallref_jet[%i]" ,bin),"jet_pt",pair_jetpt_cut[bin]);

        for(int rl_bin = 0 ; rl_bin < Nbin_R_L ; rl_bin++)
        {
            hall_jet[bin]->SetBinContent(rl_bin+1,  hallref_jet[bin]->GetBinContent(1)); 
            hcorr_jet[bin]->SetBinContent(rl_bin+1, hcorrref_jet[bin]->GetBinContent(1));
            hall_jet[bin]->SetBinError(rl_bin+1,  hallref_jet[bin]->GetBinError(1)); 
            hcorr_jet[bin]->SetBinError(rl_bin+1, hcorrref_jet[bin]->GetBinError(1));
        }

        hcorr_data[bin]->Divide(hcorr_jet[bin]);
        hall_data[bin]->Divide(hall_jet[bin]);

        ftest->cd();
        hall_jet[bin]->Write();
        hcorr_jet[bin]->Write();
        hcorr_data[bin]->Write();
        hall_data[bin]->Write();
        gROOT->cd();

        s_data->Add(hcorr_data[bin],"EX0");
        // s_data->Add(hall_data[bin]);
        l_data->AddEntry(hcorr_data[bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]),"lf");
        // l_data->AddEntry(hall_data[bin] ,Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]),"lpf");
    
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
    gPad->SetLogx(1);
    // gPad->SetLogy(1);
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print(Form("../plots/corr_norme2c_jetpt_relerrorleq%.2f.pdf",corr_rel_error));
}