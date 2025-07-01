#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_mc_chargede2c()
{
    // Open the necessary files
    TFile* fmc = new TFile((output_folder+namef_ntuple_mc_e2c).c_str());
    
    TNtuple* ntuple_mcreco     = (TNtuple*) fmc->Get((name_ntuple_mcreco).c_str());
    TNtuple* ntuple_mcreco_jet = (TNtuple*) fmc->Get((name_ntuple_mcreco_jet).c_str());
    TNtuple* ntuple_mc         = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
    TNtuple* ntuple_mc_jet     = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());
    
    TH1F* hmcreco_pp[Nbin_jet_pt]; 
    TH1F* hmcreco_pm[Nbin_jet_pt]; 
    TH1F* hmcreco_mm[Nbin_jet_pt]; 
    TH1F* hmcreco_all[Nbin_jet_pt]; 
    
    TH1F* hmc_pp[Nbin_jet_pt]; 
    TH1F* hmc_pm[Nbin_jet_pt]; 
    TH1F* hmc_mm[Nbin_jet_pt]; 
    TH1F* hmc_all[Nbin_jet_pt];
    
    TCanvas* c = new TCanvas("c","",1800,600);
    c->Draw();
    c->Divide(3,1);
    
    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    THStack* s_data[3];
    TLegend* l_data[3];
    
    for (int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        c->cd(bin+1);

        s_data[bin] = new THStack();
        l_data[bin] = new TLegend();

        hmcreco_pp[bin] = new TH1F(Form("hmcreco_pp[%i]",bin)  ,"",Nbin_R_L,rl_binning);
        hmcreco_pm[bin] = new TH1F(Form("hmcreco_pm[%i]",bin)  ,"",Nbin_R_L,rl_binning);
        hmcreco_mm[bin] = new TH1F(Form("hmcreco_mm[%i]",bin)  ,"",Nbin_R_L,rl_binning);
        hmcreco_all[bin] = new TH1F(Form("hmcreco_all[%i]",bin),"",Nbin_R_L,rl_binning);
        
        hmc_pp[bin] = new TH1F(Form("hmc_pp[%i]" ,bin)  ,"",Nbin_R_L,rl_binning);
        hmc_pm[bin] = new TH1F(Form("hmc_pm[%i]" ,bin)  ,"",Nbin_R_L,rl_binning);
        hmc_mm[bin] = new TH1F(Form("hmc_mm[%i]" ,bin)  ,"",Nbin_R_L,rl_binning);
        hmc_all[bin] = new TH1F(Form("hmc_all[%i]" ,bin),"",Nbin_R_L,rl_binning);
        
        // Project into the histograms
        ntuple_mcreco->Project(Form("hmcreco_pp[%i]",bin),"R_L",e2c_jetpt_cut_weightpt_pp[bin]);
        ntuple_mcreco->Project(Form("hmcreco_pm[%i]",bin),"R_L",e2c_jetpt_cut_weightpt_pm[bin]);
        ntuple_mcreco->Project(Form("hmcreco_mm[%i]",bin),"R_L",e2c_jetpt_cut_weightpt_mm[bin]);
        ntuple_mcreco->Project(Form("hmcreco_all[%i]",bin),"R_L",e2c_jetpt_cut_weightpt[bin]);
        hmcreco_pp[bin]->Add(hmcreco_mm[bin]);
        hmcreco_pp[bin]->Divide(hmcreco_all[bin]);
        hmcreco_pm[bin]->Divide(hmcreco_all[bin]);

        ntuple_mc->Project(Form("hmc_pp[%i]" ,bin),"R_L",e2c_jetpt_cut_weightpt_pp[bin]);
        ntuple_mc->Project(Form("hmc_pm[%i]" ,bin),"R_L",e2c_jetpt_cut_weightpt_pm[bin]);
        ntuple_mc->Project(Form("hmc_mm[%i]" ,bin),"R_L",e2c_jetpt_cut_weightpt_mm[bin]);
        ntuple_mc->Project(Form("hmc_all[%i]" ,bin),"R_L",e2c_jetpt_cut_weightpt[bin]);
        hmc_pp[bin]->Add(hmc_mm[bin]);
        hmc_pp[bin]->Divide(hmc_all[bin]);
        hmc_pm[bin]->Divide(hmc_all[bin]);

        set_histogram_style(hmcreco_pp[bin], corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
        set_histogram_style(hmcreco_pm[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size);
        set_histogram_style(hmc_pp[bin]    , std_marker_color_jet_pt[bin] , std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
        set_histogram_style(hmc_pm[bin]    , std_marker_color_jet_pt[bin] , std_line_width, corr_marker_style_jet_pt[bin], std_marker_size);

        // s_data[bin]->Add(hmcreco_pp[bin],"E");
        // s_data[bin]->Add(hmcreco_pm[bin],"E");
        s_data[bin]->Add(hmc_pp[bin],"E");
        s_data[bin]->Add(hmc_pm[bin],"E");
        l_data[bin]->SetHeader(Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]));
        // l_data[bin]->AddEntry(hmcreco_pp[bin],"MCReco Eq. Charge","lpf");
        // l_data[bin]->AddEntry(hmcreco_pm[bin],"MCReco Op. Charge","lpf");
        l_data[bin]->AddEntry(hmc_pp[bin]    ,"MC Eq. Charge","lpf");
        l_data[bin]->AddEntry(hmc_pm[bin]    ,"MC Op. Charge","lpf");
        s_data[bin]->Draw("NOSTACK");
        s_data[bin]->SetTitle(";R_{L};");
        s_data[bin]->SetMinimum(0.2);
        s_data[bin]->SetMaximum(0.8);
        // gPad->SetLogx(1);
        // gPad->SetLogy(1);
        l_data[bin]->Draw("SAME");
    }
    
    // tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print("./plots/mc_mcreco_chargede2c_rlleqjetradius_biggerbins_mc.pdf");
}
