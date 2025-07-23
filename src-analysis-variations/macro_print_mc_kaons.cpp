#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_mc_kaons()
{
    // Open the necessary files
    TFile* fmc = new TFile((output_folder + namef_ntuple_mc_e2c).c_str());
    
    TNtuple* ntuple_mcreco     = (TNtuple*) fmc->Get((name_ntuple_mcreco).c_str());
    TNtuple* ntuple_mcreco_jet = (TNtuple*) fmc->Get((name_ntuple_mcreco_jet).c_str());
    TNtuple* ntuple_mc         = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
    TNtuple* ntuple_mc_jet     = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());
    
    TH1F* hmcreco_kaon[nbin_jet_pt]; 
    TH1F* hmcreco_nokaon[nbin_jet_pt]; 
    TH1F* hmcreco_all[nbin_jet_pt]; 
    
    TH1F* hmc_kaon[nbin_jet_pt]; 
    TH1F* hmc_nokaon[nbin_jet_pt]; 
    TH1F* hmc_all[nbin_jet_pt];
    
    TCanvas* c = new TCanvas("c","",1800,600);
    c->Draw();
    c->Divide(3,1);
    
    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    THStack* s_data[3];
    TLegend* l_data[3];
    
    for (int bin = 0 ; bin < nbin_jet_pt ; bin++)
    {
        c->cd(bin + 1);

        s_data[bin] = new THStack();
        l_data[bin] = new TLegend();

        hmcreco_kaon[bin]   = new TH1F(Form("hmcreco_kaon[%i]",bin)  ,"",nbin_rl,rl_binning);
        hmcreco_nokaon[bin] = new TH1F(Form("hmcreco_nokaon[%i]",bin),"",nbin_rl,rl_binning);
        hmcreco_all[bin]    = new TH1F(Form("hmcreco_all[%i]",bin)   ,"",nbin_rl,rl_binning);
        
        hmc_kaon[bin]   = new TH1F(Form("hmc_kaon[%i]" ,bin)  ,"",nbin_rl,rl_binning);
        hmc_nokaon[bin] = new TH1F(Form("hmc_nokaon[%i]" ,bin),"",nbin_rl,rl_binning);
        hmc_all[bin]    = new TH1F(Form("hmc_all[%i]" ,bin)   ,"",nbin_rl,rl_binning);
        
        // Project into the histograms
        ntuple_mcreco->Project(Form("hmcreco_kaon[%i]",bin),"R_L",e2c_jet_pt_cut_weightpt_kaon[bin]);
        ntuple_mcreco->Project(Form("hmcreco_nokaon[%i]",bin),"R_L",e2c_jet_pt_cut_weightpt_nokaon[bin]);
        ntuple_mcreco->Project(Form("hmcreco_all[%i]",bin),"R_L",e2c_jet_pt_cut_weightpt[bin]);
    
        // hmcreco_kaon[bin]->Divide(hmcreco_all[bin]);
        // hmcreco_nokaon[bin]->Divide(hmcreco_all[bin]);

        ntuple_mc->Project(Form("hmc_kaon[%i]" ,bin)  ,"R_L",e2c_jet_pt_cut_weightpt_kaon[bin]);
        ntuple_mc->Project(Form("hmc_nokaon[%i]" ,bin),"R_L",e2c_jet_pt_cut_weightpt_nokaon[bin]);
        ntuple_mc->Project(Form("hmc_all[%i]" ,bin)   ,"R_L",e2c_jet_pt_cut_weightpt[bin]);
        // hmc_kaon[bin]->Divide(hmc_all[bin]);
        // hmc_nokaon[bin]->Divide(hmc_all[bin]);

        set_histogram_style(hmcreco_kaon[bin], corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
        set_histogram_style(hmcreco_nokaon[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size);
        set_histogram_style(hmc_kaon[bin]    , std_marker_color_jet_pt[bin] , std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
        set_histogram_style(hmc_nokaon[bin]    , std_marker_color_jet_pt[bin] , std_line_width, corr_marker_style_jet_pt[bin], std_marker_size);

        s_data[bin]->Add(hmcreco_kaon[bin],"E");
        s_data[bin]->Add(hmcreco_nokaon[bin],"E");
        // s_data[bin]->Add(hmc_kaon[bin],"E");
        // s_data[bin]->Add(hmc_nokaon[bin],"E");
        l_data[bin]->SetHeader(Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
        l_data[bin]->AddEntry(hmcreco_kaon[bin]  ,"MCReco Kaon","lpf");
        l_data[bin]->AddEntry(hmcreco_nokaon[bin],"MCReco No Kaons","lpf");
        // l_data[bin]->AddEntry(hmc_kaon[bin]  ,"MC Kaon","lpf");
        // l_data[bin]->AddEntry(hmc_nokaon[bin],"MC No Kaons","lpf");
        s_data[bin]->Draw("NOSTACK");
        s_data[bin]->SetTitle(";R_{L};");
        // s_data[bin]->SetMinimum(0.2);
        // s_data[bin]->SetMaximum(0.8);
        gPad->SetLogx(1);
        // gPad->SetLogy(1);
        l_data[bin]->Draw("SAME");
    }
    
    // tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    // c->Print("./plots/mc_mcreco_kaonse2c_rlleqjetradius_mcreco.pdf");
}
