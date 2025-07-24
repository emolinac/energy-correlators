#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_mc_e2c()
{
    // Open the necessary files
    TFile* fmc = new TFile((output_folder + namef_ntuple_mc_e2c).c_str());
    
    TNtuple* ntuple_mcreco     = (TNtuple*) fmc->Get((name_ntuple_mcreco).c_str());
    TNtuple* ntuple_mcreco_jet = (TNtuple*) fmc->Get((name_ntuple_mcreco_jet).c_str());
    TNtuple* ntuple_mc         = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
    TNtuple* ntuple_mc_jet     = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());
    
    TH1F* hmcreco[nbin_jet_pt]; 
    TH1F* hmcreco_jet[nbin_jet_pt]; 
    
    TH1F* hmc[nbin_jet_pt]; 
    TH1F* hmc_jet[nbin_jet_pt]; 
    
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

        hmcreco_jet[bin] = new TH1F(Form("hmcreco_jet[%i]" ,bin),"",1,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
        hmcreco[bin]     = new TH1F(Form("hmcreco[%i]",bin),"",nbin_rl,rl_binning);
        
        hmc_jet[bin] = new TH1F(Form("hmc_jet[%i]" ,bin),"",1,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
        hmc[bin]     = new TH1F(Form("hmc[%i]" ,bin)       ,"",nbin_rl,rl_binning);
        
        set_histogram_style(hmcreco[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size);
        set_histogram_style(hmc[bin]         , std_marker_color_jet_pt[bin] , std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);

        // Project into the histograms
        ntuple_mcreco->Project(Form("hmcreco[%i]",bin),"R_L",e2c_jet_pt_cut_weightpt[bin]);
        ntuple_mc->Project(Form("hmc[%i]" ,bin),"R_L",e2c_jet_pt_cut_weightpt[bin]);
        
        ntuple_mcreco_jet->Project(Form("hmcreco_jet[%i]" ,bin), "jet_pt",pair_jet_pt_cut[bin]);
        ntuple_mc_jet->Project(Form("hmc_jet[%i]" ,bin), "jet_pt",pair_jet_pt_cut[bin]);
        
        hmcreco[bin]->Scale(1./hmcreco_jet[bin]->Integral());
        hmc[bin]->Scale(1./hmc_jet[bin]->Integral());

        s_data[bin]->Add(hmcreco[bin],"E");
        s_data[bin]->Add(hmc[bin],"E");
        l_data[bin]->SetHeader(Form("%.1f<p_{T,jet}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
        l_data[bin]->AddEntry(hmcreco[bin],"MCReco","lpf");
        l_data[bin]->AddEntry(hmc[bin]    ,"MC","lpf");
        s_data[bin]->Draw("NOSTACK");
        s_data[bin]->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        gPad->SetLogx(1);
        // gPad->SetLogy(1);
        l_data[bin]->Draw("SAME");
    }
    
    // tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print("./plots/mc_mcreco_e2c_fullrange.pdf");
}