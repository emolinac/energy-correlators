#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_mc_at_eec()
{
    // Open the necessary files
    TFile* fmc = new TFile((output_folder + namef_ntuple_mc_at_eec).c_str());
    
    // TNtuple* ntuple_mcreco     = (TNtuple*) fmc->Get((name_ntuple_mcreco).c_str());
    // TNtuple* ntuple_mcreco_jet = (TNtuple*) fmc->Get((name_ntuple_mcreco_jet).c_str());
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
    
    for (int bin = 0 ; bin < nbin_z_pt ; bin++)
    {
        c->cd(bin + 1);

        s_data[bin] = new THStack();
        l_data[bin] = new TLegend();

        hmc_jet[bin] = new TH1F(Form("hmc_jet[%i]" ,bin),"",1,z_pt_binning[bin],z_pt_binning[bin + 1]); 
        hmc[bin]     = new TH1F(Form("hmc[%i]" ,bin)    ,"",nbin_rl_nominal,rl_binning_at);
        
        set_histogram_style(hmc[bin], std_marker_color_jet_pt[bin] , std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);

        // Project into the histograms
        ntuple_mc->Project(Form("hmc[%i]" ,bin),"R_L",eec_zpt_cut_weightpt[bin]);
        ntuple_mc_jet->Project(Form("hmc_jet[%i]" ,bin),"z_pt",pair_zpt_cut[bin]);
        
        hmc[bin]->Scale(1./hmc_jet[bin]->Integral());

        // s_data[bin]->Add(hmcreco[bin],"E");
        s_data[bin]->Add(hmc[bin],"E");
        l_data[bin]->SetHeader(Form("%.1f<p^{Z}_{t}<%.1f (GeV)",z_pt_binning[bin],z_pt_binning[bin + 1]));
        l_data[bin]->AddEntry(hmc[bin]    ,"MC","lpf");
        s_data[bin]->Draw("NOSTACK");
        s_data[bin]->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        s_data[bin]->SetMinimum(1E-6);
        s_data[bin]->SetMaximum(3E-4);
        l_data[bin]->Draw("SAME");
    }
    
    // c->Print("./plots/mc_mcreco_eec_at_fullrange.pdf");
}