#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_jes(const int nbin = 40, double ptratio_min = 0.4 , double ptratio_max = 1.6)
{
    // Open the necessary files
    TFile* f = new TFile((output_folder + namef_ntuple_jes_jer).c_str());
    
    // Get the corresponding Ntuples
    TNtuple* ntuple_jes_data = (TNtuple*) f->Get((name_ntuple_jes_data).c_str());
    TNtuple* ntuple_jes_reco = (TNtuple*) f->Get((name_ntuple_jes_reco).c_str());
    
    // Define the necessary histograms to calculate purity
    TH1F* hjes_data[Nbin_jet_pt]; 
    TH1F* hjes_reco[Nbin_jet_pt];  
    TH1F* hjes_ratio[Nbin_jet_pt];  

    THStack* hs[Nbin_jet_pt];
    TLegend* l[Nbin_jet_pt];
    TCanvas* c = new TCanvas("c","",1500,500);
    c->Draw();
    c->Divide(3,1);
    
    for (int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        hs[bin] = new THStack();
        l[bin]  = new TLegend(0.6,0.8,0.75,0.9,Form(" %.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]));
        hjes_data[bin] = new TH1F(Form("hjes_data[%i]",bin) ,"",nbin , ptratio_min , ptratio_max); 
        hjes_reco[bin] = new TH1F(Form("hjes_reco[%i]",bin) ,"",nbin , ptratio_min , ptratio_max); 
    
        // set_histogram_style(hjes_data[bin],corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[0], std_marker_size);
        // set_histogram_style(hjes_reco[bin],std_marker_color_jet_pt[bin] , std_line_width, std_marker_style_jet_pt[0] , std_marker_size);
        set_histogram_style(hjes_data[bin],1, std_line_width, corr_marker_style_jet_pt[0], std_marker_size);
        set_histogram_style(hjes_reco[bin],2, std_line_width, corr_marker_style_jet_pt[0], std_marker_size);

        ntuple_jes_data->Project(Form("hjes_data[%i]",bin) ,"jet_pt/z_pt",pair_jetpt_cut[bin]);
        ntuple_jes_reco->Project(Form("hjes_reco[%i]" ,bin),"jet_pt/z_pt",pair_jetpt_cut[bin]);

        hjes_data[bin]->Scale(1./hjes_data[bin]->Integral());   
        hjes_reco[bin]->Scale(1./hjes_reco[bin]->Integral());

        c->cd(bin+1);
        
        hs[bin]->Add(hjes_data[bin]);
        hs[bin]->Add(hjes_reco[bin]);
        hs[bin]->SetTitle(";p^{jet}_{t}/p^{Z}_{t};");
        hs[bin]->Draw("NOSTACK");
        l[bin]->AddEntry(hjes_data[bin],"Normalized Data","p");
        l[bin]->AddEntry(hjes_reco[bin],"Normalized Reco","p");
        l[bin]->Draw("SAME");
    }

    c->Print(Form("./plots/jet_jes_ptratios.pdf"));

    for (int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        hjes_ratio[bin] = new TH1F(Form("hjes_ratio[%i]",bin) ,"",nbin , ptratio_min , ptratio_max);
        hjes_ratio[bin]->Divide(hjes_data[bin],hjes_reco[bin]);
        set_histogram_style(hjes_ratio[bin],corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[0], std_marker_size);
        
        c->cd(bin+1);
        hjes_ratio[bin]->Draw();
        hjes_ratio[bin]->SetTitle(";p^{jet}_{t}/p^{Z}_{t};");
        l[bin]->Clear();
        l[bin]->SetHeader(Form(" %.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]));
        l[bin]->AddEntry(hjes_ratio[bin],"Data/Reco","p");
        l[bin]->Draw("SAME");
    }

    c->Print(Form("./plots/jet_jes_datareco_ratio.pdf"));
}