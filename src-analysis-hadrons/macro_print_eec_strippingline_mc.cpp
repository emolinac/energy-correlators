#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils.cpp"
#include "../include/utils.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

void macro_print_eec_strippingline_mc()
{
        TFile* fin_28r1 = new TFile((output_folder + namef_ntuple_mc_eec_28r1).c_str());
        TFile* fin_28r2 = new TFile((output_folder + namef_ntuple_mc_eec_28r2).c_str());
        
        TNtuple* ntuple_mcreco_28r1 = (TNtuple*) fin_28r1->Get(name_ntuple_mcreco.c_str());
        TNtuple* ntuple_mcreco_28r2 = (TNtuple*) fin_28r2->Get(name_ntuple_mcreco.c_str());

        TH1F* h_eec_mc_28r1 = new TH1F("h_eec_mc_28r1","",nbin_rl_nominal,rl_nominal_binning);
        TH1F* h_eec_mc_28r2 = new TH1F("h_eec_mc_28r2","",nbin_rl_nominal,rl_nominal_binning);
        
        TH1F* h_w_mc_28r1 = new TH1F("h_w_mc_28r1","",nbin_weight,weight_binning);
        TH1F* h_w_mc_28r2 = new TH1F("h_w_mc_28r2" ,"",nbin_weight,weight_binning);
        
        set_histogram_style(h_eec_mc_28r1, 797 ,std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eec_mc_28r2, 868 ,std_line_width, std_marker_style, std_marker_size);
        
        set_histogram_style(h_w_mc_28r1, 797 ,std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_w_mc_28r2, 868 ,std_line_width, std_marker_style, std_marker_size);
        
        TCanvas* c = new TCanvas("c","",800,600);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        THStack* s_eec = new THStack();
        THStack* s_w   = new THStack();
        
        TLegend* l_eec = new TLegend();
        TLegend* l_w   = new TLegend();
        
        ntuple_mcreco_28r1->Project("h_eec_mc_28r1","R_L","weight_pt");
        ntuple_mcreco_28r2->Project("h_eec_mc_28r2","R_L","weight_pt");
        
        ntuple_mcreco_28r1->Project("h_w_mc_28r1","weight_pt","");
        ntuple_mcreco_28r2->Project("h_w_mc_28r2","weight_pt","");
        
        h_eec_mc_28r1->Scale(1./h_eec_mc_28r1->Integral(),"width");
        h_eec_mc_28r2->Scale(1./h_eec_mc_28r2->Integral(),"width");
        
        h_w_mc_28r1->Scale(1./h_w_mc_28r1->Integral(),"width");
        h_w_mc_28r2->Scale(1./h_w_mc_28r2->Integral(),"width");
        
        s_eec->Add(h_eec_mc_28r1);
        s_eec->Add(h_eec_mc_28r2);
        
        s_w->Add(h_w_mc_28r1);
        s_w->Add(h_w_mc_28r2);
        
        l_eec->AddEntry(h_eec_mc_28r1,"28r1","lpf");
        l_eec->AddEntry(h_eec_mc_28r2,"28r2","lpf");
        
        l_w->AddEntry(h_w_mc_28r1,"28r1","lpf");
        l_w->AddEntry(h_w_mc_28r2,"28r2","lpf");
        
        s_eec->Draw("NOSTACK");
        s_eec->SetTitle(";R_{L};Normalized EEC Distributions");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        l_eec->Draw("SAME");

        c->Print("./plots/eec_28rX_mcreco.pdf");

        s_w->Draw("NOSTACK");
        s_w->SetTitle(";w_{p_{T}};Normalized Pair Distributions");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        l_w->Draw("SAME");

        c->Print("./plots/weight_28rX_mcreco.pdf");
}