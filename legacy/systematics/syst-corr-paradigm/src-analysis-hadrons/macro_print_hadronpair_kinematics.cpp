#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_hadronpair_kinematics()
{
    TFile* findata = new TFile((output_folder + namef_ntuple_eec_corr).c_str());
    TFile* fin     = new TFile((output_folder + namef_ntuple_mc_eec).c_str());
    
    TNtuple* ntuple_data   = (TNtuple*) findata->Get(name_ntuple_data.c_str());
    TNtuple* ntuple_mc     = (TNtuple*) fin->Get(name_ntuple_mc.c_str());
    TNtuple* ntuple_mcreco = (TNtuple*) fin->Get(name_ntuple_mcreco.c_str());

    TH1F* h_rl_data   = new TH1F("h_rl_data"   ,"",50,rl_absmin,rl_absmax);
    TH1F* h_rl_mc     = new TH1F("h_rl_mc"     ,"",50,rl_absmin,rl_absmax);
    TH1F* h_rl_mcreco = new TH1F("h_rl_mcreco" ,"",50,rl_absmin,rl_absmax);
    
    TH1F* h_weight_data   = new TH1F("h_weight_data"   ,"",50,weight_min,weight_max);
    TH1F* h_weight_mc     = new TH1F("h_weight_mc"     ,"",50,weight_min,weight_max);
    TH1F* h_weight_mcreco = new TH1F("h_weight_mcreco" ,"",50,weight_min,weight_max);

    set_histogram_style(h_rl_data   , 875 , std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_rl_mc     , 797 , std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_rl_mcreco , 868 , std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_weight_data  , 875 , std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_weight_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(h_weight_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);

    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    THStack* s_rl     = new THStack();
    TLegend* l_rl     = new TLegend();
    THStack* s_weight = new THStack();
    TLegend* l_weight = new TLegend();

    ntuple_data->Project("h_rl_data","R_L");
    ntuple_mc->Project("h_rl_mc","R_L");
    ntuple_mcreco->Project("h_rl_mcreco","R_L");
    
    ntuple_data->Project("h_weight_data","weight_pt");
    ntuple_mc->Project("h_weight_mc","weight_pt");
    ntuple_mcreco->Project("h_weight_mcreco","weight_pt");
    
    h_rl_data->Scale(1./h_rl_data->Integral());
    h_rl_mc->Scale(1./h_rl_mc->Integral());
    h_rl_mcreco->Scale(1./h_rl_mcreco->Integral());

    h_weight_data->Scale(1./h_weight_data->Integral());
    h_weight_mc->Scale(1./h_weight_mc->Integral());
    h_weight_mcreco->Scale(1./h_weight_mcreco->Integral());

    s_rl->Add(h_rl_data);
    s_rl->Add(h_rl_mc);
    s_rl->Add(h_rl_mcreco);
    s_weight->Add(h_weight_data);
    s_weight->Add(h_weight_mc);
    s_weight->Add(h_weight_mcreco);
    
    l_rl->AddEntry(h_rl_data,"data","lpf");
    l_rl->AddEntry(h_rl_mc,"mc","lpf");
    l_rl->AddEntry(h_rl_mcreco,"mcreco","lpf");
    l_weight->AddEntry(h_weight_data,"data","lpf");
    l_weight->AddEntry(h_weight_mc,"mc","lpf");
    l_weight->AddEntry(h_weight_mcreco,"mcreco","lpf");
    
    s_rl->Draw("NOSTACK");
    s_rl->SetTitle(";R_{L};Normalized Distributions");
    gPad->SetLogx(0);
    gPad->SetLogy(1);
    l_rl->Draw("SAME");

    c->Print("./plots/rl_kinematics.pdf");

    s_weight->Draw("NOSTACK");
    s_weight->SetTitle(";w;Normalized Distributions");
    gPad->SetLogx(0);
    gPad->SetLogy(1);
    l_weight->Draw("SAME"); 

    c->Print("./plots/weight_kinematics.pdf");
}