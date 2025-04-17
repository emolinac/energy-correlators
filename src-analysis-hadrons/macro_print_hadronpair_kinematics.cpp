#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_hadronpair_kinematics()
{
    TFile* findata = new TFile((output_folder+namef_ntuple_e2c_corr).c_str());
    TFile* fin     = new TFile((output_folder+namef_ntuple_mc_e2c).c_str());
    
    TNtuple* ntuple_data   = (TNtuple*) findata->Get(name_ntuple_corrjet.c_str());
    TNtuple* ntuple_mc     = (TNtuple*) fin->Get(name_ntuple_mc.c_str());
    TNtuple* ntuple_mcreco = (TNtuple*) fin->Get(name_ntuple_mcreco.c_str());

    TH1F* h_rl_data   = new TH1F("h_rl_data"   ,"",Nbin_R_L,rl_binning);
    TH1F* h_rl_mc     = new TH1F("h_rl_mc"     ,"",Nbin_R_L,rl_binning);
    TH1F* h_rl_mcreco = new TH1F("h_rl_mcreco" ,"",Nbin_R_L,rl_binning);
    
    TH1F* h_weight_data   = new TH1F("h_weight_data"   ,"",Nbin_weight,weight_binning);
    TH1F* h_weight_mc     = new TH1F("h_weight_mc"     ,"",Nbin_weight,weight_binning);
    TH1F* h_weight_mcreco = new TH1F("h_weight_mcreco" ,"",Nbin_weight,weight_binning);

    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

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
    s_rl->SetTitle(";R_{L};A.U.");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    l_rl->Draw("SAME");

    // tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    // c->Print("./plots/weight_jetpt.pdf");
    s_weight->Draw("NOSTACK");
    s_weight->SetTitle(";w;A.U.");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    l_weight->Draw("SAME");

}