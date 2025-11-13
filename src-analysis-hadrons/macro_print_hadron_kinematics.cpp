#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.cpp"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

void macro_print_hadron_kinematics()
{
        TFile* fin     = new TFile((output_folder + namef_ntuple_hadron).c_str());
        
        TNtuple* ntuple_data   = (TNtuple*) fin->Get(name_ntuple_data.c_str());
        TNtuple* ntuple_mc     = (TNtuple*) fin->Get(name_ntuple_mc.c_str());
        TNtuple* ntuple_mcreco = (TNtuple*) fin->Get(name_ntuple_mcreco.c_str());

        TH1F* h_p_data   = new TH1F("h_p_data"   ,"",50,track_p_min,track_p_max);
        TH1F* h_p_mc     = new TH1F("h_p_mc"     ,"",50,track_p_min,track_p_max);
        TH1F* h_p_mcreco = new TH1F("h_p_mcreco" ,"",50,track_p_min,track_p_max);
        
        TH1F* h_pt_data   = new TH1F("h_pt_data"   ,"",50,0,100);
        TH1F* h_pt_mc     = new TH1F("h_pt_mc"     ,"",50,0,100);
        TH1F* h_pt_mcreco = new TH1F("h_pt_mcreco" ,"",50,0,100);
        
        TH1F* h_eta_data   = new TH1F("h_eta_data"   ,"",50,eta_min,eta_max);
        TH1F* h_eta_mc     = new TH1F("h_eta_mc"     ,"",50,eta_min,eta_max);
        TH1F* h_eta_mcreco = new TH1F("h_eta_mcreco" ,"",50,eta_min,eta_max);

        set_histogram_style(h_p_data    , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_p_mc      , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_p_mcreco  , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_pt_data   , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_pt_mc     , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_pt_mcreco , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eta_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eta_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eta_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);

        TCanvas* c = new TCanvas("c","",800,600);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        THStack* s_p   = new THStack();
        TLegend* l_p   = new TLegend();
        THStack* s_pt  = new THStack();
        TLegend* l_pt  = new TLegend();
        THStack* s_eta = new THStack();
        TLegend* l_eta = new TLegend();

        ntuple_data->Project("h_p_data","h_p");
        ntuple_mc->Project("h_p_mc","h_p");
        ntuple_mcreco->Project("h_p_mcreco","h_p");
        
        ntuple_data->Project("h_pt_data","h_pt");
        ntuple_mc->Project("h_pt_mc","h_pt");
        ntuple_mcreco->Project("h_pt_mcreco","h_pt");
        
        ntuple_data->Project("h_eta_data","h_eta");
        ntuple_mc->Project("h_eta_mc","h_eta");
        ntuple_mcreco->Project("h_eta_mcreco","h_eta");
        
        h_p_data->Scale(1./h_p_data->Integral());
        h_p_mc->Scale(1./h_p_mc->Integral());
        h_p_mcreco->Scale(1./h_p_mcreco->Integral());

        h_pt_data->Scale(1./h_pt_data->Integral());
        h_pt_mc->Scale(1./h_pt_mc->Integral());
        h_pt_mcreco->Scale(1./h_pt_mcreco->Integral());
        
        h_eta_data->Scale(1./h_eta_data->Integral());
        h_eta_mc->Scale(1./h_eta_mc->Integral());
        h_eta_mcreco->Scale(1./h_eta_mcreco->Integral());

        s_p->Add(h_p_data);
        s_p->Add(h_p_mc);
        s_p->Add(h_p_mcreco);
        s_pt->Add(h_pt_data);
        s_pt->Add(h_pt_mc);
        s_pt->Add(h_pt_mcreco);
        s_eta->Add(h_eta_data);
        s_eta->Add(h_eta_mc);
        s_eta->Add(h_eta_mcreco);
        
        l_p->AddEntry(h_p_data,"data","lpf");
        l_p->AddEntry(h_p_mc,"mc","lpf");
        l_p->AddEntry(h_p_mcreco,"mcreco","lpf");
        l_pt->AddEntry(h_pt_data,"data","lpf");
        l_pt->AddEntry(h_pt_mc,"mc","lpf");
        l_pt->AddEntry(h_pt_mcreco,"mcreco","lpf");
        l_eta->AddEntry(h_eta_data,"data","lpf");
        l_eta->AddEntry(h_eta_mc,"mc","lpf");
        l_eta->AddEntry(h_eta_mcreco,"mcreco","lpf");
        
        s_p->Draw("NOSTACK");
        s_p->SetTitle(";p(GeV);Normalized Distributions");
        gPad->SetLogx(0);
        gPad->SetLogy(1);
        l_p->Draw("SAME");

        c->Print("./plots/hadron_p_kinematics.pdf");

        s_pt->Draw("NOSTACK");
        s_pt->SetTitle(";p_{t}(GeV);Normalized Distributions");
        gPad->SetLogx(0);
        gPad->SetLogy(1);
        l_pt->Draw("SAME");

        c->Print("./plots/hadron_pt_kinematics.pdf");

        s_eta->Draw("NOSTACK");
        s_eta->SetTitle(";#eta;Normalized Distributions");
        gPad->SetLogx(0);
        gPad->SetLogy(0);
        l_eta->Draw("SAME"); 

        c->Print("./plots/hadron_eta_kinematics.pdf");
}