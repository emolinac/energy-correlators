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

void macro_print_hadron_kinematics()
{
        TFile* fin = new TFile((output_folder + namef_ntuple_hadron_jet).c_str());
        
        TNtuple* ntuple_data   = (TNtuple*) fin->Get(name_ntuple_data.c_str());
        TNtuple* ntuple_mc     = (TNtuple*) fin->Get(name_ntuple_mc.c_str());
        TNtuple* ntuple_mcreco = (TNtuple*) fin->Get(name_ntuple_mcreco.c_str());

        TH1F* h_pt_data   = new TH1F("h_pt_data"   ,"",50,0,100);
        TH1F* h_pt_mc     = new TH1F("h_pt_mc"     ,"",50,0,100);
        TH1F* h_pt_mcreco = new TH1F("h_pt_mcreco" ,"",50,0,100);
        
        TH1F* h_p_data   = new TH1F("h_p_data"   ,"",50,track_p_min,track_p_max);
        TH1F* h_p_mc     = new TH1F("h_p_mc"     ,"",50,track_p_min,track_p_max);
        TH1F* h_p_mcreco = new TH1F("h_p_mcreco" ,"",50,track_p_min,track_p_max);
        
        TH1F* h_y_data   = new TH1F("h_y_data"  ,"",15,lhcb_eta_min,lhcb_eta_max);
        TH1F* h_y_mc     = new TH1F("h_y_mc"    ,"",15,lhcb_eta_min,lhcb_eta_max);
        TH1F* h_y_mcreco = new TH1F("h_y_mcreco","",15,lhcb_eta_min,lhcb_eta_max);

        TH1F* h_eta_data   = new TH1F("h_eta_data"  ,"",15,lhcb_eta_min,lhcb_eta_max);
        TH1F* h_eta_mc     = new TH1F("h_eta_mc"    ,"",15,lhcb_eta_min,lhcb_eta_max);
        TH1F* h_eta_mcreco = new TH1F("h_eta_mcreco","",15,lhcb_eta_min,lhcb_eta_max);

        TH1F* h_phi_data   = new TH1F("h_phi_data"  ,"",30,-TMath::Pi(),TMath::Pi());
        TH1F* h_phi_mc     = new TH1F("h_phi_mc"    ,"",30,-TMath::Pi(),TMath::Pi());
        TH1F* h_phi_mcreco = new TH1F("h_phi_mcreco","",30,-TMath::Pi(),TMath::Pi());

        set_histogram_style(h_pt_data   , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_pt_mc     , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_pt_mcreco , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_p_data    , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_p_mc      , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_p_mcreco  , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eta_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eta_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eta_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_y_data    , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_y_mc      , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_y_mcreco  , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_phi_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_phi_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_phi_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);

        ntuple_data->Project("h_pt_data", "h_pt","");
        ntuple_mc->Project("h_pt_mc", "h_pt","");
        ntuple_mcreco->Project("h_pt_mcreco", "h_pt","");
        
        ntuple_data->Project("h_p_data", "h_p","");
        ntuple_mc->Project("h_p_mc", "h_p","");
        ntuple_mcreco->Project("h_p_mcreco", "h_p","");

        ntuple_data->Project("h_eta_data","h_eta","");
        ntuple_mc->Project("h_eta_mc","h_eta","");
        ntuple_mcreco->Project("h_eta_mcreco","h_eta","");
        
        ntuple_data->Project("h_y_data","h_y","");
        ntuple_mc->Project("h_y_mc","h_y","");
        ntuple_mcreco->Project("h_y_mcreco","h_y","");
        
        ntuple_data->Project("h_phi_data","h_phi","");
        ntuple_mc->Project("h_phi_mc","h_phi","");
        ntuple_mcreco->Project("h_phi_mcreco","h_phi","");
        
        h_pt_data->Scale(1./h_pt_data->Integral());
        h_pt_mc->Scale(1./h_pt_mc->Integral());
        h_pt_mcreco->Scale(1./h_pt_mcreco->Integral());
        h_p_data->Scale(1./h_p_data->Integral());
        h_p_mc->Scale(1./h_p_mc->Integral());
        h_p_mcreco->Scale(1./h_p_mcreco->Integral());
        h_eta_data->Scale(1./h_eta_data->Integral());
        h_eta_mc->Scale(1./h_eta_mc->Integral());
        h_eta_mcreco->Scale(1./h_eta_mcreco->Integral());
        h_y_data->Scale(1./h_y_data->Integral());
        h_y_mc->Scale(1./h_y_mc->Integral());
        h_y_mcreco->Scale(1./h_y_mcreco->Integral());
        h_phi_data->Scale(1./h_phi_data->Integral());
        h_phi_mc->Scale(1./h_phi_mc->Integral());
        h_phi_mcreco->Scale(1./h_phi_mcreco->Integral());
        
        TCanvas* c = new TCanvas("c","",800,600);
        c->Draw();

        THStack* h_pt  = new THStack();
        h_pt->Add(h_pt_data);
        h_pt->Add(h_pt_mc);
        h_pt->Add(h_pt_mcreco);
        h_pt->SetTitle(";p_{T}(GeV);Normalized Distributions");
        
        THStack* h_p  = new THStack();
        h_p->Add(h_p_data);
        h_p->Add(h_p_mc);
        h_p->Add(h_p_mcreco);
        h_p->SetTitle(";p(GeV);Normalized Distributions");
        
        THStack* heta = new THStack();
        heta->Add(h_eta_data);
        heta->Add(h_eta_mc);
        heta->Add(h_eta_mcreco);
        heta->SetTitle(";#eta_{jet};Normalized Distributions");

        THStack* hy = new THStack();
        hy->Add(h_y_data);
        hy->Add(h_y_mc);
        hy->Add(h_y_mcreco);
        hy->SetTitle(";y_{jet};Normalized Distributions");

        THStack* hphi = new THStack();
        hphi->Add(h_phi_data);
        hphi->Add(h_phi_mc);
        hphi->Add(h_phi_mcreco);
        hphi->SetTitle(";#phi_{jet};Normalized Distributions");

        TLegend* l = new TLegend();
        l->AddEntry(h_pt_data,"data","lpf");
        l->AddEntry(h_pt_mc,"mc","lpf");
        l->AddEntry(h_pt_mcreco,"mcreco","lpf");

        h_pt->Draw("NOSTACK");
        gPad->SetLogy(1);
        l->Draw("SAME");
        c->Print("./plots/hadron_pt_kinematics.pdf");

        h_p->Draw("NOSTACK");
        gPad->SetLogy(1);
        l->Draw("SAME");
        c->Print("./plots/hadron_p_kinematics.pdf");

        heta->Draw("NOSTACK");
        gPad->SetLogy(0);
        l->Draw("SAME");
        c->Print("./plots/hadron_eta_kinematics.pdf");

        hy->Draw("NOSTACK");
        gPad->SetLogy(0);
        l->Draw("SAME");
        c->Print("./plots/hadron_y_kinematics.pdf");

        hphi->Draw("NOSTACK");
        gPad->SetLogy(0);
        l->Draw("SAME");
        c->Print("./plots/hadron_phi_kinematics.pdf");
}