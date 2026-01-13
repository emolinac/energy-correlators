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

void macro_print_jet_kinematics()
{
        TFile* fin = new TFile((output_folder + namef_ntuple_hadron_jet).c_str());
        
        TNtuple* ntuple_jet_data   = (TNtuple*) fin->Get(name_ntuple_data_jet.c_str());
        TNtuple* ntuple_jet_mc     = (TNtuple*) fin->Get(name_ntuple_mc_jet.c_str());
        TNtuple* ntuple_jet_mcreco = (TNtuple*) fin->Get(name_ntuple_mcreco_jet.c_str());

        TH1F* h_jet_pt_data   = new TH1F("h_jet_pt_data"   ,"",20,jet_pt_binning[0],jet_pt_binning[nbin_jet_pt]);
        TH1F* h_jet_pt_mc     = new TH1F("h_jet_pt_mc"     ,"",20,jet_pt_binning[0],jet_pt_binning[nbin_jet_pt]);
        TH1F* h_jet_pt_mcreco = new TH1F("h_jet_pt_mcreco" ,"",20,jet_pt_binning[0],jet_pt_binning[nbin_jet_pt]);
        
        TH1F* h_jet_y_data   = new TH1F("h_jet_y_data"  ,"",12,jet_eta_min,jet_eta_max);
        TH1F* h_jet_y_mc     = new TH1F("h_jet_y_mc"    ,"",12,jet_eta_min,jet_eta_max);
        TH1F* h_jet_y_mcreco = new TH1F("h_jet_y_mcreco","",12,jet_eta_min,jet_eta_max);

        TH1F* h_jet_eta_data   = new TH1F("h_jet_eta_data"  ,"",12,jet_eta_min,jet_eta_max);
        TH1F* h_jet_eta_mc     = new TH1F("h_jet_eta_mc"    ,"",12,jet_eta_min,jet_eta_max);
        TH1F* h_jet_eta_mcreco = new TH1F("h_jet_eta_mcreco","",12,jet_eta_min,jet_eta_max);

        TH1F* h_jet_phi_data   = new TH1F("h_jet_phi_data"  ,"",30,-TMath::Pi(),TMath::Pi());
        TH1F* h_jet_phi_mc     = new TH1F("h_jet_phi_mc"    ,"",30,-TMath::Pi(),TMath::Pi());
        TH1F* h_jet_phi_mcreco = new TH1F("h_jet_phi_mcreco","",30,-TMath::Pi(),TMath::Pi());

        TH1F* h_z_pt_data   = new TH1F("h_z_pt_data"   ,"",20,0,200);
        TH1F* h_z_pt_mc     = new TH1F("h_z_pt_mc"     ,"",20,0,200);
        TH1F* h_z_pt_mcreco = new TH1F("h_z_pt_mcreco" ,"",20,0,200);
        
        TH1F* h_z_y_data   = new TH1F("h_z_y_data"  ,"",12,eta_min,eta_max);
        TH1F* h_z_y_mc     = new TH1F("h_z_y_mc"    ,"",12,eta_min,eta_max);
        TH1F* h_z_y_mcreco = new TH1F("h_z_y_mcreco","",12,eta_min,eta_max);

        TH1F* h_z_eta_data   = new TH1F("h_z_eta_data"  ,"",12,eta_min,eta_max);
        TH1F* h_z_eta_mc     = new TH1F("h_z_eta_mc"    ,"",12,eta_min,eta_max);
        TH1F* h_z_eta_mcreco = new TH1F("h_z_eta_mcreco","",12,eta_min,eta_max);

        TH1F* h_z_phi_data   = new TH1F("h_z_phi_data"  ,"",30,-TMath::Pi(),TMath::Pi());
        TH1F* h_z_phi_mc     = new TH1F("h_z_phi_mc"    ,"",30,-TMath::Pi(),TMath::Pi());
        TH1F* h_z_phi_mcreco = new TH1F("h_z_phi_mcreco","",30,-TMath::Pi(),TMath::Pi());

        set_histogram_style(h_jet_pt_data   , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_pt_mc     , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_pt_mcreco , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_eta_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_eta_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_eta_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_y_data    , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_y_mc      , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_y_mcreco  , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_phi_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_phi_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_jet_phi_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);

        set_histogram_style(h_z_pt_data   , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_pt_mc     , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_pt_mcreco , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_eta_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_eta_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_eta_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_y_data    , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_y_mc      , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_y_mcreco  , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_phi_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_phi_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_phi_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);

        ntuple_jet_data->Project("h_jet_pt_data", "jet_pt","");
        ntuple_jet_mc->Project("h_jet_pt_mc", "jet_pt","");
        ntuple_jet_mcreco->Project("h_jet_pt_mcreco", "jet_pt","");
        
        ntuple_jet_data->Project("h_jet_eta_data","jet_eta","");
        ntuple_jet_mc->Project("h_jet_eta_mc","jet_eta","");
        ntuple_jet_mcreco->Project("h_jet_eta_mcreco","jet_eta","");
        
        ntuple_jet_data->Project("h_jet_y_data","jet_y","");
        ntuple_jet_mc->Project("h_jet_y_mc","jet_y","");
        ntuple_jet_mcreco->Project("h_jet_y_mcreco","jet_y","");
        
        ntuple_jet_data->Project("h_jet_phi_data","jet_phi","");
        ntuple_jet_mc->Project("h_jet_phi_mc","jet_phi","");
        ntuple_jet_mcreco->Project("h_jet_phi_mcreco","jet_phi","");
        
        ntuple_jet_data->Project("h_z_pt_data","z_pt","");
        ntuple_jet_mc->Project("h_z_pt_mc","z_pt","");
        ntuple_jet_mcreco->Project("h_z_pt_mcreco","z_pt","");
        
        ntuple_jet_data->Project("h_z_eta_data","z_eta","");
        ntuple_jet_mc->Project("h_z_eta_mc","z_eta","");
        ntuple_jet_mcreco->Project("h_z_eta_mcreco","z_eta","");

        ntuple_jet_data->Project("h_z_y_data","z_y","");
        ntuple_jet_mc->Project("h_z_y_mc","z_y","");
        ntuple_jet_mcreco->Project("h_z_y_mcreco","z_y","");

        ntuple_jet_data->Project("h_z_phi_data","z_phi","");
        ntuple_jet_mc->Project("h_z_phi_mc","z_phi","");
        ntuple_jet_mcreco->Project("h_z_phi_mcreco","z_phi","");

        h_jet_pt_data->Scale(1./h_jet_pt_data->Integral());
        h_jet_pt_mc->Scale(1./h_jet_pt_mc->Integral());
        h_jet_pt_mcreco->Scale(1./h_jet_pt_mcreco->Integral());
        h_jet_eta_data->Scale(1./h_jet_eta_data->Integral());
        h_jet_eta_mc->Scale(1./h_jet_eta_mc->Integral());
        h_jet_eta_mcreco->Scale(1./h_jet_eta_mcreco->Integral());
        h_jet_y_data->Scale(1./h_jet_y_data->Integral());
        h_jet_y_mc->Scale(1./h_jet_y_mc->Integral());
        h_jet_y_mcreco->Scale(1./h_jet_y_mcreco->Integral());
        h_jet_phi_data->Scale(1./h_jet_phi_data->Integral());
        h_jet_phi_mc->Scale(1./h_jet_phi_mc->Integral());
        h_jet_phi_mcreco->Scale(1./h_jet_phi_mcreco->Integral());
        
        h_z_pt_data->Scale(1./h_z_pt_data->Integral());
        h_z_pt_mc->Scale(1./h_z_pt_mc->Integral());
        h_z_pt_mcreco->Scale(1./h_z_pt_mcreco->Integral());
        h_z_eta_data->Scale(1./h_z_eta_data->Integral());
        h_z_eta_mc->Scale(1./h_z_eta_mc->Integral());
        h_z_eta_mcreco->Scale(1./h_z_eta_mcreco->Integral());
        h_z_y_data->Scale(1./h_z_y_data->Integral());
        h_z_y_mc->Scale(1./h_z_y_mc->Integral());
        h_z_y_mcreco->Scale(1./h_z_y_mcreco->Integral());
        h_z_phi_data->Scale(1./h_z_phi_data->Integral());
        h_z_phi_mc->Scale(1./h_z_phi_mc->Integral());
        h_z_phi_mcreco->Scale(1./h_z_phi_mcreco->Integral());

        TCanvas* c = new TCanvas("c","",800,600);
        c->Draw();

        THStack* hjet_pt  = new THStack();
        hjet_pt->Add(h_jet_pt_data);
        hjet_pt->Add(h_jet_pt_mc);
        hjet_pt->Add(h_jet_pt_mcreco);
        hjet_pt->SetTitle(";p_{T,jet}(GeV);Normalized Distributions");
        
        THStack* hjeteta = new THStack();
        hjeteta->Add(h_jet_eta_data);
        hjeteta->Add(h_jet_eta_mc);
        hjeteta->Add(h_jet_eta_mcreco);
        hjeteta->SetTitle(";#eta_{jet};Normalized Distributions");

        THStack* hjety = new THStack();
        hjety->Add(h_jet_y_data);
        hjety->Add(h_jet_y_mc);
        hjety->Add(h_jet_y_mcreco);
        hjety->SetTitle(";y_{jet};Normalized Distributions");

        THStack* hjetphi = new THStack();
        hjetphi->Add(h_jet_phi_data);
        hjetphi->Add(h_jet_phi_mc);
        hjetphi->Add(h_jet_phi_mcreco);
        hjetphi->SetTitle(";#phi_{jet};Normalized Distributions");

        THStack* hzpt  = new THStack();
        hzpt->Add(h_z_pt_data);
        hzpt->Add(h_z_pt_mc);
        hzpt->Add(h_z_pt_mcreco);
        hzpt->SetTitle(";p^{Z}_{T}(GeV);Normalized Distributions");
        
        THStack* hzeta = new THStack();
        hzeta->Add(h_z_eta_data);
        hzeta->Add(h_z_eta_mc);
        hzeta->Add(h_z_eta_mcreco);
        hzeta->SetTitle(";#eta_{Z};Normalized Distributions");

        THStack* hzy = new THStack();
        hzy->Add(h_z_y_data);
        hzy->Add(h_z_y_mc);
        hzy->Add(h_z_y_mcreco);
        hzy->SetTitle(";y_{Z};Normalized Distributions");

        THStack* hzphi = new THStack();
        hzphi->Add(h_z_phi_data);
        hzphi->Add(h_z_phi_mc);
        hzphi->Add(h_z_phi_mcreco);
        hzphi->SetTitle(";#phi_{Z};Normalized Distributions");

        TLegend* ljet_pt = new TLegend();
        ljet_pt->AddEntry(h_jet_pt_data,"data","lpf");
        ljet_pt->AddEntry(h_jet_pt_mc,"mc","lpf");
        ljet_pt->AddEntry(h_jet_pt_mcreco,"mcreco","lpf");

        TLegend* ljeteta = new TLegend();
        ljeteta->AddEntry(h_jet_eta_data,"data","lpf");
        ljeteta->AddEntry(h_jet_eta_mc,"mc","lpf");
        ljeteta->AddEntry(h_jet_eta_mcreco,"mcreco","lpf");

        hjet_pt->Draw("NOSTACK");
        gPad->SetLogy(1);
        ljet_pt->Draw("SAME");
        c->Print("./plots/jet_pt_kinematics.pdf");

        hzpt->Draw("NOSTACK");
        gPad->SetLogy(1);
        ljet_pt->Draw("SAME");
        c->Print("./plots/z_pt_kinematics.pdf");

        hjeteta->Draw("NOSTACK");
        gPad->SetLogy(0);
        ljeteta->Draw("SAME");
        c->Print("./plots/jet_eta_kinematics.pdf");

        hzeta->Draw("NOSTACK");
        gPad->SetLogy(1);
        ljeteta->Draw("SAME");
        c->Print("./plots/z_eta_kinematics.pdf");

        hjety->Draw("NOSTACK");
        gPad->SetLogy(0);
        ljeteta->Draw("SAME");
        c->Print("./plots/jet_y_kinematics.pdf");

        hzy->Draw("NOSTACK");
        gPad->SetLogy(1);
        ljeteta->Draw("SAME");
        c->Print("./plots/z_y_kinematics.pdf");

        hjetphi->Draw("NOSTACK");
        gPad->SetLogy(0);
        ljeteta->Draw("SAME");
        c->Print("./plots/jet_phi_kinematics.pdf");

        hzphi->Draw("NOSTACK");
        gPad->SetLogy(1);
        ljeteta->Draw("SAME");
        c->Print("./plots/z_phi_kinematics.pdf");
}