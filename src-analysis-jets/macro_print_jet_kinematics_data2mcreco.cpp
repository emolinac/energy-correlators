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

void macro_print_jet_kinematics_data2mcreco()
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

        TH1F* h_jet_phi_data   = new TH1F("h_jet_phi_data"  ,"",10,-TMath::Pi(),TMath::Pi());
        TH1F* h_jet_phi_mc     = new TH1F("h_jet_phi_mc"    ,"",10,-TMath::Pi(),TMath::Pi());
        TH1F* h_jet_phi_mcreco = new TH1F("h_jet_phi_mcreco","",10,-TMath::Pi(),TMath::Pi());

        TH1F* h_z_pt_data   = new TH1F("h_z_pt_data"   ,"",20,0,200);
        TH1F* h_z_pt_mc     = new TH1F("h_z_pt_mc"     ,"",20,0,200);
        TH1F* h_z_pt_mcreco = new TH1F("h_z_pt_mcreco" ,"",20,0,200);
        
        TH1F* h_z_y_data   = new TH1F("h_z_y_data"  ,"",12,eta_min,eta_max);
        TH1F* h_z_y_mc     = new TH1F("h_z_y_mc"    ,"",12,eta_min,eta_max);
        TH1F* h_z_y_mcreco = new TH1F("h_z_y_mcreco","",12,eta_min,eta_max);

        TH1F* h_z_eta_data   = new TH1F("h_z_eta_data"  ,"",12,eta_min,eta_max);
        TH1F* h_z_eta_mc     = new TH1F("h_z_eta_mc"    ,"",12,eta_min,eta_max);
        TH1F* h_z_eta_mcreco = new TH1F("h_z_eta_mcreco","",12,eta_min,eta_max);

        TH1F* h_z_phi_data   = new TH1F("h_z_phi_data"  ,"",10,-TMath::Pi(),TMath::Pi());
        TH1F* h_z_phi_mc     = new TH1F("h_z_phi_mc"    ,"",10,-TMath::Pi(),TMath::Pi());
        TH1F* h_z_phi_mcreco = new TH1F("h_z_phi_mcreco","",10,-TMath::Pi(),TMath::Pi());

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
        
        h_jet_pt_data->Divide(h_jet_pt_mcreco);
        h_jet_eta_data->Divide(h_jet_eta_mcreco);
        h_jet_y_data->Divide(h_jet_y_mcreco);
        h_jet_phi_data->Divide(h_jet_phi_mcreco);
        
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

        h_z_pt_data->Divide(h_z_pt_mcreco);
        h_z_eta_data->Divide(h_z_eta_mcreco);
        h_z_y_data->Divide(h_z_y_mcreco);
        h_z_phi_data->Divide(h_z_phi_mcreco);
        
        TCanvas* c = new TCanvas("c","",800,600);
        c->Draw();

        THStack* hjet_pt  = new THStack();
        hjet_pt->Add(h_jet_pt_data);
        hjet_pt->SetTitle(";p_{T,jet}(GeV);Data/MC(Reco)");
        
        THStack* hjeteta = new THStack();
        hjeteta->Add(h_jet_eta_data);
        hjeteta->SetTitle(";#eta_{jet};Data/MC(Reco)");

        THStack* hjety = new THStack();
        hjety->Add(h_jet_y_data);
        hjety->SetTitle(";y_{jet};Data/MC(Reco)");

        THStack* hjetphi = new THStack();
        hjetphi->Add(h_jet_phi_data);
        hjetphi->SetTitle(";#phi_{jet};Data/MC(Reco)");

        THStack* hzpt  = new THStack();
        hzpt->Add(h_z_pt_data);
        hzpt->SetTitle(";p^{Z}_{T}(GeV);Data/MC(Reco)");
        
        THStack* hzeta = new THStack();
        hzeta->Add(h_z_eta_data);
        hzeta->SetTitle(";#eta_{Z};Data/MC(Reco)");

        THStack* hzy = new THStack();
        hzy->Add(h_z_y_data);
        hzy->SetTitle(";y_{Z};Data/MC(Reco)");

        THStack* hzphi = new THStack();
        hzphi->Add(h_z_phi_data);
        hzphi->SetTitle(";#phi_{Z};Data/MC(Reco)");

        hjet_pt->Draw("NOSTACK");
        hjet_pt->SetMaximum(1.3);
        hjet_pt->SetMinimum(0.7);
        c->Print("./plots/jet_pt_kinematics_data2mcreco_ratio.pdf");

        hzpt->Draw("NOSTACK");
        hzpt->SetMaximum(1.3);
        hzpt->SetMinimum(0.7);
        c->Print("./plots/z_pt_kinematics_data2mcreco_ratio.pdf");

        hjeteta->Draw("NOSTACK");
        hjeteta->SetMaximum(1.3);
        hjeteta->SetMinimum(0.7);
        c->Print("./plots/jet_eta_kinematics_data2mcreco_ratio.pdf");

        hzeta->Draw("NOSTACK");
        hzeta->SetMaximum(1.3);
        hzeta->SetMinimum(0.7);
        c->Print("./plots/z_eta_kinematics_data2mcreco_ratio.pdf");

        hjety->Draw("NOSTACK");
        hjety->SetMaximum(1.3);
        hjety->SetMinimum(0.7);
        c->Print("./plots/jet_y_kinematics_data2mcreco_ratio.pdf");

        hzy->Draw("NOSTACK");
        hzy->SetMaximum(1.3);
        hzy->SetMinimum(0.7);
        c->Print("./plots/z_y_kinematics_data2mcreco_ratio.pdf");

        hjetphi->Draw("NOSTACK");
        hjetphi->SetMaximum(1.3);
        hjetphi->SetMinimum(0.7);
        c->Print("./plots/jet_phi_kinematics_data2mcreco_ratio.pdf");

        hzphi->Draw("NOSTACK");
        hzphi->SetMaximum(1.3);
        hzphi->SetMinimum(0.7);
        c->Print("./plots/z_phi_kinematics_data2mcreco_ratio.pdf");
}