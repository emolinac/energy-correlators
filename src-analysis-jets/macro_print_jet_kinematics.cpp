#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_jet_kinematics()
{
        TFile* fin     = new TFile((output_folder + namef_ntuple_mc_e2c).c_str());
        TFile* findata = new TFile((output_folder + namef_ntuple_e2c_corr).c_str());
        
        TNtuple* ntuple_jet_data   = (TNtuple*) findata->Get(name_ntuple_corrjet.c_str());
        TNtuple* ntuple_jet_mc     = (TNtuple*) fin->Get(name_ntuple_mc_jet.c_str());
        TNtuple* ntuple_jet_mcreco = (TNtuple*) fin->Get(name_ntuple_mcreco_jet.c_str());

        TH1F* h_pt_data   = new TH1F("h_pt_data"   ,"",20,jet_pt_binning[0],jet_pt_binning[nbin_jet_pt]);
        TH1F* h_pt_mc     = new TH1F("h_pt_mc"     ,"",20,jet_pt_binning[0],jet_pt_binning[nbin_jet_pt]);
        TH1F* h_pt_mcreco = new TH1F("h_pt_mcreco" ,"",20,jet_pt_binning[0],jet_pt_binning[nbin_jet_pt]);
        
        TH1F* h_z_pt_data   = new TH1F("h_z_pt_data"   ,"",20,0,200);
        TH1F* h_z_pt_mc     = new TH1F("h_z_pt_mc"     ,"",20,0,200);
        TH1F* h_z_pt_mcreco = new TH1F("h_z_pt_mcreco" ,"",20,0,200);
        
        TH1F* h_eta_data   = new TH1F("h_eta_data"  ,"",15,jet_eta_min,jet_eta_max);
        TH1F* h_eta_mc     = new TH1F("h_eta_mc"    ,"",15,jet_eta_min,jet_eta_max);
        TH1F* h_eta_mcreco = new TH1F("h_eta_mcreco","",15,jet_eta_min,jet_eta_max);

        TH1F* h_z_y_data   = new TH1F("h_z_y_data"  ,"",15,eta_min,eta_max);
        TH1F* h_z_y_mc     = new TH1F("h_z_y_mc"    ,"",15,eta_min,eta_max);
        TH1F* h_z_y_mcreco = new TH1F("h_z_y_mcreco","",15,eta_min,eta_max);

        set_histogram_style(h_pt_data   , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_pt_mc     , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_pt_mcreco , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eta_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eta_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eta_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);

        set_histogram_style(h_z_pt_data   , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_pt_mc     , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_pt_mcreco , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_y_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_y_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_z_y_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);

        ntuple_jet_data->Project("h_pt_data", "jet_pt","");
        ntuple_jet_mc->Project("h_pt_mc", "jet_pt","");
        ntuple_jet_mcreco->Project("h_pt_mcreco", "jet_pt","");
        
        ntuple_jet_data->Project("h_eta_data","jet_eta","");
        ntuple_jet_mc->Project("h_eta_mc","jet_eta","");
        ntuple_jet_mcreco->Project("h_eta_mcreco","jet_eta","");
        
        ntuple_jet_data->Project("h_z_pt_data","z_pt","");
        ntuple_jet_mc->Project("h_z_pt_mc","z_pt","");
        ntuple_jet_mcreco->Project("h_z_pt_mcreco","z_pt","");
        
        ntuple_jet_data->Project("h_z_y_data","z_y","");
        ntuple_jet_mc->Project("h_z_y_mc","z_y","");
        ntuple_jet_mcreco->Project("h_z_y_mcreco","z_y","");

        h_pt_data->Scale(1./h_pt_data->Integral());
        h_pt_mc->Scale(1./h_pt_mc->Integral());
        h_pt_mcreco->Scale(1./h_pt_mcreco->Integral());
        h_eta_data->Scale(1./h_eta_data->Integral());
        h_eta_mc->Scale(1./h_eta_mc->Integral());
        h_eta_mcreco->Scale(1./h_eta_mcreco->Integral());
        h_z_pt_data->Scale(1./h_z_pt_data->Integral());
        h_z_pt_mc->Scale(1./h_z_pt_mc->Integral());
        h_z_pt_mcreco->Scale(1./h_z_pt_mcreco->Integral());
        h_z_y_data->Scale(1./h_z_y_data->Integral());
        h_z_y_mc->Scale(1./h_z_y_mc->Integral());
        h_z_y_mcreco->Scale(1./h_z_y_mcreco->Integral());

        TCanvas* c = new TCanvas("c","",800,600);
        c->Draw();

        THStack* hjet_pt  = new THStack();
        hjet_pt->Add(h_pt_data);
        hjet_pt->Add(h_pt_mc);
        hjet_pt->Add(h_pt_mcreco);
        hjet_pt->SetTitle(";p^{jet}_{t}(GeV);Normalized Distributions");
        
        THStack* hjeteta = new THStack();
        hjeteta->Add(h_eta_data);
        hjeteta->Add(h_eta_mc);
        hjeteta->Add(h_eta_mcreco);
        hjeteta->SetTitle(";#eta^{jet};Normalized Distributions");

        THStack* hzpt  = new THStack();
        hzpt->Add(h_z_pt_data);
        hzpt->Add(h_z_pt_mc);
        hzpt->Add(h_z_pt_mcreco);
        hzpt->SetTitle(";p^{Z}_{t}(GeV);Normalized Distributions");
        
        THStack* hzeta = new THStack();
        hzeta->Add(h_z_y_data);
        hzeta->Add(h_z_y_mc);
        hzeta->Add(h_z_y_mcreco);
        hzeta->SetTitle(";y^{Z};Normalized Distributions");

        TLegend* ljet_pt = new TLegend();
        ljet_pt->AddEntry(h_pt_data,"data","lpf");
        ljet_pt->AddEntry(h_pt_mc,"mc","lpf");
        ljet_pt->AddEntry(h_pt_mcreco,"mcreco","lpf");

        TLegend* ljeteta = new TLegend();
        ljeteta->AddEntry(h_eta_data,"data","lpf");
        ljeteta->AddEntry(h_eta_mc,"mc","lpf");
        ljeteta->AddEntry(h_eta_mcreco,"mcreco","lpf");

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
        c->Print("./plots/z_y_kinematics.pdf");

}