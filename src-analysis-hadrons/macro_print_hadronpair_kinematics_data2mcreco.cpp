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

void macro_print_hadronpair_kinematics_data2mcreco()
{
        TFile* fin = new TFile((output_folder + namef_ntuple_hadron_jet).c_str());
        
        TNtuple* ntuple_data   = (TNtuple*) fin->Get(name_ntuple_data_pair.c_str());
        TNtuple* ntuple_mc     = (TNtuple*) fin->Get(name_ntuple_mc_pair.c_str());
        TNtuple* ntuple_mcreco = (TNtuple*) fin->Get(name_ntuple_mcreco_pair.c_str());

        TH1F* h_eec_data   = new TH1F("h_eec_data"   ,"",nbin_rl_nominal,rl_nominal_binning);
        TH1F* h_eec_mc     = new TH1F("h_eec_mc"     ,"",nbin_rl_nominal,rl_nominal_binning);
        TH1F* h_eec_mcreco = new TH1F("h_eec_mcreco" ,"",nbin_rl_nominal,rl_nominal_binning);
        
        TH1F* h_rl_data   = new TH1F("h_rl_data"   ,"",nbin_rl_nominal,rl_nominal_binning);
        TH1F* h_rl_mc     = new TH1F("h_rl_mc"     ,"",nbin_rl_nominal,rl_nominal_binning);
        TH1F* h_rl_mcreco = new TH1F("h_rl_mcreco" ,"",nbin_rl_nominal,rl_nominal_binning);
        
        TH1F* h_weight_data   = new TH1F("h_weight_data"   ,"",nbin_weight,weight_binning);
        TH1F* h_weight_mc     = new TH1F("h_weight_mc"     ,"",nbin_weight,weight_binning);
        TH1F* h_weight_mcreco = new TH1F("h_weight_mcreco" ,"",nbin_weight,weight_binning);

        set_histogram_style(h_eec_data     , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eec_mc       , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eec_mcreco   , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_rl_data      , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_rl_mc        , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_rl_mcreco    , 868 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_weight_data  , 875 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_weight_mc    , 797 , std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_weight_mcreco, 868 , std_line_width, std_marker_style, std_marker_size);

        TCanvas* c = new TCanvas("c","",800,600);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        THStack* s_eec    = new THStack();
        TLegend* l_eec    = new TLegend();
        THStack* s_rl     = new THStack();
        TLegend* l_rl     = new TLegend();
        THStack* s_weight = new THStack();
        TLegend* l_weight = new TLegend();

        ntuple_data->Project("h_eec_data","R_L","weight_pt");
        ntuple_mc->Project("h_eec_mc","R_L","weight_pt");
        ntuple_mcreco->Project("h_eec_mcreco","R_L","weight_pt");
        
        ntuple_data->Project("h_rl_data","R_L");
        ntuple_mc->Project("h_rl_mc","R_L");
        ntuple_mcreco->Project("h_rl_mcreco","R_L");
        
        ntuple_data->Project("h_weight_data","weight_pt");
        ntuple_mc->Project("h_weight_mc","weight_pt");
        ntuple_mcreco->Project("h_weight_mcreco","weight_pt");
        
        h_eec_data->Scale(1./h_eec_data->Integral(),"width");
        h_eec_mc->Scale(1./h_eec_mc->Integral(),"width");
        h_eec_mcreco->Scale(1./h_eec_mcreco->Integral(),"width");

        h_eec_data->Divide(h_eec_mcreco);

        h_rl_data->Scale(1./h_rl_data->Integral());
        h_rl_mc->Scale(1./h_rl_mc->Integral());
        h_rl_mcreco->Scale(1./h_rl_mcreco->Integral());

        h_rl_data->Divide(h_rl_mcreco);

        h_weight_data->Scale(1./h_weight_data->Integral());
        h_weight_mc->Scale(1./h_weight_mc->Integral());
        h_weight_mcreco->Scale(1./h_weight_mcreco->Integral());

        h_weight_data->Divide(h_weight_mcreco);

        s_eec->Add(h_eec_data);
        s_rl->Add(h_rl_data);
        s_weight->Add(h_weight_data);
        
        s_eec->Draw("NOSTACK");
        s_eec->SetTitle(";R_{L};Data/MC(Reco)");
        s_eec->SetMaximum(1.3);
        s_eec->SetMinimum(0.7);
        gPad->SetLogx(1);
        
        c->Print("./plots/kinematics_eec_data2mcreco_ratio.pdf");

        s_rl->Draw("NOSTACK");
        s_rl->SetTitle(";R_{L};Data/MC(Reco)");
        s_rl->SetMaximum(1.3);
        s_rl->SetMinimum(0.7);
        gPad->SetLogx(1);
        
        c->Print("./plots/kinematics_rl_data2mcreco_ratio.pdf");

        s_weight->Draw("NOSTACK");
        s_weight->SetTitle(";w;Data/MC(Reco)");
        s_weight->SetMaximum(1.3);
        s_weight->SetMinimum(0.7);
        gPad->SetLogx(1);
        
        c->Print("./plots/kinematics_weights_data2mcreco_ratio.pdf");
}