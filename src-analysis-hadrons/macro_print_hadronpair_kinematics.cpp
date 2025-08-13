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

        // TH1F* h_eec_data   = new TH1F("h_eec_data"   ,"",30,rl_min,0.5);
        // TH1F* h_eec_mc     = new TH1F("h_eec_mc"     ,"",30,rl_min,0.5);
        // TH1F* h_eec_mcreco = new TH1F("h_eec_mcreco" ,"",30,rl_min,0.5);
        
        TH1F* h_eec_data   = new TH1F("h_eec_data"   ,"",nbin_rl_nominal,rl_nominal_binning);
        TH1F* h_eec_mc     = new TH1F("h_eec_mc"     ,"",nbin_rl_nominal,rl_nominal_binning);
        TH1F* h_eec_mcreco = new TH1F("h_eec_mcreco" ,"",nbin_rl_nominal,rl_nominal_binning);
        
        TH1F* h_rl_data   = new TH1F("h_rl_data"   ,"",50,rl_absmin,rl_absmax);
        TH1F* h_rl_mc     = new TH1F("h_rl_mc"     ,"",50,rl_absmin,rl_absmax);
        TH1F* h_rl_mcreco = new TH1F("h_rl_mcreco" ,"",50,rl_absmin,rl_absmax);
        
        TH1F* h_weight_data   = new TH1F("h_weight_data"   ,"",50,weight_min,weight_max);
        TH1F* h_weight_mc     = new TH1F("h_weight_mc"     ,"",50,weight_min,weight_max);
        TH1F* h_weight_mcreco = new TH1F("h_weight_mcreco" ,"",50,weight_min,weight_max);

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
        
        // h_eec_data->Scale(1./h_eec_data->Integral());
        // h_eec_mc->Scale(1./h_eec_mc->Integral());
        // h_eec_mcreco->Scale(1./h_eec_mcreco->Integral());

        h_eec_data->Scale(1./h_eec_data->Integral(),"width");
        h_eec_mc->Scale(1./h_eec_mc->Integral(),"width");
        h_eec_mcreco->Scale(1./h_eec_mcreco->Integral(),"width");

        h_rl_data->Scale(1./h_rl_data->Integral());
        h_rl_mc->Scale(1./h_rl_mc->Integral());
        h_rl_mcreco->Scale(1./h_rl_mcreco->Integral());

        h_weight_data->Scale(1./h_weight_data->Integral());
        h_weight_mc->Scale(1./h_weight_mc->Integral());
        h_weight_mcreco->Scale(1./h_weight_mcreco->Integral());

        s_eec->Add(h_eec_data);
        s_eec->Add(h_eec_mc);
        s_eec->Add(h_eec_mcreco);
        s_rl->Add(h_rl_data);
        s_rl->Add(h_rl_mc);
        s_rl->Add(h_rl_mcreco);
        s_weight->Add(h_weight_data);
        s_weight->Add(h_weight_mc);
        s_weight->Add(h_weight_mcreco);
        
        l_eec->AddEntry(h_eec_data,"data","lpf");
        l_eec->AddEntry(h_eec_mc,"mc","lpf");
        l_eec->AddEntry(h_eec_mcreco,"mcreco","lpf");
        l_rl->AddEntry(h_rl_data,"data","lpf");
        l_rl->AddEntry(h_rl_mc,"mc","lpf");
        l_rl->AddEntry(h_rl_mcreco,"mcreco","lpf");
        l_weight->AddEntry(h_weight_data,"data","lpf");
        l_weight->AddEntry(h_weight_mc,"mc","lpf");
        l_weight->AddEntry(h_weight_mcreco,"mcreco","lpf");
        
        s_eec->Draw("NOSTACK");
        s_eec->SetTitle(";R_{L};Normalized EEC Distributions");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        l_eec->Draw("SAME");

        c->Print("./plots/kinematics_eec.pdf");

        s_rl->Draw("NOSTACK");
        s_rl->SetTitle(";R_{L};Normalized Distributions");
        gPad->SetLogx(0);
        gPad->SetLogy(1);
        l_rl->Draw("SAME");

        c->Print("./plots/kinematics_rl.pdf");

        s_weight->Draw("NOSTACK");
        s_weight->SetTitle(";w;Normalized Distributions");
        gPad->SetLogx(0);
        gPad->SetLogy(1);
        l_weight->Draw("SAME"); 

        c->Print("./plots/kinematics_weights.pdf");
}