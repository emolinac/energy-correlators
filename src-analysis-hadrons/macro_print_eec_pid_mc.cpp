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

void macro_print_eec_pid_mc()
{
        TFile* fin = new TFile((output_folder + namef_ntuple_mc_eec).c_str());
        
        TNtuple* ntuple_mc     = (TNtuple*) fin->Get(name_ntuple_mc.c_str());
        TNtuple* ntuple_mcreco = (TNtuple*) fin->Get(name_ntuple_mcreco.c_str());

        TH1F* h_eec_mc_pi = new TH1F("h_eec_mc_pi","",nbin_rl_nominal,rl_nominal_binning);
        TH1F* h_eec_mc_k  = new TH1F("h_eec_mc_k" ,"",nbin_rl_nominal,rl_nominal_binning);
        TH1F* h_eec_mc_p  = new TH1F("h_eec_mc_p" ,"",nbin_rl_nominal,rl_nominal_binning);
        
        TH1F* h_w_mc_pi = new TH1F("h_w_mc_pi","",nbin_weight,weight_binning);
        TH1F* h_w_mc_k  = new TH1F("h_w_mc_k" ,"",nbin_weight,weight_binning);
        TH1F* h_w_mc_p  = new TH1F("h_w_mc_p" ,"",nbin_weight,weight_binning);
        
        set_histogram_style(h_eec_mc_pi, 797 ,std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eec_mc_k , 868 ,std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_eec_mc_p , 875 ,std_line_width, std_marker_style, std_marker_size);
        
        set_histogram_style(h_w_mc_pi, 797 ,std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_w_mc_k , 868 ,std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_w_mc_p , 875 ,std_line_width, std_marker_style, std_marker_size);
        
        TCanvas* c = new TCanvas("c","",800,600);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        THStack* s_eec = new THStack();
        THStack* s_w   = new THStack();
        TLegend* l_eec = new TLegend();
        TLegend* l_w   = new TLegend();
        
        ntuple_mc->Project("h_eec_mc_pi","R_L","weight_pt*(TMath::Abs(h1_pid)==211||TMath::Abs(h2_pid)==211)");
        ntuple_mc->Project("h_eec_mc_k" ,"R_L","weight_pt*(TMath::Abs(h1_pid)==321||TMath::Abs(h2_pid)==321)");
        ntuple_mc->Project("h_eec_mc_p" ,"R_L","weight_pt*(TMath::Abs(h1_pid)==2212||TMath::Abs(h2_pid)==2212)");
        
        ntuple_mc->Project("h_w_mc_pi","weight_pt","(TMath::Abs(h1_pid)==211||TMath::Abs(h2_pid)==211)");
        ntuple_mc->Project("h_w_mc_k" ,"weight_pt","(TMath::Abs(h1_pid)==321||TMath::Abs(h2_pid)==321)");
        ntuple_mc->Project("h_w_mc_p" ,"weight_pt","(TMath::Abs(h1_pid)==2212||TMath::Abs(h2_pid)==2212)");
        
        h_eec_mc_pi->Scale(1./h_eec_mc_pi->Integral(),"width");
        h_eec_mc_k->Scale(1./h_eec_mc_k->Integral(),"width");
        h_eec_mc_p->Scale(1./h_eec_mc_p->Integral(),"width");

        h_w_mc_pi->Scale(1./h_w_mc_pi->Integral(),"width");
        h_w_mc_k->Scale(1./h_w_mc_k->Integral(),"width");
        h_w_mc_p->Scale(1./h_w_mc_p->Integral(),"width");

        s_eec->Add(h_eec_mc_pi);
        s_eec->Add(h_eec_mc_k);
        s_eec->Add(h_eec_mc_p);
        
        s_w->Add(h_w_mc_pi);
        s_w->Add(h_w_mc_k);
        s_w->Add(h_w_mc_p);
        
        l_eec->AddEntry(h_eec_mc_pi,"#pi pairs","lpf");
        l_eec->AddEntry(h_eec_mc_k ,"K pairs","lpf");
        l_eec->AddEntry(h_eec_mc_p ,"p pairs","lpf");
        
        l_w->AddEntry(h_w_mc_pi,"#pi pairs","lpf");
        l_w->AddEntry(h_w_mc_k ,"K pairs","lpf");
        l_w->AddEntry(h_w_mc_p ,"p pairs","lpf");
        
        s_eec->Draw("NOSTACK");
        s_eec->SetTitle(";R_{L};Normalized EEC Distributions");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        l_eec->Draw("SAME");

        c->Print("./plots/eec_pid_mc.pdf");

        s_w->Draw("NOSTACK");
        s_w->SetTitle(";w_{p_{T}};Normalized Pair Distributions");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        l_w->Draw("SAME");

        c->Print("./plots/weight_pid_mc.pdf");
}