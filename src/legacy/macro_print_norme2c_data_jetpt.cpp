#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_norme2c_data_jetpt()
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple = (TNtuple*) f->Get((name_ntuple_data).c_str());
    
    // Determine binning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    TH1F* h[Nbin_jet_pt];
    THStack* s = new THStack();
    TLegend* l = new TLegend(.6,.2,.9,.4);
    
    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        h[jet_pt_bin] = new TH1F(Form("h[%i]",jet_pt_bin),"",Nbin_R_L,R_L_min,R_L_max);
        h[jet_pt_bin]->Sumw2();

        ntuple->Draw(Form("R_L>>h[%i]",jet_pt_bin),e2c_jetpt_cut[jet_pt_bin],"goff");
        
        h[jet_pt_bin]->Scale(1./h[jet_pt_bin]->Integral());
        set_histogram_style(h[jet_pt_bin] , corr_marker_color_jet_pt[jet_pt_bin] , std_line_width, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size);

        s->Add(h[jet_pt_bin]);
        l->AddEntry(h[jet_pt_bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lpf");    
    }

    // Create Canvas and draw in it
    TCanvas* c = new TCanvas("","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    s->Draw("NOSTACK");
    s->SetTitle(";R_{L};Norm. E2C");
    //s->SetMaximum(1);

    gPad->SetLogx(1);
    //gPad->SetLogy(1);

    l->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print("./plots/norme2c_rl_data_jetpt.pdf");
}