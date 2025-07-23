#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_weightpt_rl_jetpt()
{
    TFile* fin = new TFile((output_folder + namef_ntuple_e2c_corr).c_str());

    TNtuple* ntuple = (TNtuple*) fin->Get(name_ntuple_data.c_str());

    TH2D* h[Nbin_jet_pt];

    TCanvas* c = new TCanvas("c","",1500,500);
    c->Draw();
    c->Divide(3,1);
    
    for (int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        h[bin] = new TH2D(Form("h%i",bin),"",300,0.01,1.2,70000,10E-7,.3);
        h[bin]->SetTitle(Form("%.0f<p^{jet}_{t}<%.0f;R_{L};w",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
        ntuple->Project(Form("h%i",bin),"weight_pt:R_L",pair_jetpt_cut[bin]);

        c->cd(bin + 1);
        h[bin]->Draw("col");
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        h[bin]->Smooth();
    }

    c->Print("./plots/weightpt_rl_jetpt_distributions.pdf");
}