#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_ndtrs()
{
        // Open the necessary files
        TFile* f    = new TFile((output_folder + namef_ntuple_eec_paircorr).c_str());
        TFile* fsim = new TFile((output_folder + namef_ntuple_eec_paircorrections).c_str());
        
        // Get the corresponding Ntuples
        TNtuple* ntuple      = (TNtuple*) f->Get((name_ntuple_corrjet).c_str());
        
        TNtuple* ntuple_reco = (TNtuple*) fsim->Get((name_ntuple_mcreco_jet).c_str());
        TNtuple* ntuple_mc   = (TNtuple*) fsim->Get((name_ntuple_mc_jet).c_str());

        TH1F* h[nbin_jet_pt];
        THStack* hs = new THStack();
        TLegend* l = new TLegend();

        // Data
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                h[bin] = new TH1F(Form("h[%i]",bin),"",25,0,25);

                ntuple->Project(Form("h[%i]",bin),"jet_ndtr",pair_jet_pt_cut[bin]);
                
                hs->Add(h[bin],"E X0");
                l->AddEntry(h[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)  <Dtr> = %.0f",jet_pt_binning[bin],jet_pt_binning[bin + 1],h[bin]->GetMean()),"lpf");
                
                set_histogram_style(h[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();
        hs->Draw("NOSTACK");
        hs->SetTitle("Data;#Dtrs per Jet;");
        l->Draw("SAME");
        gPad->SetLogy(1);
        c->Print("./plots/jet_ndtrs.pdf");

        // MCReco
        hs = new THStack();
        l = new TLegend();
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                // h[bin] = new TH1F(Form("h[%i]",bin),"",25,0,25);
                h[bin]->Reset();

                ntuple_reco->Project(Form("h[%i]",bin),"jet_ndtr",pair_jet_pt_cut[bin]);
                
                hs->Add(h[bin],"E X0");
                l->AddEntry(h[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)  <Dtr> = %.0f",jet_pt_binning[bin],jet_pt_binning[bin + 1],h[bin]->GetMean()),"lpf");
                
                set_histogram_style(h[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        }

        hs->Draw("NOSTACK");
        hs->SetTitle("MCReco;#Dtrs per Jet;");
        l->Draw("SAME");
        gPad->SetLogy(1);
        c->Print("./plots/jet_ndtrs_mcreco.pdf");

        // MC
        hs = new THStack();
        l = new TLegend();
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                // h[bin] = new TH1F(Form("h[%i]",bin),"",25,0,25);
                h[bin]->Reset();

                ntuple_mc->Project(Form("h[%i]",bin),"jet_ndtr",pair_jet_pt_cut[bin]);
                
                hs->Add(h[bin],"E X0");
                l->AddEntry(h[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)  <Dtr> = %.0f",jet_pt_binning[bin],jet_pt_binning[bin + 1],h[bin]->GetMean()),"lpf");
                
                set_histogram_style(h[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        }

        hs->Draw("NOSTACK");
        hs->SetTitle("MC;#Dtrs per Jet;");
        l->Draw("SAME");
        gPad->SetLogy(1);
        c->Print("./plots/jet_ndtrs_mc.pdf");
}