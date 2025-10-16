#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

const int    nbin_jet_eta = 3;
const double jet_eta_binning[] = {2.5, 3, 3.5, 4};

void macro_print_weightpt_resolution_median_jeteta(int jet_pt_bin = 0)
{
        TFile* f    = new TFile((output_folder + namef_ntuple_eec_paircorrections_ct).c_str());

        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());
        float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth, jet_eta_reco, jet_eta_truth;
        set_unfolding_ntuple_branches(ntuple, &R_L_reco, &R_L_truth, &jet_pt_reco, &jet_pt_truth, &weight_pt_reco, &weight_pt_truth);
        ntuple->SetBranchAddress("jet_eta",       &jet_eta_reco);
        ntuple->SetBranchAddress("jet_eta_truth", &jet_eta_truth);
        
        TH1F* h_weightpt_resolution[nbin_jet_eta];
        TH1F* h_weightpt_resolution_cumulative[nbin_jet_eta];

        for (int bin = 0 ; bin < nbin_jet_eta ; bin++) {
                h_weightpt_resolution[bin]            = new TH1F(Form("h_weightpt_resolution%i",bin)           , "", 1000,0,2.5);
                h_weightpt_resolution_cumulative[bin] = new TH1F(Form("h_weightpt_resolution_cumulative%i",bin), "", 1000,0,2.5);
                h_weightpt_resolution[bin]->Sumw2();
                h_weightpt_resolution_cumulative[bin]->Sumw2();
                set_histogram_style(h_weightpt_resolution[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size-1);
                set_histogram_style(h_weightpt_resolution_cumulative[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size-1);
        }
        
        for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++) {
                ntuple->GetEntry(evt);

                if (weight_pt_truth != -999) {
                        for (int jet_eta_bin = 0 ; jet_eta_bin < nbin_jet_eta ; jet_eta_bin++) {
                                if ((jet_eta_reco > jet_eta_binning[jet_eta_bin] && jet_eta_reco < jet_eta_binning[jet_eta_bin + 1]) && 
                                    (jet_pt_reco > jet_pt_binning[jet_pt_bin] && jet_pt_reco < jet_pt_binning[jet_pt_bin + 1]))
                                        h_weightpt_resolution[jet_eta_bin]->Fill(weight_pt_reco/weight_pt_truth);        
                        }
                }
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        THStack* s = new THStack();
        TLegend* l = new TLegend(1-0.31-gPad->GetRightMargin(),gPad->GetBottomMargin()+0.01,0.99-gPad->GetRightMargin(),gPad->GetBottomMargin()+0.21);
        
        for (int bin = 0 ; bin < nbin_jet_eta ; bin++) {
                h_weightpt_resolution_cumulative[bin] = (TH1F*) h_weightpt_resolution[bin]->GetCumulative();
                s->Add(h_weightpt_resolution_cumulative[bin],"AP E1");
                l->AddEntry(h_weightpt_resolution_cumulative[bin],
                            Form("%.1f<#eta_{jet}<%.1f : Median = %.4f",jet_eta_binning[bin],jet_eta_binning[bin + 1],
                            get_median_from_cumulative(h_weightpt_resolution_cumulative[bin])),"lpf");
        }
        
        s->Draw("NOSTACK");
        s->SetTitle(";w_{reco}/w_{truth};Cumulative distribution");
        l->Draw("SAME");

        c->Print(Form("./plots/weightpt_resolution_cumulative_median_jet_eta_jetpt%i.pdf",jet_pt_bin));

        s = new THStack();
        l = new TLegend(1-0.31-gPad->GetRightMargin(),1-gPad->GetTopMargin()-0.21,0.99-gPad->GetRightMargin(),0.99-gPad->GetTopMargin());
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s->Add(h_weightpt_resolution[bin],"E");
                l->AddEntry(h_weightpt_resolution[bin],
                            Form("%.1f<#eta_{jet}<%.1f: Median = %.4f",jet_eta_binning[bin],jet_eta_binning[bin + 1],
                            get_median_from_cumulative(h_weightpt_resolution_cumulative[bin])),"lf");
        }
        
        s->Draw("NOSTACK");
        s->SetTitle(";w_{reco}/w_{truth};");
        l->Draw("SAME");

        c->Print(Form("./plots/weightpt_resolution_median_jet_eta_jetpt%i.pdf",jet_pt_bin));
}
