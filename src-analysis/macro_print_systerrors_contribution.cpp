#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_systerrors_contribution(int niter = 4, int niter_jet = 4)
{
        TFile* fnominal = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niter%i_niterjet%i--get-nominal.root",niter,niter_jet)).c_str());

        // Define visual immediately
        TCanvas* c = new TCanvas("c", "", 1920, 1480);
        c->Draw();

        THStack* s[nbin_jet_pt];
        TLegend* l[nbin_jet_pt];

        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                l[i] = new TLegend(1 - 0.32 - gPad->GetRightMargin(), 1 - 0.35 - gPad->GetTopMargin(),1 - gPad->GetRightMargin(), 1 - 0.03 - gPad->GetTopMargin());
                s[i] = new THStack(Form("hs%i",i),Form("hs%i",i));
        }

        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_eec_syst[nbin_jet_pt]; 
        
        // Nominal operations
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_eec[bin]      = (TH1F*) fnominal->Get(Form("hcorr_eec%i",bin));
                hcorr_eec_syst[bin] = new TH1F(Form("hcorr_eec_syst%i",bin),"", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                
                set_nominal_error_histo(hcorr_eec[bin], hcorr_eec_syst[bin]);
                
                set_histogram_style(hcorr_eec[bin]               , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_eec_syst[bin]          , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                
                hcorr_eec[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                
                s[bin]->Add(hcorr_eec_syst[bin]);
                
                l[bin]->SetHeader(Form("%.0f<p_{T,jet}(GeV)<%.0f",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l[bin]->AddEntry(hcorr_eec_syst[bin],"Statistical Error","f");
        }

        // Include the systematics in the whole deal
        const int nsyst = sizeof(available_systematics)/sizeof(available_systematics[0]);
        TFile* fsyst[nsyst];
        TH1F* hsyst_eec[nbin_jet_pt];

        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                if (gSystem->AccessPathName((output_folder + systematic_namef[available_systematics[syst_index]]).c_str()))
                        continue;
                
                fsyst[syst_index] = new TFile((output_folder + systematic_namef[available_systematics[syst_index]]).c_str());
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hsyst_eec[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_eec%i",bin));
                        set_histogram_style(hsyst_eec[bin], corr_marker_color_jet_pt[syst_index], std_line_width, corr_marker_style_jet_pt[syst_index], std_marker_size+2);
                        s[bin]->Add(hsyst_eec[bin],"HIST PL");
                        l[bin]->AddEntry(hsyst_eec[bin],systematic_name[available_systematics[syst_index]].c_str(),"p");
                }
        }
        
        // Print EECS
        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                s[i]->Draw("NOSTACK");
                s[i]->SetTitle(";R_{L};Relative Uncertainty");
                s[i]->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
                s[i]->SetMaximum(0.45);
                l[i]->Draw("SAME");
                gPad->SetLogx(1);
                gPad->SetLogy(0);
                
                c->Print(Form("./plots/correec_systematics_contribution_jetpt%i.pdf",i));
        }
}
