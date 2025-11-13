#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.cpp"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

void macro_print_systerrors_contribution(int niter = 4, int niter_jet = 4)
{
        TFile* fnominal = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niter%i_niterjet%i--get-nominal.root",niter,niter_jet)).c_str());

        // Define visual immediately
        TCanvas* c = new TCanvas("c", "", 1920, 1480);
        c->Draw();

        THStack* s[nbin_jet_pt];
        TLegend* l[nbin_jet_pt];
        THStack* s_eqch[nbin_jet_pt];
        TLegend* l_eqch[nbin_jet_pt];
        THStack* s_neqch[nbin_jet_pt];
        TLegend* l_neqch[nbin_jet_pt];

        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                l[i]       = new TLegend(1 - 0.32 - gPad->GetRightMargin(), 1 - 0.35 - gPad->GetTopMargin(),1 - gPad->GetRightMargin(), 1 - 0.03 - gPad->GetTopMargin());
                s[i]       = new THStack(Form("hs%i",i),Form("hs%i",i));
                l_eqch[i]  = new TLegend(1 - 0.32 - gPad->GetRightMargin(), 1 - 0.35 - gPad->GetTopMargin(),1 - gPad->GetRightMargin(), 1 - 0.03 - gPad->GetTopMargin());
                s_eqch[i]  = new THStack(Form("hseqch%i",i),Form("hseqch%i",i));
                l_neqch[i] = new TLegend(1 - 0.32 - gPad->GetRightMargin(), 1 - 0.35 - gPad->GetTopMargin(),1 - gPad->GetRightMargin(), 1 - 0.03 - gPad->GetTopMargin());
                s_neqch[i] = new THStack(Form("hsneqch%i",i),Form("hsneqch%i",i));
        }

        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_eec_syst[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge_syst[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge_syst[nbin_jet_pt]; 
        
        // Nominal operations
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_eec[bin]           = (TH1F*) fnominal->Get(Form("hcorr_eec%i",bin));
                hcorr_eec_eqcharge[bin]  = (TH1F*) fnominal->Get(Form("hcorr_eqcheec%i",bin));
                hcorr_eec_neqcharge[bin] = (TH1F*) fnominal->Get(Form("hcorr_neqcheec%i",bin));

                hcorr_eec_syst[bin]           = (TH1F*) hcorr_eec[bin]->Clone(Form("hcorr_eec_syst%i",bin));
                hcorr_eec_eqcharge_syst[bin]  = (TH1F*) hcorr_eec_eqcharge[bin]->Clone(Form("hcorr_eec_eqcharge_syst%i",bin));
                hcorr_eec_neqcharge_syst[bin] = (TH1F*) hcorr_eec_neqcharge[bin]->Clone(Form("hcorr_eec_neqcharge_syst%i",bin));
                
                set_histogram_style(hcorr_eec[bin]               , corr_marker_color_jet_pt[bin+3], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_eec_syst[bin]          , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_eec_eqcharge[bin]      , corr_marker_color_jet_pt[bin+3], std_line_width, std_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_eec_eqcharge_syst[bin] , corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_eec_neqcharge[bin]     , corr_marker_color_jet_pt[bin+3], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_eec_neqcharge_syst[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                
                hcorr_eec[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin+3], 0.3);
                hcorr_eec_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_eqcharge[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin+3], 0.3);
                hcorr_eec_eqcharge_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_neqcharge[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin+3], 0.3);
                hcorr_eec_neqcharge_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                
                l[bin]->SetHeader(Form("%.0f<p_{T,jet}(GeV)<%.0f",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l[bin]->AddEntry(hcorr_eec_syst[bin],"Systematic Error","f");
                l_eqch[bin]->SetHeader(Form("%.0f<p_{T,jet}(GeV)<%.0f",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l_eqch[bin]->AddEntry(hcorr_eec_eqcharge_syst[bin],"Systematic Error","f");
                l_neqch[bin]->SetHeader(Form("%.0f<p_{T,jet}(GeV)<%.0f",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l_neqch[bin]->AddEntry(hcorr_eec_neqcharge_syst[bin],"Systematic Error","f");
        }

        // Calculate the total systematic error
        const int nsyst = sizeof(available_systematics)/sizeof(available_systematics[0]);
        TFile* fsyst[nsyst];
        TH1F* hsyst_eec[nbin_jet_pt];
        TH1F* hsyst_eec_eqcharge[nbin_jet_pt];
        TH1F* hsyst_eec_neqcharge[nbin_jet_pt];

        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                if (gSystem->AccessPathName((output_folder + systematic_namef[available_systematics[syst_index]]).c_str()))
                        continue;
                
                fsyst[syst_index] = new TFile((output_folder + systematic_namef[available_systematics[syst_index]]).c_str());
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hsyst_eec[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_eec%i",bin));
                        set_histo_with_systematics(hsyst_eec[bin], hcorr_eec[bin], hcorr_eec_syst[bin], syst_index, false);
                        
                        hsyst_eec_eqcharge[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_eqcheec%i",bin));
                        set_histo_with_systematics(hsyst_eec_eqcharge[bin], hcorr_eec_eqcharge[bin], hcorr_eec_eqcharge_syst[bin], syst_index, false);
                        
                        hsyst_eec_neqcharge[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_neqcheec%i",bin));
                        set_histo_with_systematics(hsyst_eec_neqcharge[bin], hcorr_eec_neqcharge[bin], hcorr_eec_neqcharge_syst[bin], syst_index, false);
                }
        }

        // Extract the statistical error from the total error
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                substract_stat_error(hcorr_eec[bin]          , hcorr_eec_syst[bin]);
                substract_stat_error(hcorr_eec_eqcharge[bin] , hcorr_eec_eqcharge_syst[bin]);
                substract_stat_error(hcorr_eec_neqcharge[bin], hcorr_eec_neqcharge_syst[bin]);

                set_histoa_errors_as_histob_content(hcorr_eec_syst[bin]          ,hcorr_eec_syst[bin]);
                set_histoa_errors_as_histob_content(hcorr_eec_eqcharge_syst[bin] ,hcorr_eec_eqcharge_syst[bin]);
                set_histoa_errors_as_histob_content(hcorr_eec_neqcharge_syst[bin],hcorr_eec_neqcharge_syst[bin]);
        }

        // In order to get the systematic total in the background I have to insert that first and then the other systematics
        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                s[i]->Add(hcorr_eec_syst[i],"HIST");
                s_eqch[i]->Add(hcorr_eec_eqcharge_syst[i],"HIST");
                s_neqch[i]->Add(hcorr_eec_neqcharge_syst[i],"HIST");
        }

        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                if (gSystem->AccessPathName((output_folder + systematic_namef[available_systematics[syst_index]]).c_str()))
                        continue;
                
                fsyst[syst_index] = new TFile((output_folder + systematic_namef[available_systematics[syst_index]]).c_str());
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hsyst_eec[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_eec%i",bin));
                        set_histogram_style(hsyst_eec[bin], corr_marker_color_jet_pt[syst_index+3], std_line_width, corr_marker_style_jet_pt[syst_index], std_marker_size+2);
                        s[bin]->Add(hsyst_eec[bin],"HIST PL");
                        l[bin]->AddEntry(hsyst_eec[bin],systematic_name[available_systematics[syst_index]].c_str(),"p");

                        hsyst_eec_eqcharge[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_eqcheec%i",bin));
                        set_histogram_style(hsyst_eec_eqcharge[bin], corr_marker_color_jet_pt[syst_index+3], std_line_width, corr_marker_style_jet_pt[syst_index], std_marker_size+2);
                        s_eqch[bin]->Add(hsyst_eec_eqcharge[bin],"HIST PL");
                        l_eqch[bin]->AddEntry(hsyst_eec_eqcharge[bin],systematic_name[available_systematics[syst_index]].c_str(),"p");

                        hsyst_eec_neqcharge[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_neqcheec%i",bin));
                        set_histogram_style(hsyst_eec_neqcharge[bin], corr_marker_color_jet_pt[syst_index+3], std_line_width, corr_marker_style_jet_pt[syst_index], std_marker_size+2);
                        s_neqch[bin]->Add(hsyst_eec_neqcharge[bin],"HIST PL");
                        l_neqch[bin]->AddEntry(hsyst_eec_neqcharge[bin],systematic_name[available_systematics[syst_index]].c_str(),"p");
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

                s_eqch[i]->Draw("NOSTACK");
                s_eqch[i]->SetTitle(";R_{L};Relative Uncertainty");
                s_eqch[i]->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
                s_eqch[i]->SetMaximum(0.45);
                l_eqch[i]->Draw("SAME");
                gPad->SetLogx(1);
                gPad->SetLogy(0);
                
                c->Print(Form("./plots/correqchargedeec_systematics_contribution_jetpt%i.pdf",i));

                s_neqch[i]->Draw("NOSTACK");
                s_neqch[i]->SetTitle(";R_{L};Relative Uncertainty");
                s_neqch[i]->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
                s_neqch[i]->SetMaximum(0.45);
                l_neqch[i]->Draw("SAME");
                gPad->SetLogx(1);
                gPad->SetLogy(0);
                
                c->Print(Form("./plots/corrneqchargedeec_systematics_contribution_jetpt%i.pdf",i));
        }
}
