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

void macro_print_measurements_with_systematics(int niter = 4, int niter_jet = 4)
{
        TFile* fnominal = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niter%i_niterjet%i--get-nominal.root",niter,niter_jet)).c_str());

        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_eec_syst[nbin_jet_pt]; 
        TH1F* hcorr_tau[nbin_jet_pt]; 
        TH1F* hcorr_tau_syst[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge_syst[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge_syst[nbin_jet_pt]; 
        
        // Fill the nominal histograms
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hcorr_eec[bin]      = (TH1F*) fnominal->Get(Form("hcorr_eec%i",bin));
                hcorr_eec_syst[bin] = (TH1F*) hcorr_eec[bin]->Clone(Form("hcorr_eec_syst%i",bin));

                hcorr_tau[bin]      = (TH1F*) fnominal->Get(Form("hcorr_tau%i",bin));
                hcorr_tau_syst[bin] = (TH1F*) hcorr_tau[bin]->Clone(Form("hcorr_tau_syst%i",bin));

                hcorr_eec_eqcharge[bin]      = (TH1F*) fnominal->Get(Form("hcorr_eqcheec%i",bin));
                hcorr_eec_eqcharge_syst[bin] = (TH1F*) hcorr_eec_eqcharge[bin]->Clone(Form("hcorr_eqcheec_syst%i",bin));
                
                hcorr_eec_neqcharge[bin]      = (TH1F*) fnominal->Get(Form("hcorr_neqcheec%i",bin));
                hcorr_eec_neqcharge_syst[bin] = (TH1F*) hcorr_eec_neqcharge[bin]->Clone(Form("hcorr_neqcheec_syst%i",bin));
                
                set_histogram_style(hcorr_eec[bin]               , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_eec_syst[bin]          , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_tau[bin]               , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_tau_syst[bin]          , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_eec_eqcharge[bin]      , corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size+2);
                set_histogram_style(hcorr_eec_eqcharge_syst[bin] , corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size+2);
                set_histogram_style(hcorr_eec_neqcharge[bin]     , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hcorr_eec_neqcharge_syst[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                
                hcorr_eec[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_tau[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_tau_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_eqcharge[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_eqcharge_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_neqcharge[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_neqcharge_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
        }

        // Include the systematics in the whole deal
        const int nsyst = sizeof(available_systematics)/sizeof(available_systematics[0]);
        TFile* fsyst[nsyst];
        TH1F* hsyst_eec[nbin_jet_pt];
        TH1F* hsyst_tau[nbin_jet_pt];
        TH1F* hsyst_eec_eqcharge[nbin_jet_pt];
        TH1F* hsyst_eec_neqcharge[nbin_jet_pt];

        std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                if (gSystem->AccessPathName((output_folder + systematic_namef[available_systematics[syst_index]]).c_str()))
                        continue;
                
                fsyst[syst_index] = new TFile((output_folder + systematic_namef[available_systematics[syst_index]]).c_str());
                
                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hsyst_eec[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_eec%i",bin));

                        set_histo_with_systematics(hsyst_eec[bin], hcorr_eec[bin], hcorr_eec_syst[bin], syst_index);

                        if (bin != nbin_jet_pt-1) 
                                std::cout<<" & ";
                        else
                                std::cout<<" \\\\ ";
                }
                
                std::cout<<std::endl;

                delete fsyst[syst_index];
        }

        std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                if (gSystem->AccessPathName((output_folder + systematic_namef[available_systematics[syst_index]]).c_str()))
                        continue;
                
                fsyst[syst_index] = new TFile((output_folder + systematic_namef[available_systematics[syst_index]]).c_str());
                
                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hsyst_tau[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_tau%i",bin));

                        set_histo_with_systematics(hsyst_tau[bin], hcorr_tau[bin], hcorr_tau_syst[bin], syst_index);

                        if (bin != nbin_jet_pt-1) 
                                std::cout<<" & ";
                        else
                                std::cout<<" \\\\ ";
                }
                
                std::cout<<std::endl;

                delete fsyst[syst_index];
        }

        std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                if (gSystem->AccessPathName((output_folder + systematic_namef[available_systematics[syst_index]]).c_str()))
                        continue;
                
                fsyst[syst_index] = new TFile((output_folder + systematic_namef[available_systematics[syst_index]]).c_str());
                
                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";

                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hsyst_eec_eqcharge[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_eqcheec%i",bin));

                        set_histo_with_systematics(hsyst_eec_eqcharge[bin], hcorr_eec_eqcharge[bin], hcorr_eec_eqcharge_syst[bin], syst_index);

                        if (bin != nbin_jet_pt-1) 
                                std::cout<<" & ";
                        else
                                std::cout<<" \\\\ ";
                }
                
                std::cout<<std::endl;

                delete fsyst[syst_index];
        }

        std::cout<<"Source & $20<p_{T,jet}<30$ & $30<p_{T,jet}<50$ & $50<p_{T,jet}<100$ \\\\"<<std::endl;
        std::cout<<"\\hline"<<std::endl;
        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                if (gSystem->AccessPathName((output_folder + systematic_namef[available_systematics[syst_index]]).c_str()))
                        continue;
                
                fsyst[syst_index] = new TFile((output_folder + systematic_namef[available_systematics[syst_index]]).c_str());
                
                std::cout<<systematic_name[available_systematics[syst_index]]<<" & ";
                
                for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                        hsyst_eec_neqcharge[bin] = (TH1F*) fsyst[syst_index]->Get(Form("relerror_neqcheec%i",bin));

                        set_histo_with_systematics(hsyst_eec_neqcharge[bin], hcorr_eec_neqcharge[bin], hcorr_eec_neqcharge_syst[bin], syst_index);

                        if (bin != nbin_jet_pt-1) 
                                std::cout<<" & ";
                        else
                                std::cout<<" \\\\ ";
                }
                
                std::cout<<std::endl;

                delete fsyst[syst_index];
        }
        
        // Extract the statistical uncertainty from the plots with systematic uncertainties
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                substract_stat_error(hcorr_eec[bin], hcorr_eec_syst[bin]);
                substract_stat_error(hcorr_tau[bin], hcorr_tau_syst[bin]);
                substract_stat_error(hcorr_eec_eqcharge[bin], hcorr_eec_eqcharge_syst[bin]);
                substract_stat_error(hcorr_eec_neqcharge[bin], hcorr_eec_neqcharge_syst[bin]);
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1480);
        c->Draw();

        TLatex* tex = new TLatex();
        set_lhcb_watermark_properties(tex);

        // Adding content with errors
        TLatex* lhcbprint = new TLatex();
        TLatex latex;
        latex.SetTextAlign(22); // center alignment
        latex.SetTextSize(text_size_correction_plots);
        latex.SetTextColor(kBlack);

        gStyle->SetPaintTextFormat("4.2f");
        
        THStack* s_data = new THStack();
        TLegend* l_data = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());
        
        // Print the EECs as a function of tau
        double tau_binning[nbin_jet_pt][nbin_rl_nominal + 1];
        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                double avge_pt2_jet = (jet_pt_binning[i + 1] + jet_pt_binning[i])/2.;
                get_tau_binning_from_eec_binning(tau_binning[i], rl_nominal_binning, avge_pt2_jet);
        }

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_tau[bin],"E X0");
                s_data->Add(hcorr_tau_syst[bin],"E2");
                l_data->AddEntry(hcorr_tau[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"p");
        }
        
        TH1F* frame = gPad->DrawFrame(tau_binning[0][0], 0.0, tau_binning[2][nbin_rl_nominal], 0.09);
        s_data->Draw("NOSTACK SAME");
        frame->SetTitle(";R_{L}#LT p_{T,jet} #GT;#Sigma_{EEC}(R_{L})#times ln(#LT p_{T,jet} #GT)/#LT p_{T,jet} #GT");
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
        draw_lhcb_tag(lhcbprint);

        c->Print(Form("./plots/corrtau_unf-niter%i_incsyst_newparadigm.pdf",niter));

        // Print EECS
        s_data = new THStack();
        l_data = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_eec[bin],"E X0");
                s_data->Add(hcorr_eec_syst[bin],"E2");
                l_data->AddEntry(hcorr_eec[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lpf");
        }

        s_data->Draw("NOSTACK");
        s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        s_data->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
        s_data->SetMaximum(1.45);
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);

        draw_lhcb_tag(lhcbprint);

        c->Print(Form("./plots/correec_unf-niter%i_incsyst_newparadigm.pdf",niter));

        // Print charged EECs
        TLine* line = new TLine(rl_nominal_binning[0],0.5,rl_nominal_binning[nbin_rl_nominal],0.5);
        line->SetLineWidth(1);

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data = new THStack();
                l_data = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());
                
                s_data->Add(hcorr_eec_eqcharge[bin] ,"E X0");
                s_data->Add(hcorr_eec_eqcharge_syst[bin],"E2");
                s_data->Add(hcorr_eec_neqcharge[bin] ,"E X0");
                s_data->Add(hcorr_eec_neqcharge_syst[bin],"E2");
                
                l_data->SetHeader(Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l_data->AddEntry(hcorr_eec_eqcharge[bin], "eq. charge correlations","lfp");
                l_data->AddEntry(hcorr_eec_neqcharge[bin],"op. charge correlations","lfp");

                s_data->Draw("NOSTACK");
                s_data->SetTitle(";R_{L};");
                s_data->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
                s_data->SetMaximum(1.2);
                s_data->SetMinimum(0.1);
                l_data->Draw("SAME");
                line->Draw("SAME");
                gPad->SetLogx(1);
                gPad->SetLogy(0);
                
                draw_lhcb_tag(lhcbprint);

                c->Print(Form("./plots/corrchargedeec_jetptbin%i_unf-niter%i_incsyst_newparadigm.pdf",bin,niter));
        }
}
