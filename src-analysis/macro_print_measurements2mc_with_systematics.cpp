#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_measurements2mc_with_systematics(int niter = 4, int niter_jet = 4)
{
        TFile* fnominal = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niter%i_niterjet%i--get-nominal.root",niter,niter_jet)).c_str());
        TFile* fmc      = new TFile((output_folder + "ntuple_mc_eec.root").c_str());

        TH1F* hcorr_eec[nbin_jet_pt]; 
        TH1F* hcorr_eec_syst[nbin_jet_pt]; 
        TH1F* hcorr_tau[nbin_jet_pt]; 
        TH1F* hcorr_tau_syst[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge_syst[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge_syst[nbin_jet_pt]; 

        TNtuple* ntuple_mc     = (TNtuple*) fmc->Get(name_ntuple_mc.c_str());
        TNtuple* ntuple_mc_jet = (TNtuple*) fmc->Get(name_ntuple_mc_jet.c_str());
        
        TH1F* hmc_eec[nbin_jet_pt]; 
        TH1F* hmc_eec_eqcharge[nbin_jet_pt]; 
        TH1F* hmc_eec_neqcharge[nbin_jet_pt]; 
        
        TH2D* hmc_eec_2d      = new TH2D("hmc_eec_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hmc_eqcheec_2d  = new TH2D("hmc_eqcheec_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH2D* hmc_neqcheec_2d = new TH2D("hmc_neqcheec_2d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
        TH3D* hmc_eec_3d      = new TH3D("hmc_eec_3d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* hmc_eqcheec_3d  = new TH3D("hmc_eqcheec_3d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3D* hmc_neqcheec_3d = new TH3D("hmc_neqcheec_3d","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning, nbin_weight, weight_binning);

        TH1F* hdatamcratio_eec[nbin_jet_pt]; 
        TH1F* hdatamcratio_eec_eqcharge[nbin_jet_pt]; 
        TH1F* hdatamcratio_eec_neqcharge[nbin_jet_pt]; 

        // Get the MC histos
        ntuple_mc->Project("hmc_eec_3d"     , "weight_pt:jet_pt:R_L");
        ntuple_mc->Project("hmc_eqcheec_3d" , "weight_pt:jet_pt:R_L","eq_charge > 0");
        ntuple_mc->Project("hmc_neqcheec_3d", "weight_pt:jet_pt:R_L","eq_charge < 0");

        apply_unfolded_weights(hmc_eec_3d, hmc_eec_2d);
        apply_unfolded_weights(hmc_eqcheec_3d, hmc_eqcheec_2d);
        apply_unfolded_weights(hmc_neqcheec_3d, hmc_neqcheec_2d);

        std::cout<<"Apply weights in MC"<<std::endl;

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                int nominal_jet_pt_bin = bin + 3;

                hmc_eec[bin]           = new TH1F(Form("hmc_eec%i",bin)          ,"", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hmc_eec_eqcharge[bin]  = new TH1F(Form("hmc_eec_eqcharge%i",bin) ,"", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hmc_eec_neqcharge[bin] = new TH1F(Form("hmc_eec_neqcharge%i",bin),"", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);

                project_nominal_phase_space(hmc_eec_2d     , hmc_eec[bin]          , nominal_jet_pt_bin);
                project_nominal_phase_space(hmc_eqcheec_2d , hmc_eec_eqcharge[bin] , nominal_jet_pt_bin);
                project_nominal_phase_space(hmc_neqcheec_2d, hmc_eec_neqcharge[bin], nominal_jet_pt_bin);

                std::cout<<"Project nominal space"<<std::endl;

                double njets_mc = ntuple_mc_jet->GetEntries(pair_jet_pt_cut[bin]);
                hmc_eec[bin]->Scale(1./njets_mc,"width");
                hmc_eec_eqcharge[bin]->Scale(1./njets_mc,"width");
                hmc_eec_neqcharge[bin]->Scale(1./njets_mc,"width");
                
                hmc_eec_eqcharge[bin]->Divide(hmc_eec[bin]);
                hmc_eec_neqcharge[bin]->Divide(hmc_eec[bin]);

                std::cout<<"NOrmalizing"<<std::endl;

                set_histogram_style(hmc_eec[bin]           , corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hmc_eec_eqcharge[bin]  , corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin], std_marker_size+2);
                set_histogram_style(hmc_eec_neqcharge[bin] , corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin], std_marker_size+2);
        }
        
        std::cout<<"MC done"<<std::endl;

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

        std::cout<<"Data done"<<std::endl;

        // Include the systematics in the whole deal
        const int nsyst = sizeof(available_systematics)/sizeof(available_systematics[0]);
        TFile* fsyst[nsyst];
        TH1F* hsyst_eec[nbin_jet_pt];
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
        
        std::cout<<"Drawing starts"<<std::endl;

        // NO subtraction of stat errors in this code because we want to take the ratio with the total error of the data
        TCanvas* c = new TCanvas("c", "", 1920, 1680);
        TPad* pu = new TPad("pu","",0,0.4,1,1);
        TPad* pd = new TPad("pd","",0,0.,1,0.4);

        c->Draw();
        pu->Draw();
        pu->cd();
        gPad->SetBottomMargin(0.);

        // Adding content with errors
        TLatex* lhcbprint = new TLatex();
        
        THStack* su = new THStack();
        THStack* sd = new THStack();
        TLegend* l_data = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());

        // Print EECS
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                su->Add(hcorr_eec_syst[bin],"E2 P");
                su->Add(hmc_eec[bin],"EP");
                l_data->AddEntry(hcorr_eec_syst[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lpf");
        }

        su->Draw("NOSTACK");
        su->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
        su->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
        su->SetMaximum(1.45);
        su->SetMinimum(0.01);
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);

        draw_lhcb_tag(lhcbprint);

        c->cd();
        pd->Draw();
        pd->cd();
        gPad->SetTopMargin(0);
        gPad->SetBottomMargin(0.3);
        
        TLine* line = new TLine(unfolding_rl_nominal_binning[1], 1, unfolding_rl_nominal_binning[nbin_rl_nominal_unfolding-1], 1);
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                hdatamcratio_eec[bin] = new TH1F(Form("hdatamcratio_eec%i",bin),"",nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hdatamcratio_eec[bin]->Divide(hcorr_eec_syst[bin],hmc_eec[bin],1,1);
                
                set_histogram_style(hdatamcratio_eec[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                hdatamcratio_eec[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);

                sd->Add(hdatamcratio_eec[bin],"E2 P");
        }

        sd->Draw("NOSTACK");
        sd->SetTitle(";R_{L};Data/Theory");
        sd->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
        sd->SetMaximum(1.5);
        sd->SetMinimum(0.5);
        gPad->SetLogx(1);
        gPad->SetLogy(0);

        line->Draw("SAME");
        
        c->Print(Form("./plots/correec_unf-niter%i_2dunf_incsyst_data2mc.pdf",niter));

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                c->cd();
                pu->Draw();
                pu->cd();
                gPad->SetBottomMargin(0.);

                su = new THStack();
                sd = new THStack();
                l_data = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());


                su->Add(hcorr_eec_eqcharge_syst[bin],"E2 P");
                su->Add(hmc_eec_eqcharge[bin],"HISTO C");
                su->Add(hcorr_eec_neqcharge_syst[bin],"E2 P");
                su->Add(hmc_eec_neqcharge[bin],"HISTO C");
                
                l_data->SetHeader(Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l_data->AddEntry(hcorr_eec_eqcharge_syst[bin], "eq. charge correlations","lfp");
                l_data->AddEntry(hcorr_eec_neqcharge_syst[bin],"op. charge correlations","lfp");

                su->Draw("NOSTACK");
                su->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
                su->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
                su->SetMaximum(1.45);
                su->SetMinimum(0.1);
                l_data->Draw("SAME");
                gPad->SetLogx(1);
                gPad->SetLogy(0);

                draw_lhcb_tag(lhcbprint);

                c->cd();
                pd->Draw();
                pd->cd();
                gPad->SetTopMargin(0);
                gPad->SetBottomMargin(0.3);
                
                hdatamcratio_eec_eqcharge[bin] = new TH1F(Form("hdatamcratio_eec_eqcharge%i",bin),"",nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hdatamcratio_eec_eqcharge[bin]->Divide(hcorr_eec_eqcharge_syst[bin],hmc_eec_eqcharge[bin],1,1);
                
                hdatamcratio_eec_neqcharge[bin] = new TH1F(Form("hdatamcratio_eec_neqcharge%i",bin),"",nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hdatamcratio_eec_neqcharge[bin]->Divide(hcorr_eec_neqcharge_syst[bin],hmc_eec_neqcharge[bin],1,1);
                
                set_histogram_style(hdatamcratio_eec_eqcharge[bin], corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin], std_marker_size+2);
                hdatamcratio_eec_eqcharge[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);

                set_histogram_style(hdatamcratio_eec_neqcharge[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+2);
                hdatamcratio_eec_neqcharge[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);

                sd->Add(hdatamcratio_eec_eqcharge[bin],"E2 P");
                sd->Add(hdatamcratio_eec_neqcharge[bin],"E2 P");

                sd->Draw("NOSTACK");
                sd->SetTitle(";R_{L};Data/Theory");
                sd->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
                sd->SetMaximum(1.5);
                sd->SetMinimum(0.5);
                gPad->SetLogx(1);
                gPad->SetLogy(0);

                line->Draw("SAME");

                c->Print(Form("./plots/corrchargedeec_unf-niter%i_2dunf_incsyst_data2mc_jetptbin%i.pdf",niter,bin));
        }
}
