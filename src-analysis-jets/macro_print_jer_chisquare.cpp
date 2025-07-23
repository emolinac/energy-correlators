#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_jer_chisquare(const int nbin = 50, double ptratio_min = 0.4 , double ptratio_max = 1.6, bool do_print = true)
{
        // Open the necessary files
        TFile* f = new TFile((output_folder + namef_ntuple_jes_jer).c_str());
        
        // Get the corresponding Ntuples
        TNtuple* ntuple_jes_data = (TNtuple*) f->Get((name_ntuple_jes_data).c_str());
        TNtuple* ntuple_jes_reco = (TNtuple*) f->Get((name_ntuple_jes_reco).c_str());

        float jet_pt_reco, z_pt_reco, jet_jec_cor_reco;
        ntuple_jes_reco->SetBranchAddress("jet_pt"     ,&jet_pt_reco);
        ntuple_jes_reco->SetBranchAddress("z_pt"       ,&z_pt_reco);
        ntuple_jes_reco->SetBranchAddress("jet_jec_cor",&jet_jec_cor_reco);
        
        // Define the necessary histograms to calculate purity
        TH1F* hdata_nojec[nbin_jet_pt]; 
        TH1F* hreco_nojec[nbin_jet_pt];  
        TH1F* hreco_newjec[nbin_jet_pt];  
        TH1F* halphastar_balance[nbin_jet_pt];  
        TH1F* halphastar_chisquare[nbin_jet_pt];

        TGraph* chisquare_graph[nbin_jet_pt];

        THStack* hs[nbin_jet_pt];
        TLegend* l[nbin_jet_pt];
        TCanvas* c = new TCanvas("c","",1500,500);
        c->Draw();
        c->Divide(3,1);
        
        const double alpha_star_init = 0.;
        const double alpha_star_end  = 0.5;
        const int    alpha_star_bins = 30;
        double alpha_star_step = (alpha_star_end - alpha_star_init)/alpha_star_bins;

        TRandom3* rndm = new TRandom3();
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                c->cd(bin + 1);
                
                hs[bin] = new THStack();
                l[bin]  = new TLegend(gPad->GetLeftMargin()+0.01,0.8,gPad->GetLeftMargin()+0.26,0.9,Form(" %.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]));

                hdata_nojec[bin]          = new TH1F(Form("hdata_nojec[%i]",bin)         ,"", nbin           , ptratio_min    , ptratio_max); 
                hreco_nojec[bin]          = new TH1F(Form("hreco_nojec[%i]",bin)         ,"", nbin           , ptratio_min    , ptratio_max); 
                hreco_newjec[bin]         = new TH1F(Form("hreco_newjec[%i]",bin)        ,"", nbin           , ptratio_min    , ptratio_max); 
                halphastar_balance[bin]   = new TH1F(Form("halphastar_balance[%i]",bin)  ,"", alpha_star_bins , alpha_star_init , alpha_star_end); 
                halphastar_chisquare[bin] = new TH1F(Form("halphastar_chisquare[%i]",bin),"", alpha_star_bins , alpha_star_init , alpha_star_end); 
        
                set_histogram_style(halphastar_balance[bin]  , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size);
                set_histogram_style(halphastar_chisquare[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size);
                
                // Undo the JEC
                ntuple_jes_data->Project(Form("hdata_nojec[%i]",bin),"(jet_pt/z_pt)/jet_jec_cor",pair_jetpt_cut[bin]);
                ntuple_jes_reco->Project(Form("hreco_nojec[%i]",bin),"(jet_pt/z_pt)/jet_jec_cor",pair_jetpt_cut[bin]);
                // ntuple_jes_data->Project(Form("hdata_nojec[%i]",bin),"(jet_pt/z_pt)/jet_jec_cor",Form("(jet_pt/jet_jec_cor)>%f&&(jet_pt/jet_jec_cor)<%f",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                // ntuple_jes_reco->Project(Form("hreco_nojec[%i]",bin),"(jet_pt/z_pt)/jet_jec_cor",Form("(jet_pt/jet_jec_cor)>%f&&(jet_pt/jet_jec_cor)<%f",jet_pt_binning[bin],jet_pt_binning[bin + 1]));

                for (int alpha_star_bin = 0 ; alpha_star_bin < alpha_star_bins ; alpha_star_bin++) {
                        double alpha_star = alpha_star_init + alpha_star_bin*alpha_star_step;

                        for (int entry = 0 ; entry < ntuple_jes_reco->GetEntries() ; entry++) {
                                ntuple_jes_reco->GetEntry(entry);

                                if (jet_pt_reco < jet_pt_binning[bin] || jet_pt_reco > jet_pt_binning[bin + 1]) 
                                        continue;

                                hreco_newjec[bin]->Fill(rndm->Gaus(1,alpha_star)*jet_pt_reco/z_pt_reco/jet_jec_cor_reco);
                        }

                        double delta_mean = abs(hreco_newjec[bin]->GetStdDev() - hdata_nojec[bin]->GetStdDev());
                        double chisquare  = hreco_newjec[bin]->Chi2Test(hdata_nojec[bin],"CHI2");

                        halphastar_balance[bin]->SetBinContent(alpha_star_bin + 1, delta_mean);
                        halphastar_balance[bin]->SetBinError(alpha_star_bin + 1,sqrt(hreco_newjec[bin]->GetStdDevError()*hreco_newjec[bin]->GetStdDevError() + hdata_nojec[bin]->GetStdDevError()*hdata_nojec[bin]->GetStdDevError()));
                        halphastar_chisquare[bin]->SetBinContent(alpha_star_bin + 1, chisquare);
                }

                double alpha_star_min = halphastar_balance[bin]->GetBinCenter(halphastar_balance[bin]->GetMinimumBin());

                halphastar_balance[bin]->Draw();
                halphastar_balance[bin]->SetMinimum(-0.006);
                halphastar_balance[bin]->SetMaximum(0.036);
                halphastar_balance[bin]->SetTitle(";#alpha;");

                l[bin]->Clear();
                l[bin]->SetHeader(Form(" %.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l[bin]->AddEntry(halphastar_balance[bin],Form("#Delta #sigma(p^{jet}_{t}/p^{Z}_{t}), #alpha^{*}=%.4f",alpha_star_min),"p");
                l[bin]->Draw("SAME");
        }

        if (do_print) 
                c->Print(Form("./plots/jer_alpha_balance.pdf"));

        std::cout<<"COPY INTO sys-jes-jer.h ---> const double syst_jer_array[] = { ";
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                c->cd(bin + 1);

                double alpha_star_min = halphastar_chisquare[bin]->GetBinCenter(halphastar_chisquare[bin]->GetMinimumBin());
                
                std::cout<<alpha_star_min; if (bin<nbin_jet_pt-1) std::cout<<", ";

                halphastar_chisquare[bin]->Draw("E");
                halphastar_chisquare[bin]->SetMinimum(20);
                halphastar_chisquare[bin]->SetMaximum(120);
                halphastar_chisquare[bin]->SetTitle(";#alpha;");
                
                l[bin]->Clear();
                l[bin]->SetHeader(Form(" %.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                l[bin]->AddEntry(halphastar_chisquare[bin],Form("#chi^{2}(p^{jet}_{t}/p^{Z}_{t}), #alpha^{*}=%.4f",alpha_star_min),"p");
                l[bin]->Draw("SAME");
        }

        std::cout<<"};"<<std::endl;

        if (do_print) 
                c->Print(Form("./plots/jer_alpha_chisquare.pdf"));
}