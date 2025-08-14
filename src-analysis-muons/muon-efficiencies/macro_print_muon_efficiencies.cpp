#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_muon_efficiencies()
{
        TFile* fefficiency_muon_2016_id  = new TFile((muons_folder + "IDEff_Data_2016.root").c_str());
        TFile* fefficiency_muon_2016_trk = new TFile((muons_folder + "TRKEff_Data_2016.root").c_str());
        TFile* fefficiency_muon_2016_trg = new TFile((muons_folder + "TRGEff_Data_2016.root").c_str());
        
        TH2D* h2_muon_2016_ideff_data  = (TH2D*) fefficiency_muon_2016_id->Get("Hist_ALL_2016_ETA_PT_Eff");
        TH2D* h2_muon_2016_trkeff_data = (TH2D*) fefficiency_muon_2016_trk->Get("Hist_ALL_2016_ETA_PT_Eff");
        TH2D* h2_muon_2016_trgeff_data = (TH2D*) fefficiency_muon_2016_trg->Get("Hist_ALL_2016_ETA_PT_Eff");

        h2_muon_2016_ideff_data->SetTitle("Muon ID Efficiency;#eta;p_{T}(#mu)");
        h2_muon_2016_trkeff_data->SetTitle("Muon Tracking Efficiency;#eta;p_{T}(#mu)");
        h2_muon_2016_trgeff_data->SetTitle("Muon Trigger Efficiency;#eta;p_{T}(#mu)");

        TCanvas* c = new TCanvas("c", "", 1920, 1080);
        c->Draw();
        
        // Adding content with errors
        TLatex latex;
        latex.SetTextAlign(22); // center alignment
        latex.SetTextSize(text_size_correction_plots);
        latex.SetTextColor(kBlack);
        gStyle->SetPaintTextFormat("4.2f");
        
        h2_muon_2016_ideff_data->Draw("col");
        for (int i = 1; i <= h2_muon_2016_ideff_data->GetNbinsX(); ++i) {
                for (int j = 1; j <= h2_muon_2016_ideff_data->GetNbinsY(); ++j) {
                        double x = h2_muon_2016_ideff_data->GetXaxis()->GetBinCenter(i);
                        double y = h2_muon_2016_ideff_data->GetYaxis()->GetBinCenter(j);
                        double content = h2_muon_2016_ideff_data->GetBinContent(i, j);
                        double error = h2_muon_2016_ideff_data->GetBinError(i, j);
                        // Draw content and error in the format "content ± error"
                        latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                }
        }

        c->Print("./plots/muon_ideff.pdf");

        h2_muon_2016_trkeff_data->Draw("col");
        for (int i = 1; i <= h2_muon_2016_trkeff_data->GetNbinsX(); ++i) {
                for (int j = 1; j <= h2_muon_2016_trkeff_data->GetNbinsY(); ++j) {
                        double x = h2_muon_2016_trkeff_data->GetXaxis()->GetBinCenter(i);
                        double y = h2_muon_2016_trkeff_data->GetYaxis()->GetBinCenter(j);
                        double content = h2_muon_2016_trkeff_data->GetBinContent(i, j);
                        double error = h2_muon_2016_trkeff_data->GetBinError(i, j);
                        // Draw content and error in the format "content ± error"
                        latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                }
        }

        c->Print("./plots/muon_trkeff.pdf");

        h2_muon_2016_trgeff_data->Draw("col");
        for (int i = 1; i <= h2_muon_2016_trgeff_data->GetNbinsX(); ++i) {
                for (int j = 1; j <= h2_muon_2016_trgeff_data->GetNbinsY(); ++j) {
                        double x = h2_muon_2016_trgeff_data->GetXaxis()->GetBinCenter(i);
                        double y = h2_muon_2016_trgeff_data->GetYaxis()->GetBinCenter(j);
                        double content = h2_muon_2016_trgeff_data->GetBinContent(i, j);
                        double error = h2_muon_2016_trgeff_data->GetBinError(i, j);
                        // Draw content and error in the format "content ± error"
                        latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                }
        }

        c->Print("./plots/muon_trgeff.pdf");
}