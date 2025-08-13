#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

double ptproduct_max = 700;
// const double ptproduct_binning[] = {0.0625, 0.611951, 1.05991, 1.60586, 2.2988, 3.25771, 4.65059, 6.85539, 10.81, 20.3782, 700};

void macro_determine_binning()
{
        // Open data file 
        TFile* f = new TFile((output_folder + namef_ntuple_eec_corr).c_str());
        
        // Get TNtuple from file
        TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_data.c_str());
        
        // Declare the histos to use
        TH1F* h_weight = new TH1F("h_weight","",100000,track_pt_min*track_pt_min,ptproduct_max);
        
        // Put the data into the histos and get the cumulative distributions
        ntuple->Project("h_weight","h1_pt*h2_pt",pair_cut);

        TH1F* h_weight_cumul = (TH1F*) h_weight->GetCumulative();
        // h_weight_cumul->Draw();

        // Determine binning in weights
        double entries_bin_weight = h_weight->Integral()/nbin_pt_product;
        std::cout<<"The number of entries for weight is "<<entries_bin_weight<<std::endl;
        std::cout<<"const double pt_product_binning[] = {"<<track_pt_min*track_pt_min<<"";
        int counter_weight = 1;
        for (int bin = 1 ; bin <= h_weight->GetNbinsX() ; bin++) {
                double q = h_weight_cumul->GetBinContent(bin);

                // Exit when determined last bin
                if (counter_weight == nbin_pt_product) {
                        std::cout<<", "<<ptproduct_max<<"};"<<std::endl;
                        break;
                }

                // Condition to determine limit
                if (q > entries_bin_weight * counter_weight) {
                        counter_weight++;

                        // Print limits
                        std::cout<<", "<<h_weight_cumul->GetBinCenter(bin-1);
                }
        }

        // // Declare the histos to use
        // TH1F* h_weight = new TH1F("h_weight","",nbin_pt_product,ptproduct_binning);
        // ntuple->Project("h_weight","h1_pt*h2_pt",pair_cut);
        // h_weight->Draw();
        
        
        // // Close file
        // f->Close();

        // const int Nbin = nbin_rl_nominal;
        // double binning[Nbin + 1];
        
        // const int nbin_log = nbin_rl_nominal;
        // double binning_log[nbin_log+1];
        // double binning_corr_log[nbin_rl_altlogbin + 1];

        // double binning_tau_log[nbin_tau_logbin + 1];

        // determine_eqsizebinning(nbin_log, rl_logmin, rl_logmax, binning);
        // determine_log10binning(nbin_log, rl_logmin, rl_logmax, binning_log);
        // determine_log10binning(nbin_tau_logbin, tau_min, tau_max, binning_tau_log);
        // determine_log10binning(nbin_rl_altlogbin, rl_logmin, rl_logmax, binning_corr_log);
        
        // // rl binning
        // std::cout<<"const double charged_rl_chargedeec_binning[]              = {rl_logmin";
        
        // for (int i = 1 ; i < Nbin ; i++)
        //         std::cout<<", "<<binning[i];

        // std::cout<<", rl_logmax};"<<std::endl;

        // // rl log bin
        // std::cout<<"const double rl_nominal_binning[]           = {rl_logmin";
        
        // for (int i = 1 ; i < nbin_log ; i++)
        //         std::cout<<", "<<binning_log[i];

        // std::cout<<", rl_logmax};"<<std::endl;

        // // rl alternative logbinning
        // std::cout<<"const double rl_altlogbinning[]           = {rl_logmin";
        
        // for (int i = 1 ; i < nbin_rl_altlogbin ; i++)
        //         std::cout<<", "<<binning_corr_log[i];

        // std::cout<<", rl_logmax};"<<std::endl;

        // // tau logbinning
        // std::cout<<"const double tau_nominal_binning[]          = {tau_min";
        
        // for (int i = 1 ; i < nbin_tau_logbin ; i++)
        //         std::cout<<", "<<binning_tau_log[i];
        
        // std::cout<<", tau_max};"<<std::endl;

        return 0;
}