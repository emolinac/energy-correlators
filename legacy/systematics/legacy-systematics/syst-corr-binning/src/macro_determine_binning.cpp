#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_determine_binning()
{
    // Open data file 
    TFile* f = new TFile((output_folder + namef_ntuple_e2c_corr).c_str());
    
    // Get TNtuple from file
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_data.c_str());
    
    // Declare the histos to use
    TH1F* h_weight = new TH1F("h_weight","",100000000,weight_min,weight_max);

    // Put the data into the histos and get the cumulative distributions
    ntuple->Project("h_weight","weight_pt",pair_cut);

    TH1F* h_weight_cumul = (TH1F*) h_weight->GetCumulative();

    // Determine binning in weights
    double entries_bin_weight = h_weight->Integral()/nbin_weight;
    std::cout<<"The number of entries for weight is "<<entries_bin_weight<<std::endl;
    std::cout<<"const double weight_binning[] = {"<<weight_min<<"";
    int counter_weight = 1;
    for (int bin = 1 ; bin <= h_weight->GetNbinsX() ; bin++)
    {
        double q = h_weight_cumul->GetBinContent(bin);

        // Exit when determined last bin
        if (counter_weight==nbin_weight) 
        {
            std::cout<<", "<<weight_max<<"};"<<std::endl;
            break;
        }
        // Condition to determine limit
        if (q>entries_bin_weight*counter_weight) 
        {
            counter_weight++;

            // Print limits
            std::cout<<", "<<h_weight_cumul->GetBinCenter(bin-1);
        }
    }
    
    // Close file
    f->Close();

    const int Nbin = nbin_rl;
    double binning[Nbin + 1];
    
    const int nbin_log = nbin_rl_nominal;
    double binning_log[nbin_log+1];
    double binning_corr_log[nbin_log*2+1];

    double binning_tau_log[nbin_tau_logbin + 1];

    determine_eqsizebinning(Nbin, rl_min, rl_max, binning);
    determine_log10binning(nbin_log, rl_logmin, rl_logmax, binning_log);
    determine_log10binning(nbin_tau_logbin, tau_min, tau_max, binning_tau_log);

    determine_log10binning(nbin_log*2, rl_logmin, rl_logmax, binning_corr_log);
    std::cout<<"const double rl_binning[]              = {rl_min";
    for (int i = 1 ; i < Nbin ; i++)
    {
        std::cout<<", "<<binning[i];
    }
    std::cout<<", rl_max};"<<std::endl;

    std::cout<<"const double rl_nominal_binning[]           = {rl_logmin";
    for (int i = 1 ; i < nbin_log ; i++)
    {
        std::cout<<", "<<binning_log[i];
    }
    std::cout<<", rl_logmax};"<<std::endl;

    std::cout<<"const double tau_nominal_binning[]          = {tau_min";
    for (int i = 1 ; i < nbin_tau_logbin ; i++)
    {
        std::cout<<", "<<binning_tau_log[i];
    }
    std::cout<<", tau_max};"<<std::endl;

    return 0;
}