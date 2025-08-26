#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_determine_binning()
{
        double binning[nbin_chargedeec_nominal + 1];
        double binning_log[nbin_rl_nominal+1];
        double binning_corr_log[nbin_rl_altlogbin + 1];
        double binning_weight[nbin_weight + 1];
        double binning_ptprod[nbin_weight + 1];

        double binning_tau_log[nbin_tau_logbin + 1];

        determine_eqsizebinning(nbin_chargedeec_nominal, rl_logmin, rl_logmax, binning);
        determine_log10binning(nbin_weight, weight_absmin, weight_absmax, binning_weight);
        determine_log10binning(nbin_weight, ptprod_absmin, ptprod_absmax, binning_ptprod);
        determine_log10binning(nbin_rl_nominal, rl_logmin, rl_logmax, binning_log);
        determine_log10binning(nbin_tau_logbin, tau_min, tau_max, binning_tau_log);
        determine_log10binning(nbin_rl_altlogbin, rl_logmin, rl_logmax, binning_corr_log);
        
        // Pt prod
        std::cout<<"const double ptprod_binning[] = {ptprod_absmin";
        
        for (int i = 1 ; i < nbin_weight ; i++)
                std::cout<<", "<<binning_ptprod[i];

        std::cout<<", ptprod_absmax};"<<std::endl;

        // Weight
        std::cout<<"const double weight_binning[] = {weight_absmin";
        
        for (int i = 1 ; i < nbin_weight ; i++)
                std::cout<<", "<<binning_weight[i];

        std::cout<<", weight_absmax};"<<std::endl;

        // rl binning
        std::cout<<"const double charged_rl_chargedeec_binning[]              = {rl_logmin";
        
        for (int i = 1 ; i < nbin_chargedeec_nominal ; i++)
                std::cout<<", "<<binning[i];

        std::cout<<", rl_logmax};"<<std::endl;

        // rl log bin
        std::cout<<"const double rl_nominal_binning[]           = {rl_logmin";
        
        for (int i = 1 ; i < nbin_rl_nominal ; i++)
                std::cout<<", "<<binning_log[i];

        std::cout<<", rl_logmax};"<<std::endl;

        // rl alternative logbinning
        std::cout<<"const double rl_altlogbinning[]           = {rl_logmin";
        
        for (int i = 1 ; i < nbin_rl_altlogbin ; i++)
                std::cout<<", "<<binning_corr_log[i];

        std::cout<<", rl_logmax};"<<std::endl;

        // tau logbinning
        std::cout<<"const double tau_nominal_binning[]          = {tau_min";
        
        for (int i = 1 ; i < nbin_tau_logbin ; i++)
                std::cout<<", "<<binning_tau_log[i];
        
        std::cout<<", tau_max};"<<std::endl;

        return 0;
}