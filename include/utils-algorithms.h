#ifndef UTILS_ALGORITHMS_H
#define UTILS_ALGORITHMS_H

#include "TRandom3.h"

double get_hwhm(TH1F* h)
{
        // As written, this function works for single peak distributions
        double half_max = h->GetMaximum()/2.;
    
        int halfwidth_bin;
        for (int bin = 1 ; bin <= h->GetNbinsX() ; bin++)
                if (h->GetBinContent(bin) >= half_max) {
                        halfwidth_bin = bin; 
                        break;
                }
        
        
        // Return the Half width at half maximum value
        return abs(h->GetBinCenter(h->GetMaximumBin()) - h->GetBinCenter(halfwidth_bin));
}

void determine_log10binning(int Nbins, double x_i, double x_f, double* binning)
{
        double log_x_i = log10(x_i);
        double log_x_f = log10(x_f);
        double delta = (log_x_f-log_x_i)/Nbins;

        for (int i = 0 ; i <=Nbins ; i++) 
                binning[i] = pow(10,log_x_i+delta*i);

        return;
}

void determine_eqsizebinning(int Nbins, double x_i, double x_f, double* binning)
{
        double delta = (x_f - x_i)/Nbins;

        for (int i = 0;i <= Nbins;i++)
                binning[i] = x_i + delta*i;
        
        return;
}


void regularize_correction_factors(TH2F* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h->GetNbinsY() ; j++) {
                        double bin_content = h->GetBinContent(i, j);

                        if (bin_content > 1)
                                h->SetBinContent(i, j, 1.);
                }
        }
}

void regularize_correction_factors(TH2D* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h->GetNbinsY() ; j++) {
                        double bin_content = h->GetBinContent(i, j);

                        if (bin_content > 1)
                                h->SetBinContent(i, j, 1.);
                }
        }
}

void regularize_correction_factors(TH3F* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h->GetNbinsY() ; j++) {
                        for (int k = 1 ; k <= h->GetNbinsZ() ; k++){
                                double bin_content = h->GetBinContent(i, j, k);

                                if (bin_content > 1)
                                        h->SetBinContent(i, j, k, 1.);
                        }
                }
        }
}

void regularize_correction_factors(TH3D* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h->GetNbinsY() ; j++) {
                        for (int k = 1 ; k <= h->GetNbinsZ() ; k++){
                                double bin_content = h->GetBinContent(i, j, k);

                                if (bin_content > 1)
                                        h->SetBinContent(i, j, k, 1.);
                        }
                }
        }
}

void set_histo_with_systematics(TH1F* hdeviations, TH1F* hnominal, TH1F* hsystematic, std::string err_type = "normal", bool print_table = true)
{
        double total_err = 0, total = 0;

        for (int hbin = 1 ; hbin <= hdeviations->GetNbinsX() ; hbin++)
        {
                double dev     = hdeviations->GetBinContent(hbin);
                double dev_err = hdeviations->GetBinError(hbin);

                total += hnominal->GetBinContent(hbin);

                // Demand more than one sigma to be considered
                if ((dev+dev_err > 1 && dev < 1) || (dev-dev_err < 1 && dev > 1)) 
                        continue;

                double syst_error;
                double syst_error_percentage;
                
                if (err_type=="uniform") {
                        syst_error            = abs(1. - dev)*hnominal->GetBinContent(hbin)/sqrt(12.);
                        syst_error_percentage = abs(1. - dev)/sqrt(12.);
                } else {
                        syst_error            = abs(1. - dev)*hnominal->GetBinContent(hbin);
                        syst_error_percentage = abs(1. - dev);
                }
                
                hsystematic->SetBinError(hbin, sqrt(syst_error*syst_error + hnominal->GetBinError(hbin)*hnominal->GetBinError(hbin)));

                total_err += hnominal->GetBinContent(hbin)*syst_error_percentage;
        }

        if (!print_table)
                return;

        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        std::cout.precision(2);

        if (total > 0)
                std::cout<<total_err*100./total;
        else if (total == 0)
                std::cout<<0;
}

void set_histo_sqrt_content(TH1F* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                h->SetBinContent(i, sqrt(h->GetBinContent(i)));
                h->SetBinError(i, sqrt(h->GetBinError(i)));
        }
}

void set_histo_null_errors(TH1F* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++)
                h->SetBinError(i, 0);
}

void set_shift_histo(TH2F* href, TH2F* hshift, TRandom3* rndm)
{
        for (int xbin = 1 ; xbin <= href->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= href->GetNbinsY() ; ybin++) {
                        double shift_window = href->GetBinError(href->GetBin(xbin, ybin));
                        double refdata      = href->GetBinContent(href->GetBin(xbin, ybin));
                        double shift        = rndm->Gaus(refdata, shift_window)/refdata;
                        hshift->SetBinContent(xbin, ybin, shift);
                }
        }
}

void set_shift_histo(TH2D* href, TH2D* hshift, TRandom3* rndm)
{
        for (int xbin = 1 ; xbin <= href->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= href->GetNbinsY() ; ybin++) {
                        double shift_window = href->GetBinError(href->GetBin(xbin, ybin));
                        double refdata      = href->GetBinContent(href->GetBin(xbin, ybin));
                        double shift        = rndm->Gaus(refdata, shift_window)/refdata;

                        hshift->SetBinContent(xbin, ybin, shift);
                        hshift->SetBinError(xbin, ybin, shift_window/refdata);
                }
        }
}

void set_shift_histo(TH3D* href, TH3D* hshift, TRandom3* rndm)
{
        for (int xbin = 1 ; xbin <= href->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= href->GetNbinsY() ; ybin++) {
                        for (int zbin = 1 ; zbin <= href->GetNbinsZ() ; zbin++) {
                                double shift_window = href->GetBinError(href->GetBin(xbin, ybin, zbin));
                                double refdata      = href->GetBinContent(href->GetBin(xbin, ybin, zbin));
                                double shift        = rndm->Gaus(refdata, shift_window)/refdata;

                                if (std::isnan(shift)){
                                        hshift->SetBinContent(xbin, ybin, zbin, 0);
                                        hshift->SetBinError(xbin, ybin, zbin, 0);

                                        continue;
                                }
                                
                                hshift->SetBinContent(xbin, ybin, zbin, shift);
                                hshift->SetBinError(xbin, ybin, zbin, shift_window/refdata);
                        }
                }
        }
}

void set_data_ntuple_branches(TNtuple* ntuple, float* event_weight, float* R_L, float* jet_pt, float* weight_pt, float* efficiency, float* purity, float* efficiency_relerror, float* purity_relerror)
{
        ntuple->SetBranchAddress("event_weight", event_weight);
        ntuple->SetBranchAddress("R_L", R_L);
        ntuple->SetBranchAddress("jet_pt", jet_pt);
        ntuple->SetBranchAddress("weight_pt", weight_pt);
        ntuple->SetBranchAddress("efficiency", efficiency);
        ntuple->SetBranchAddress("efficiency_relerror", efficiency_relerror);
        ntuple->SetBranchAddress("purity", purity);
        ntuple->SetBranchAddress("purity_relerror", purity_relerror);
}

void set_data_ntuple_branches(TNtuple* ntuple, float* event_weight, float* R_L, float* jet_pt, float* weight_pt, float* efficiency, float* purity, float* efficiency_relerror, float* purity_relerror, float* eq_charge)
{
        ntuple->SetBranchAddress("event_weight", event_weight);
        ntuple->SetBranchAddress("R_L", R_L);
        ntuple->SetBranchAddress("jet_pt", jet_pt);
        ntuple->SetBranchAddress("weight_pt", weight_pt);
        ntuple->SetBranchAddress("efficiency", efficiency);
        ntuple->SetBranchAddress("efficiency_relerror", efficiency_relerror);
        ntuple->SetBranchAddress("purity", purity);
        ntuple->SetBranchAddress("purity_relerror", purity_relerror);
        ntuple->SetBranchAddress("eq_charge", eq_charge);
}

void set_unfolding_ntuple_branches(TNtuple* ntuple, float* R_L_reco, float* R_L_truth, float* jet_pt_reco, float* jet_pt_truth, float* weight_pt_reco, float* weight_pt_truth)
{
        ntuple->SetBranchAddress("jet_pt", jet_pt_reco);
        ntuple->SetBranchAddress("jet_pt_truth", jet_pt_truth);
        ntuple->SetBranchAddress("R_L", R_L_reco);
        ntuple->SetBranchAddress("R_L_truth", R_L_truth);
        ntuple->SetBranchAddress("weight_pt", weight_pt_reco);
        ntuple->SetBranchAddress("weight_pt_truth", weight_pt_truth);
}

void set_unfolding_jet_ntuple_branches(TNtuple* ntuple, float* jet_pt_reco, float* jet_pt_truth)
{
        ntuple->SetBranchAddress("jet_pt", jet_pt_reco);
        ntuple->SetBranchAddress("jet_pt_truth", jet_pt_truth);
}

#endif