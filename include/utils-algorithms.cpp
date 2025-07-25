#ifndef UTILS_ALGORITHMS_H
#define UTILS_ALGORITHMS_H

double get_hwhm(TH1F* h);

void determine_log10binning(int Nbins, double x_i, double x_f, double* binning);

void determine_eqsizebinning(int Nbins, double x_i, double x_f, double* binning);

void regularize_correction_factors(TH2F* h);

void regularize_correction_factors(TH2D* h);

void regularize_correction_factors(TH3F* h);

void regularize_correction_factors(TH3D* h);

void set_histo_with_systematics(TH1F* hdeviations, TH1F* hnominal, TH1F* hsystematic, std::string err_type);

void set_histo_sqrt_content(TH1F* h);

void set_histo_null_errors(TH1F* h);

void set_shift_histo(TH2F* href, TH2F* hshift, TRandom3* rndm);

void set_shift_histo(TH2D* href, TH2D* hshift, TRandom3* rndm);

void set_data_ntuple_branches(TNtuple* ntuple, float* R_L, float* jet_pt, float* weight_pt, float* efficiency, float* purity, float* efficiency_relerror, float* purity_relerror);

void set_unfolding_ntuple_branches(TNtuple* ntuple, float* R_L_reco, float* R_L_truth, float* jet_pt_reco, float* jet_pt_truth, float* weight_pt_reco, float* weight_pt_truth);

#endif