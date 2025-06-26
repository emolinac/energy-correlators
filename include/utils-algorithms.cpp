#ifndef UTILS_ALGORITHMS_H
#define UTILS_ALGORITHMS_H

void determine_log10binning(int Nbins, double x_i, double x_f, double* binning);

void determine_eqsizebinning(int Nbins, double x_i, double x_f, double* binning);

double get_hwhm(TH1F* h);

void get_histo_with_systematics(TH1F* hdeviations, TH1F* hnominal, TH1F* hsystematic);

#endif