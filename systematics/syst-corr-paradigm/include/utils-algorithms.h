#ifndef UTILS_ALGORITHMS_H
#define UTILS_ALGORITHMS_H

void determine_log10binning(int Nbins, double x_i, double x_f, double* binning)
{
    double log_x_i = log10(x_i);
    double log_x_f = log10(x_f);
    double delta = (log_x_f-log_x_i)/Nbins;

    for(int i = 0 ; i <=Nbins ; i++)
    {
        binning[i] = pow(10,log_x_i+delta*i);
    }

    return;
}

void determine_eqsizebinning(int Nbins, double x_i, double x_f, double* binning)
{
    double delta = (x_f - x_i)/Nbins;

    for(int i = 0 ; i <=Nbins ; i++)
    {
        binning[i] = x_i + delta*i;
    }

    return;
}

double get_hwhm(TH1F* h)
{
    // As written, this function works for single peak distributions
    double half_max = h->GetMaximum()/2.;
 
    int halfwidth_bin;
    for(int bin = 1 ; bin <= h->GetNbinsX() ; bin++)
    {
        if (h->GetBinContent(bin)>=half_max){halfwidth_bin = bin; break;}
    }
    
    // Return the Half width at half maximum value
    return abs(h->GetBinCenter(h->GetMaximumBin()) - h->GetBinCenter(halfwidth_bin));
}

#endif