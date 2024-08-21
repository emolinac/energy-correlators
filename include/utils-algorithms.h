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

#endif