#ifndef ANALYSIS_FUNCTIONS_H
#define ANALYSIS_FUNCTIONS_H

double weight(double h1_E, double h2_E, double jet_E)
{
    return h1_E*h2_E/jet_E/jet_E;
}

double X_L(double h1_eta, double h2_eta, double h1_phi, double h2_phi)
{
    double delta_eta = h1_eta - h2_eta;
    double delta_phi = 0;

    if(h1_phi*h2_phi>0) delta_phi = h1_phi - h2_phi;
    else
    {
        if(TMath::Abs(h1_phi - h2_phi)>TMath::Pi()) delta_phi = 2*TMath::Pi() - TMath::Abs(h1_phi - h2_phi);
        else delta_phi = h1_phi - h2_phi;
    }
    
    return sqrt(delta_eta*delta_eta + delta_phi*delta_phi);
}

#endif