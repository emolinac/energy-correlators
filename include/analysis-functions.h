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
        // This is the case in which both particles are in different sides of the plane
        double h1h2_angle = TMath::Abs(h1_phi) + TMath::Abs(h2_phi);

        if(h1h2_angle>TMath::Pi()) delta_phi = 2*TMath::Pi() - h1h2_angle;
        else delta_phi = h1h2_angle;
    }
    
    return sqrt(delta_eta*delta_eta + delta_phi*delta_phi);
}

#endif