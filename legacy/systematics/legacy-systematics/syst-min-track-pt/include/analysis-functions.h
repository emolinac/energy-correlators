#ifndef ANALYSIS_FUNCTIONS_H
#define ANALYSIS_FUNCTIONS_H

double weight(double h1_E, double h2_E, double jet_E)
{
    return h1_E*h2_E/jet_E/jet_E;
}

double R_L(double h1_y, double h2_y, double h1_phi, double h2_phi)
{
    double delta_y = h1_y - h2_y;
    double delta_phi = 0;

    if (h1_phi*h2_phi>0) delta_phi = h1_phi - h2_phi;
    else
    {
        // This is the case in which both particles are in different sides of the plane
        double h1h2_angle = TMath::Abs(h1_phi) + TMath::Abs(h2_phi);

        if (h1h2_angle>TMath::Pi()) delta_phi = 2*TMath::Pi() - h1h2_angle;
        else delta_phi = h1h2_angle;
    }
    
    return sqrt(delta_y*delta_y + delta_phi*delta_phi);
}

double delta_phi(double h1_phi, double h2_phi)
{
    double delta_phi = 0;

    if (h1_phi*h2_phi>0) delta_phi = h1_phi - h2_phi;
    else
    {
        // This is the case in which both particles are in different sides of the plane
        double h1h2_angle = TMath::Abs(h1_phi) + TMath::Abs(h2_phi);

        if (h1h2_angle>TMath::Pi()) delta_phi = 2*TMath::Pi() - h1h2_angle;
        else delta_phi = h1h2_angle;
    }
    
    return delta_phi;
}

double rapidity(double energy, double pz)
{
   return 0.5*log( (energy+pz) / (energy-pz) );
}

#endif