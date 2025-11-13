#ifndef ANALYSIS_FUNCTIONS_H
#define ANALYSIS_FUNCTIONS_H

double weight(double h1_E, double h2_E, double jet_E)
{
        return h1_E*h2_E/jet_E/jet_E;
}

#endif