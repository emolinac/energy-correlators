#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_determine_binning()
{
    // Open data file 
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());
    
    // Get TNtuple from file
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_data.c_str());
    
    // Declare the histos to use
    TH1F* h_jet_pt = new TH1F("h_jet_pt","",10000,jet_pt_min,jet_pt_max);

    // Put the data into the histos and get the cumulative distributions
    ntuple->Project("h_jet_pt","jet_pt",pair_data_cut);

    TH1F* h_jet_pt_cumul = (TH1F*) h_jet_pt->GetCumulative();

    // Determine binning in jet pt
    double entries_bin_jet_pt = h_jet_pt->Integral()/Nbin_jet_pt;
    std::cout<<"The number of entries for jet_pt is "<<entries_bin_jet_pt<<std::endl;
    std::cout<<"Binning in jet_pt : {jet_pt_min";
    int counter_jet_pt = 1;
    for(int bin = 1 ; bin <= h_jet_pt->GetNbinsX() ; bin++)
    {
        double q = h_jet_pt_cumul->GetBinContent(bin);

        // Exit when determined last bin
        if(counter_jet_pt==Nbin_jet_pt) 
        {
            std::cout<<", jet_pt_max}"<<std::endl;
            break;
        }
        // Condition to determine limit
        if(q>entries_bin_jet_pt*counter_jet_pt) 
        {
            counter_jet_pt++;

            // Print limits
            std::cout<<", "<<h_jet_pt_cumul->GetBinCenter(bin-1);
        }
    }
    
    // Close file
    f->Close();

    std::cout<<"Given a range of angular distance from "<<R_L_min<<" to "<<R_L_max<<std::endl;
    std::cout<<"Given an angular resolution of "<<R_L_res<<" we would need "<<(R_L_max-R_L_min)/R_L_res<<std::endl;

    return 0;
}