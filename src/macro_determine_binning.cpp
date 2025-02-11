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
    TH1F* h_jet_pt = new TH1F("h_jet_pt","",100000,jet_pt_min,jet_pt_max);
    TH1F* h_weight = new TH1F("h_weight","",100000,0,.4);
    TH1F* h_jet_e  = new TH1F("h_jet_e" ,"",100000,100,3000);

    // Put the data into the histos and get the cumulative distributions
    ntuple->Project("h_jet_pt","jet_pt",pair_cut);
    ntuple->Project("h_weight","weight",pair_cut);
    ntuple->Project("h_jet_e" ,"jet_e" ,pair_cut);

    TH1F* h_jet_pt_cumul = (TH1F*) h_jet_pt->GetCumulative();
    TH1F* h_jet_e_cumul  = (TH1F*) h_jet_e->GetCumulative();
    TH1F* h_weight_cumul = (TH1F*) h_weight->GetCumulative();

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

// Determine binning in jet pt
    double entries_bin_jet_e = h_jet_e->Integral()/Nbin_jet_e;
    std::cout<<"The number of entries for jet_e is "<<entries_bin_jet_e<<std::endl;
    std::cout<<"Binning in jet_e : {jet_e_min";
    int counter_jet_e = 1;
    for(int bin = 1 ; bin <= h_jet_e->GetNbinsX() ; bin++)
    {
        double q = h_jet_e_cumul->GetBinContent(bin);

        // Exit when determined last bin
        if(counter_jet_e==Nbin_jet_e) 
        {
            std::cout<<", jet_e_max}"<<std::endl;
            break;
        }
        // Condition to determine limit
        if(q>entries_bin_jet_e*counter_jet_e) 
        {
            counter_jet_e++;

            // Print limits
            std::cout<<", "<<h_jet_e_cumul->GetBinCenter(bin-1);
        }
    }

    // Determine binning in weights
    double entries_bin_weight = h_weight->Integral()/Nbin_weight;
    std::cout<<"The number of entries for weight is "<<entries_bin_weight<<std::endl;
    std::cout<<"Binning in weight : {0";
    int counter_weight = 1;
    for(int bin = 1 ; bin <= h_weight->GetNbinsX() ; bin++)
    {
        double q = h_weight_cumul->GetBinContent(bin);

        // Exit when determined last bin
        if(counter_weight==Nbin_weight) 
        {
            std::cout<<", 0.4}"<<std::endl;
            break;
        }
        // Condition to determine limit
        if(q>entries_bin_weight*counter_weight) 
        {
            counter_weight++;

            // Print limits
            std::cout<<", "<<h_weight_cumul->GetBinCenter(bin-1);
        }
    }
    
    // Close file
    f->Close();
    double binning[Nbin_R_L+1];
    determine_eqsizebinning(Nbin_R_L, R_L_min, R_L_max, binning);

    std::cout<<"Binning in R_L : {R_L_min";
    for(int i = 1 ; i < Nbin_R_L ; i++)
    {
        std::cout<<", "<<binning[i];
    }
    std::cout<<", R_L_max};"<<std::endl;

    return 0;
}