#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"

void macro_print_responsematrix_weights()
{
  TFile* f = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());
  TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

  float weight, weight_truth;
  ntuple->SetBranchAddress("weight",&weight);
  ntuple->SetBranchAddress("weight_truth",&weight_truth);

  // Create histograms with the respective true and matched reco 
//   double em_weightunfolding_binning[Nbin_weight+1];
//   determine_log10binning(Nbin_weight,weight_min,weight_max,em_weightunfolding_binning);    
//   TH2F* hresp = new TH2F("hresp","",Nbin_weight, em_weightunfolding_binning,Nbin_weight, em_weightunfolding_binning);

  TH2F* hresp = new TH2F("hresp","",Nbin_weight,weight_binning,Nbin_weight,weight_binning);

  for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
  {
    // Access entry of ntuple
    ntuple->GetEntry(evt);
    if(weight_truth!=-999) hresp->Fill(weight, weight_truth);
  }

  // Draw response matrix
  hresp->Draw("COL TEXT");
  hresp->SetTitle("Response matrix of weight;Reco;Truth");
  gPad->SetLogy(1);
  gPad->SetLogx(1);
}