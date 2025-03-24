#include <iostream>
#include "TZJetsData.h"
#include "TZJetsData.C"
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH3.h"
#include "analysis-constants.h"
#include "analysis-cuts.h"
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"

int main()
{
  // Open correction files
  TFile* fid  = new TFile((muons_folder+"IDEff_Data_2016.root").c_str());
  TFile* ftrg = new TFile((muons_folder+"TRGEff_Data_2016.root").c_str());
  TFile* ftrk = new TFile((muons_folder+"TRKEff_Data_2016.root").c_str());
  
  // Create output file
  TFile* fout = new TFile((output_folder+namef_ntuple_muon).c_str(),"RECREATE");
  
  // Declare the TTrees to be used to build the ntuples
  TZJetsData* datatree = new TZJetsData();
  
  // Create Ntuples
  TNtuple* ntuple_muons = new TNtuple(name_ntuple_muon.c_str(),"",ntuple_muons_vars);
  ntuple_muons->SetAutoSave(0);

  TH2D* h2_muon_ideff_data  = (TH2D*) fid->Get("Hist_ALL_2016_ETA_PT_Eff");
  TH2D* h2_muon_trkeff_data = (TH2D*) ftrg->Get("Hist_ALL_2016_ETA_PT_Eff");
  TH2D* h2_muon_trgeff_data = (TH2D*) ftrk->Get("Hist_ALL_2016_ETA_PT_Eff");

  // Create necessary 4vectors
  TLorentzVector* Jet_4vector = new TLorentzVector();
  TLorentzVector* Z0_4vector  = new TLorentzVector();
  TLorentzVector* mum_4vector = new TLorentzVector();
  TLorentzVector* mup_4vector = new TLorentzVector();
  
  int eventNum;
  unsigned long long last_eventNum = 0;
  int events = 0;
  bool maxjetpT_found = false;
    
  // Define array carrying the variables
  float vars[Nvars_muons];
  
  // Fill the data TNtuple
  for(int evt = 0 ; evt < datatree->fChain->GetEntries() ; evt++)
  {
    // Access entry of tree
    datatree->GetEntry(evt);
    if(evt%10000==0)
    {
      double percentage = 100.*evt/datatree->fChain->GetEntries();
      std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
    }
    // Access entry of tree
    datatree->GetEntry(evt);

    if (evt != 0)
    {
      if (datatree->eventNumber != last_eventNum) maxjetpT_found = false;
      if (last_eventNum == datatree->eventNumber) continue;
    }

    last_eventNum = datatree->eventNumber;
    if (maxjetpT_found) continue;

    // Apply PV cut
    if(datatree->nPV!=1) continue;

    // Apply trigger cut
    bool mum_trigger = (datatree->mum_L0MuonEWDecision_TOS==1&&datatree->mum_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree->mum_Hlt2EWSingleMuonVHighPtDecision_TOS==1);
    bool mup_trigger = (datatree->mup_L0MuonEWDecision_TOS==1&&datatree->mup_Hlt1SingleMuonHighPTDecision_TOS==1&&datatree->mup_Hlt2EWSingleMuonVHighPtDecision_TOS==1);

    if(!mum_trigger&&!mup_trigger) continue;
    
    // Set Jet-associated 4 vectors and apply cuts
    Jet_4vector->SetPxPyPzE(datatree->Jet_PX/1000.,datatree->Jet_PY/1000.,datatree->Jet_PZ/1000.,datatree->Jet_PE/1000.);
    if(!apply_jet_cuts(Jet_4vector->Eta(),Jet_4vector->Pt())) continue;
    
    mum_4vector->SetPxPyPzE(datatree->mum_PX/1000.,datatree->mum_PY/1000.,datatree->mum_PZ/1000.,datatree->mum_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector),mum_4vector->Pt(),mum_4vector->Eta())) continue;
    //if(datatree->mum_TRACK_PCHI2<muon_trackprob_min) continue;
    
    mup_4vector->SetPxPyPzE(datatree->mup_PX/1000.,datatree->mup_PY/1000.,datatree->mup_PZ/1000.,datatree->mup_PE/1000.);
    if(!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector),mup_4vector->Pt(),mup_4vector->Eta())) continue;
    //if(datatree->mup_TRACK_PCHI2<muon_trackprob_min) continue;
    
    Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),mup_4vector->Py()+mum_4vector->Py(),mup_4vector->Pz()+mum_4vector->Pz(),mup_4vector->E() +mum_4vector->E());
    if(!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)),Z0_4vector->M())) continue;

    double mup_eff_id  = h2_muon_ideff_data->GetBinContent(h2_muon_ideff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));
    double mup_eff_trk = h2_muon_trkeff_data->GetBinContent(h2_muon_trkeff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));
    double mup_eff_trg = h2_muon_trgeff_data->GetBinContent(h2_muon_trgeff_data->FindBin(mup_4vector->Eta(),mup_4vector->Pt()));

    double mum_eff_id  = h2_muon_ideff_data->GetBinContent(h2_muon_ideff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));
    double mum_eff_trk = h2_muon_trkeff_data->GetBinContent(h2_muon_trkeff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));
    double mum_eff_trg = h2_muon_trgeff_data->GetBinContent(h2_muon_trgeff_data->FindBin(mum_4vector->Eta(),mum_4vector->Pt()));

    vars[0] = mup_4vector->Eta();
    vars[1] = mup_4vector->Pt();
    vars[2] = mup_eff_id;
    vars[3] = mup_eff_trk;
    vars[4] = mup_eff_trg;
    vars[5] = mum_4vector->Eta();
    vars[6] = mum_4vector->Pt();
    vars[7] = mum_eff_id;
    vars[8] = mum_eff_trk;
    vars[9] = mum_eff_trg;
    
    ntuple_muons->Fill(vars);
  }

  fout->cd();
  ntuple_muons->Write();
  fout->Close();
  
  std::cout<<std::endl;

  return 0;
}
