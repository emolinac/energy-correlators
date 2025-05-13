#include <iostream>
#include "TZJetsMC.h"
#include "TZJetsMC.C"
#include "TZJetsMCReco.h"
#include "TZJetsMCReco.C"
#include "TZJetsData.h"
#include "TZJetsData.C"
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH3.h"
#include "analysis-constants.h"
#include "analysis-binning.h"
#include "analysis-cuts.h"
#include "analysis-functions.h"
#include "directories.h"
#include "names.h"

int main()
{
  // Open correction files
  TFile* fpurity     = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());
  TFile* fefficiency = new TFile((output_folder+namef_ntuple_e2c_efficiency).c_str());
  
  // Create output file
  TFile* fout = new TFile((output_folder+namef_ntuple_e2c_corr).c_str(),"RECREATE");
  
  // Declare the TTrees to be used to build the ntuples
  TZJetsData* datatree = new TZJetsData();
  
  // Create Ntuples
  TNtuple* ntuple_purity          = (TNtuple*) fpurity->Get((name_ntuple_purity.c_str()));
  TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_mc.c_str()));
  TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_reco.c_str()));
  TNtuple* ntuple_data            = new TNtuple(name_ntuple_data.c_str(),"All Data",ntuple_corrdata_vars); 
  ntuple_data->SetAutoSave(0);

  // Calculate corrections
  TH3F* hsigp   = new TH3F("hsigp"  ,"",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,Nbin_jetpt_corrections,corrections_jetpt_binning);
  TH3F* hallp   = new TH3F("hallp"  ,"",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,Nbin_jetpt_corrections,corrections_jetpt_binning);
  TH3F* hpurity = new TH3F("hpurity","",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,Nbin_jetpt_corrections,corrections_jetpt_binning);
  hsigp->Sumw2();
  hallp->Sumw2();
  hpurity->Sumw2();

  TH3F* hsigeff     = new TH3F("hsigeff"    ,"",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,Nbin_jetpt_corrections,corrections_jetpt_binning);
  TH3F* halleff     = new TH3F("halleff"    ,"",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,Nbin_jetpt_corrections,corrections_jetpt_binning);
  TH3F* hefficiency = new TH3F("hefficiency","",ic_p_nbins,ic_p_binning,sl_eta_nbins,sl_eta_binning,Nbin_jetpt_corrections,corrections_jetpt_binning);
  hsigeff->Sumw2();
  halleff->Sumw2();
  hefficiency->Sumw2();

  ntuple_purity->Project("hsigp","jet_pt:h_eta:h_p",single_signal_cut);
  ntuple_purity->Project("hallp","jet_pt:h_eta:h_p",pair_cut         );
  hpurity->Divide(hsigp,hallp,1,1,"B");

  // ntuple_efficiency_reco->Project("hsigeff","jet_pt:h_eta:h_p",single_signal_cut);
  ntuple_efficiency_reco->Project("hsigeff","jet_pt_truth:h_eta_truth:h_p_truth",single_signal_cut);
  ntuple_efficiency_mc->Project("halleff","jet_pt:h_eta:h_p",pair_cut         );
  hefficiency->Divide(hsigeff,halleff,1,1,"B");

  // Create necessary 4vectors
  TLorentzVector* Jet_4vector = new TLorentzVector();
  TLorentzVector* Z0_4vector  = new TLorentzVector();
  TLorentzVector* mum_4vector = new TLorentzVector();
  TLorentzVector* mup_4vector = new TLorentzVector();
  TLorentzVector* h1_4vector  = new TLorentzVector();
  TLorentzVector* h2_4vector  = new TLorentzVector();
  int eventNum;
  unsigned long long last_eventNum = 0;
  int events = 0;
  bool maxjetpT_found = false;
    
  // Define array carrying the variables
  float vars[Nvars_corrdata];

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
        
    // Loop over hadron 1
    for(int h1_index = 0 ; h1_index < datatree->Jet_NDtr ; h1_index++)
    {
      // Skip non-hadronic particles
      if(datatree->Jet_Dtr_IsMeson[h1_index]!=1&&datatree->Jet_Dtr_IsBaryon[h1_index]!=1) continue;

      h1_4vector->SetPxPyPzE(datatree->Jet_Dtr_PX[h1_index]/1000.,datatree->Jet_Dtr_PY[h1_index]/1000.,datatree->Jet_Dtr_PZ[h1_index]/1000.,datatree->Jet_Dtr_E[h1_index]/1000.);
      if(!apply_chargedtrack_cuts(datatree->Jet_Dtr_ThreeCharge[h1_index],
                                  datatree->Jet_Dtr_P[h1_index]/1000.,
                                  datatree->Jet_Dtr_PT[h1_index]/1000.,
                                  datatree->Jet_Dtr_TrackChi2[h1_index]/datatree->Jet_Dtr_TrackNDF[h1_index],
                                  datatree->Jet_Dtr_ProbNNghost[h1_index],
                                  Jet_4vector->DeltaR(*h1_4vector))) continue;

      double h1_purity     = hpurity->GetBinContent(hpurity->FindBin(datatree->Jet_Dtr_P[h1_index]/1000.,datatree->Jet_Dtr_ETA[h1_index],datatree->Jet_PT/1000.));
      double h1_efficiency = hefficiency->GetBinContent(hefficiency->FindBin(datatree->Jet_Dtr_P[h1_index]/1000.,datatree->Jet_Dtr_ETA[h1_index],datatree->Jet_PT/1000.));
      if(h1_purity>1.||h1_efficiency>1.) continue;

      // Loop over hadron 2
      for(int h2_index = h1_index+1 ; h2_index < datatree->Jet_NDtr ; h2_index++)
      {
        // Skip non-hadronic particles
        if(datatree->Jet_Dtr_IsMeson[h2_index]!=1&&datatree->Jet_Dtr_IsBaryon[h2_index]!=1) continue;

        h2_4vector->SetPxPyPzE(datatree->Jet_Dtr_PX[h2_index]/1000.,datatree->Jet_Dtr_PY[h2_index]/1000.,datatree->Jet_Dtr_PZ[h2_index]/1000.,datatree->Jet_Dtr_E[h2_index]/1000.);
        if(!apply_chargedtrack_cuts(datatree->Jet_Dtr_ThreeCharge[h2_index],
                                    datatree->Jet_Dtr_P[h2_index]/1000.,
                                    datatree->Jet_Dtr_PT[h2_index]/1000.,
                                    datatree->Jet_Dtr_TrackChi2[h2_index]/datatree->Jet_Dtr_TrackNDF[h2_index],
                                    datatree->Jet_Dtr_ProbNNghost[h2_index],
                                    Jet_4vector->DeltaR(*h2_4vector))) continue;

        double h2_purity     = hpurity->GetBinContent(hpurity->FindBin(datatree->Jet_Dtr_P[h2_index]/1000.,datatree->Jet_Dtr_ETA[h2_index],datatree->Jet_PT/1000.));
        double h2_efficiency = hefficiency->GetBinContent(hefficiency->FindBin(datatree->Jet_Dtr_P[h2_index]/1000.,datatree->Jet_Dtr_ETA[h2_index],datatree->Jet_PT/1000.));
        if(h2_purity>1.||h2_efficiency>1.) continue;

        double purity_correction = (h1_purity)*(h2_purity);

        double h1_purity_err = hpurity->GetBinError(hpurity->FindBin(datatree->Jet_Dtr_P[h1_index]/1000.,datatree->Jet_Dtr_ETA[h1_index],datatree->Jet_PT/1000.));
        double h2_purity_err = hpurity->GetBinError(hpurity->FindBin(datatree->Jet_Dtr_P[h2_index]/1000.,datatree->Jet_Dtr_ETA[h2_index],datatree->Jet_PT/1000.));
        double purity_error = sqrt(pow((h1_purity)*(h2_purity_err),2) + pow((h1_purity_err)*(h2_purity),2));

        double efficiency_correction = (h1_efficiency)*(h2_efficiency);

        double h1_efficiency_err = hefficiency->GetBinError(hefficiency->FindBin(datatree->Jet_Dtr_P[h1_index]/1000.,datatree->Jet_Dtr_ETA[h1_index],datatree->Jet_PT/1000.));
        double h2_efficiency_err = hefficiency->GetBinError(hefficiency->FindBin(datatree->Jet_Dtr_P[h2_index]/1000.,datatree->Jet_Dtr_ETA[h2_index],datatree->Jet_PT/1000.));
        double efficiency_error = sqrt(pow((h1_efficiency)*(h2_efficiency_err),2) + pow((h1_efficiency_err)*(h2_efficiency),2));

        vars[0 ] = weight(datatree->Jet_Dtr_E[h1_index], datatree->Jet_Dtr_E[h2_index], datatree->Jet_PE);
        vars[1 ] = efficiency_correction;
        vars[2 ] = purity_correction;
        vars[3 ] = efficiency_error/efficiency_correction;
        vars[4 ] = purity_error/purity_correction;
        vars[5 ] = R_L(rapidity(datatree->Jet_Dtr_E[h1_index],datatree->Jet_Dtr_PZ[h1_index]), rapidity(datatree->Jet_Dtr_E[h2_index],datatree->Jet_Dtr_PZ[h2_index]),
                                datatree->Jet_Dtr_PHI[h1_index], datatree->Jet_Dtr_PHI[h2_index]);
        vars[6 ] = datatree->Jet_Dtr_ETA[h1_index];
        vars[7 ] = datatree->Jet_Dtr_ETA[h2_index];
        vars[8 ] = rapidity(datatree->Jet_Dtr_E[h1_index],datatree->Jet_Dtr_PZ[h1_index]);
        vars[9 ] = rapidity(datatree->Jet_Dtr_E[h2_index],datatree->Jet_Dtr_PZ[h2_index]);
        vars[10] = datatree->Jet_Dtr_PHI[h1_index];
        vars[11] = datatree->Jet_Dtr_PHI[h2_index];
        vars[12] = datatree->Jet_Dtr_P[h1_index]/1000.;
        vars[13] = datatree->Jet_Dtr_P[h2_index]/1000.;
        vars[14] = datatree->Jet_Dtr_PT[h1_index]/1000.;
        vars[15] = datatree->Jet_Dtr_PT[h2_index]/1000.;
        vars[16] = datatree->Jet_PT/1000.;
        vars[17] = Jet_4vector->Eta();
        vars[18] = Jet_4vector->Phi();
        vars[19] = delta_phi(Jet_4vector->Phi(),Z0_4vector->Phi());
        vars[20] = datatree->Jet_PE/1000.;//R_L(Jet_4vector->Rapidity(),mum_4vector->Rapidity(),Jet_4vector->Phi(),mum_4vector->Phi());
        vars[21] = datatree->Jet_Dtr_E[h1_index]/1000.;//mum_4vector->Pt();
        vars[22] = datatree->Jet_Dtr_E[h2_index]/1000.;//mum_4vector->Eta();
        //vars[23] = //R_L(Jet_4vector->Rapidity(),mup_4vector->Rapidity(),Jet_4vector->Phi(),mup_4vector->Phi());
        //vars[24] = //mup_4vector->Pt();
        //vars[25] = //mup_4vector->Eta();
        
        // Fill the TNtuple
        ntuple_data->Fill(vars);
      }
    }
  }

  fout->cd();
  ntuple_data->Write();
  fout->Close();
  
  std::cout<<std::endl;

  return 0;
}
