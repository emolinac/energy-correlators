#include <iostream>
#include "TZJetsMC.h"
#include "TZJetsMC.C"
#include "TZJetsMCReco.h"
#include "TZJetsMCReco.C"
#include "TZJets2016Data.h"
#include "TZJets2016Data.C"
#include "TZJets2017Data.h"
#include "TZJets2017Data.C"
#include "TZJets2018Data.h"
#include "TZJets2018Data.C"
#include "TROOT.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TH3.h"
#include "TLatex.h"
#include "analysis-constants.h"
#include "analysis-binning.h"
#include "analysis-cuts.h"
#include "directories.h"
#include "names.h"
#include "syst-jes-jer.h"
#include "utils.h"

int main(int argc, char* argv[])
{
        bool get_nominal = false;
        bool get_jes_jer = false;
        bool get_jer     = false;
        bool get_muon    = false;

        if (argc < 2 || std::string(argv[1]) == "--help") {
                std::cout<<"You have to pass one argument to this code. Possible options are:"<<std::endl;
                std::cout<<"--get-nominal : Gets the nominal pair corrections."<<std::endl;
                std::cout<<"--get-jes-jer : Gets the pair corrections with the JES-JER variation."<<std::endl;
                std::cout<<"--get-jer     : Gets the pair corrections with the JER variation."<<std::endl;
                std::cout<<"--get-muon    : Gets the pair corrections with the muon eff varied."<<std::endl;

                return 0;
        }

        if (std::string(argv[1]) == "--get-nominal")
                get_nominal = true;
        else if (std::string(argv[1]) == "--get-jes-jer")
                get_jes_jer = true;
        else if (std::string(argv[1]) == "--get-jer")
                get_jer = true;
        else if (std::string(argv[1]) == "--get-muon")
                get_muon = true;
        else {
                std::cout<<"No valid options provided"<<std::endl;
                return 0;
        }
        
        // Create output file
        TFile* fout = new TFile((output_folder + namef_all_corrections[argv[1]]).c_str(),"RECREATE");
        gROOT->cd();
        
        // Declare the TTrees to be used to build the ntuples
        TZJets2016Data* datatree_2016 = new TZJets2016Data();
        TZJets2017Data* datatree_2017 = new TZJets2017Data();
        TZJets2018Data* datatree_2018 = new TZJets2018Data();
        
        // Open correction files
        TFile* fcorrections_pair = new TFile((output_folder + namef_reco_corrections[argv[1]]).c_str());
        TFile* fefficiency_jet   = new TFile((output_folder + namef_ntuple_truth2reco_match).c_str());
        
        TFile* fefficiency_muon_2016_id  = new TFile((muons_folder + "IDEff_Data_2016.root").c_str());
        TFile* fefficiency_muon_2016_trk = new TFile((muons_folder + "TRKEff_Data_2016.root").c_str());
        TFile* fefficiency_muon_2016_trg = new TFile((muons_folder + "TRGEff_Data_2016.root").c_str());
        TFile* fefficiency_muon_2017_id  = new TFile((muons_folder + "IDEff_Data_2017.root").c_str());
        TFile* fefficiency_muon_2017_trk = new TFile((muons_folder + "TRKEff_Data_2017.root").c_str());
        TFile* fefficiency_muon_2017_trg = new TFile((muons_folder + "TRGEff_Data_2017.root").c_str());
        TFile* fefficiency_muon_2018_id  = new TFile((muons_folder + "IDEff_Data_2018.root").c_str());
        TFile* fefficiency_muon_2018_trk = new TFile((muons_folder + "TRKEff_Data_2018.root").c_str());
        TFile* fefficiency_muon_2018_trg = new TFile((muons_folder + "TRGEff_Data_2018.root").c_str());
        
        // Create Ntuples
        TNtuple* ntuple_purity          = (TNtuple*) fcorrections_pair->Get((name_ntuple_correction_reco.c_str()));
        TNtuple* ntuple_efficiency_mc   = (TNtuple*) fcorrections_pair->Get((name_ntuple_correction_mc.c_str()));
        TNtuple* ntuple_efficiency_reco = (TNtuple*) fcorrections_pair->Get((name_ntuple_correction_reco.c_str()));
        TNtuple* ntuple_purity_jet      = (TNtuple*) fcorrections_pair->Get((name_ntuple_jet_reco2truth_match.c_str()));
        TNtuple* ntuple_efficiency_jet  = (TNtuple*) fefficiency_jet->Get((name_ntuple_jet_truth2reco_match.c_str()));
        
        // Muon corrections
        TH2D* h2_muon_2016_ideff_data  = (TH2D*) fefficiency_muon_2016_id->Get("Hist_ALL_2016_ETA_PT_Eff");
        TH2D* h2_muon_2016_trkeff_data = (TH2D*) fefficiency_muon_2016_trk->Get("Hist_ALL_2016_ETA_PT_Eff");
        TH2D* h2_muon_2016_trgeff_data = (TH2D*) fefficiency_muon_2016_trg->Get("Hist_ALL_2016_ETA_PT_Eff");
        TH2D* h2_muon_2017_ideff_data  = (TH2D*) fefficiency_muon_2017_id->Get("Hist_ALL_2017_ETA_PT_Eff");
        TH2D* h2_muon_2017_trkeff_data = (TH2D*) fefficiency_muon_2017_trk->Get("Hist_ALL_2017_ETA_PT_Eff");
        TH2D* h2_muon_2017_trgeff_data = (TH2D*) fefficiency_muon_2017_trg->Get("Hist_ALL_2017_ETA_PT_Eff");
        TH2D* h2_muon_2018_ideff_data  = (TH2D*) fefficiency_muon_2018_id->Get("Hist_ALL_2018_ETA_PT_Eff");
        TH2D* h2_muon_2018_trkeff_data = (TH2D*) fefficiency_muon_2018_trk->Get("Hist_ALL_2018_ETA_PT_Eff");
        TH2D* h2_muon_2018_trgeff_data = (TH2D*) fefficiency_muon_2018_trg->Get("Hist_ALL_2018_ETA_PT_Eff");

        // Jet corrections
        TH1F* hnum_pur_jet = new TH1F("hnum_pur_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* hden_pur_jet = new TH1F("hden_pur_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* hpurity_jet  = new TH1F("hpurity_jet" , "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        hnum_pur_jet->Sumw2();
        hden_pur_jet->Sumw2();

        TH1F* hnum_eff_jet    = new TH1F("hnum_eff_jet"   , "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* hden_eff_jet    = new TH1F("hden_eff_jet"   , "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* hefficiency_jet = new TH1F("hefficiency_jet", "", nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        hnum_eff_jet->Sumw2();
        hden_eff_jet->Sumw2();

        ntuple_purity_jet->Project("hnum_pur_jet", "jet_pt", "jet_pt_truth!=-999");
        ntuple_purity_jet->Project("hden_pur_jet", "jet_pt");
        ntuple_efficiency_jet->Project("hnum_eff_jet", "jet_pt_truth", "jet_pt!=-999");
        ntuple_efficiency_jet->Project("hden_eff_jet", "jet_pt_truth");

        hpurity_jet->Divide(hnum_pur_jet, hden_pur_jet, 1, 1, "B");
        hefficiency_jet->Divide(hnum_eff_jet, hden_eff_jet, 1, 1, "B");

        // Pair corrections
        TH3F* hnum_pur    = new TH3F("hnum_pur"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hden_pur    = new TH3F("hden_pur"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hpurity     = new TH3F("hpurity"    , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hnum_eff    = new TH3F("hnum_eff"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hden_eff    = new TH3F("hden_eff"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hefficiency = new TH3F("hefficiency", "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        
        hnum_pur->Sumw2();
        hden_pur->Sumw2();
        hnum_eff->Sumw2();
        hden_eff->Sumw2();
        
        ntuple_purity->Project("hnum_pur", "weight_pt:jet_pt:R_L", pair_matching_cut);
        ntuple_purity->Project("hden_pur", "weight_pt:jet_pt:R_L");
        ntuple_efficiency_reco->Project("hnum_eff", "weight_pt_truth:jet_pt_truth:R_L_truth", pair_matching_cut);
        ntuple_efficiency_mc->Project("hden_eff", "weight_pt:jet_pt:R_L");

        hpurity->Divide(hnum_pur, hden_pur, 1, 1, "B");
        hefficiency->Divide(hnum_eff, hden_eff, 1, 1, "B");

        regularize_correction_factors(hpurity);
        regularize_correction_factors(hefficiency);

        // Charged config pair corrections
        TH3F* hnum_pur_eqcharge    = new TH3F("hnum_pur_eqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hden_pur_eqcharge    = new TH3F("hden_pur_eqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hpurity_eqcharge     = new TH3F("hpurity_eqcharge"    , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hnum_eff_eqcharge    = new TH3F("hnum_eff_eqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hden_eff_eqcharge    = new TH3F("hden_eff_eqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hefficiency_eqcharge = new TH3F("hefficiency_eqcharge", "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        
        hnum_pur_eqcharge->Sumw2();
        hden_pur_eqcharge->Sumw2();
        hnum_eff_eqcharge->Sumw2();
        hden_eff_eqcharge->Sumw2();
        
        ntuple_purity->Project("hnum_pur_eqcharge", "weight_pt:jet_pt:R_L", pair_matching_cut + "eq_charge > 0");
        ntuple_purity->Project("hden_pur_eqcharge", "weight_pt:jet_pt:R_L", "eq_charge > 0");
        ntuple_efficiency_reco->Project("hnum_eff_eqcharge", "weight_pt_truth:jet_pt_truth:R_L_truth", pair_matching_cut + "eq_charge > 0");
        ntuple_efficiency_mc->Project("hden_eff_eqcharge", "weight_pt:jet_pt:R_L", "eq_charge > 0");

        hpurity_eqcharge->Divide(hnum_pur_eqcharge, hden_pur_eqcharge, 1, 1, "B");
        hefficiency_eqcharge->Divide(hnum_eff_eqcharge, hden_eff_eqcharge, 1, 1, "B");

        regularize_correction_factors(hpurity_eqcharge);
        regularize_correction_factors(hefficiency_eqcharge);

        TH3F* hnum_pur_neqcharge    = new TH3F("hnum_pur_neqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hden_pur_neqcharge    = new TH3F("hden_pur_neqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hpurity_neqcharge     = new TH3F("hpurity_neqcharge"    , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hnum_eff_neqcharge    = new TH3F("hnum_eff_neqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hden_eff_neqcharge    = new TH3F("hden_eff_neqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* hefficiency_neqcharge = new TH3F("hefficiency_neqcharge", "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        
        hnum_pur_neqcharge->Sumw2();
        hden_pur_neqcharge->Sumw2();
        hnum_eff_neqcharge->Sumw2();
        hden_eff_neqcharge->Sumw2();
        
        ntuple_purity->Project("hnum_pur_neqcharge", "weight_pt:jet_pt:R_L", pair_matching_cut + "eq_charge < 0");
        ntuple_purity->Project("hden_pur_neqcharge", "weight_pt:jet_pt:R_L", "eq_charge < 0");
        ntuple_efficiency_reco->Project("hnum_eff_neqcharge", "weight_pt_truth:jet_pt_truth:R_L_truth", pair_matching_cut + "eq_charge < 0");
        ntuple_efficiency_mc->Project("hden_eff_neqcharge", "weight_pt:jet_pt:R_L", "eq_charge < 0");

        hpurity_neqcharge->Divide(hnum_pur_neqcharge, hden_pur_neqcharge, 1, 1, "B");
        hefficiency_neqcharge->Divide(hnum_eff_neqcharge, hden_eff_neqcharge, 1, 1, "B");

        regularize_correction_factors(hpurity_neqcharge);
        regularize_correction_factors(hefficiency_neqcharge);

        // Print the corrections effect as a of function RL and jet pt
        if (get_nominal) {
                TCanvas* c = new TCanvas("c", "", 1920, 1080);
                c->Draw();

                gStyle->SetOptStat("");
                gStyle->SetPaintTextFormat("4.2f");        
                
                TLatex latex;
                latex.SetTextAlign(22); // center alignment
                latex.SetTextSize(text_size_correction_plots);
                latex.SetTextColor(kBlack);

                // Inclusive
                TH2F* hnum_pur_rl_jetpt    = new TH2F("hnum_pur_rl_jetpt"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hden_pur_rl_jetpt    = new TH2F("hden_pur_rl_jetpt"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hpurity_rl_jetpt     = new TH2F("hpurity_rl_jetpt"    , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hnum_eff_rl_jetpt    = new TH2F("hnum_eff_rl_jetpt"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hden_eff_rl_jetpt    = new TH2F("hden_eff_rl_jetpt"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hefficiency_rl_jetpt = new TH2F("hefficiency_rl_jetpt", "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                
                hnum_pur_rl_jetpt->Sumw2();
                hden_pur_rl_jetpt->Sumw2();
                hnum_eff_rl_jetpt->Sumw2();
                hden_eff_rl_jetpt->Sumw2();
                
                ntuple_purity->Project("hnum_pur_rl_jetpt", "jet_pt:R_L", pair_matching_cut);
                ntuple_purity->Project("hden_pur_rl_jetpt", "jet_pt:R_L");
                ntuple_efficiency_reco->Project("hnum_eff_rl_jetpt", "jet_pt_truth:R_L_truth", pair_matching_cut);
                ntuple_efficiency_mc->Project("hden_eff_rl_jetpt", "jet_pt:R_L");

                hpurity_rl_jetpt->Divide(hnum_pur_rl_jetpt, hden_pur_rl_jetpt, 1, 1, "B");
                hefficiency_rl_jetpt->Divide(hnum_eff_rl_jetpt, hden_eff_rl_jetpt, 1, 1, "B");

                hpurity_rl_jetpt->Draw("col");
                for (int i = 2; i < hpurity_rl_jetpt->GetNbinsX(); ++i) {
                        for (int j = 3; j <= hpurity_rl_jetpt->GetNbinsY(); ++j) {
                                double x = hpurity_rl_jetpt->GetXaxis()->GetBinCenter(i);
                                double y = hpurity_rl_jetpt->GetYaxis()->GetBinCenter(j);
                                double content = hpurity_rl_jetpt->GetBinContent(i, j);
                                double error = hpurity_rl_jetpt->GetBinError(i, j);

                                latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                        }
                }
                hpurity_rl_jetpt->SetTitle("Purity Correction;R_{L};p_{T,jet}(GeV)");
                hpurity_rl_jetpt->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                hpurity_rl_jetpt->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                c->Print("../src-analysis/plots/pair_purity_correction_rl_jetpt.pdf");

                hefficiency_rl_jetpt->Draw("col");
                for (int i = 2; i < hefficiency_rl_jetpt->GetNbinsX(); ++i) {
                        for (int j = 3; j <= hefficiency_rl_jetpt->GetNbinsY(); ++j) {
                                double x = hefficiency_rl_jetpt->GetXaxis()->GetBinCenter(i);
                                double y = hefficiency_rl_jetpt->GetYaxis()->GetBinCenter(j);
                                double content = hefficiency_rl_jetpt->GetBinContent(i, j);
                                double error = hefficiency_rl_jetpt->GetBinError(i, j);

                                latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                        }
                }
                hefficiency_rl_jetpt->SetTitle("Efficiency Correction;R_{L};p_{T,jet}(GeV)");
                hefficiency_rl_jetpt->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                hefficiency_rl_jetpt->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                c->Print("../src-analysis/plots/pair_efficiency_correction_rl_jetpt.pdf");

                // Eq. charged EECs
                TH2F* hnum_pur_rl_jetpt_eqcharge    = new TH2F("hnum_pur_rl_jetpt_eqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hden_pur_rl_jetpt_eqcharge    = new TH2F("hden_pur_rl_jetpt_eqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hpurity_rl_jetpt_eqcharge     = new TH2F("hpurity_rl_jetpt_eqcharge"    , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hnum_eff_rl_jetpt_eqcharge    = new TH2F("hnum_eff_rl_jetpt_eqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hden_eff_rl_jetpt_eqcharge    = new TH2F("hden_eff_rl_jetpt_eqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hefficiency_rl_jetpt_eqcharge = new TH2F("hefficiency_rl_jetpt_eqcharge", "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                
                hnum_pur_rl_jetpt_eqcharge->Sumw2();
                hden_pur_rl_jetpt_eqcharge->Sumw2();
                hnum_eff_rl_jetpt_eqcharge->Sumw2();
                hden_eff_rl_jetpt_eqcharge->Sumw2();
                
                ntuple_purity->Project("hnum_pur_rl_jetpt_eqcharge", "jet_pt:R_L", pair_matching_cut + "eq_charge > 0");
                ntuple_purity->Project("hden_pur_rl_jetpt_eqcharge", "jet_pt:R_L", "eq_charge > 0");
                ntuple_efficiency_reco->Project("hnum_eff_rl_jetpt_eqcharge", "jet_pt_truth:R_L_truth", pair_matching_cut + "eq_charge > 0");
                ntuple_efficiency_mc->Project("hden_eff_rl_jetpt_eqcharge", "jet_pt:R_L", "eq_charge > 0");

                hpurity_rl_jetpt_eqcharge->Divide(hnum_pur_rl_jetpt_eqcharge, hden_pur_rl_jetpt_eqcharge, 1, 1, "B");
                hefficiency_rl_jetpt_eqcharge->Divide(hnum_eff_rl_jetpt_eqcharge, hden_eff_rl_jetpt_eqcharge, 1, 1, "B");

                hpurity_rl_jetpt_eqcharge->Draw("col");
                for (int i = 2; i < hpurity_rl_jetpt_eqcharge->GetNbinsX(); ++i) {
                        for (int j = 3; j <= hpurity_rl_jetpt_eqcharge->GetNbinsY(); ++j) {
                                double x = hpurity_rl_jetpt_eqcharge->GetXaxis()->GetBinCenter(i);
                                double y = hpurity_rl_jetpt_eqcharge->GetYaxis()->GetBinCenter(j);
                                double content = hpurity_rl_jetpt_eqcharge->GetBinContent(i, j);
                                double error = hpurity_rl_jetpt_eqcharge->GetBinError(i, j);

                                latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                        }
                }

                hpurity_rl_jetpt_eqcharge->SetTitle("Purity Correction Eq. Charge;R_{L};p_{T,jet}(GeV)");
                hpurity_rl_jetpt_eqcharge->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                hpurity_rl_jetpt_eqcharge->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                c->Print("../src-analysis/plots/pair_purity_correction_rl_jetpt_eqcharge.pdf");

                hefficiency_rl_jetpt_eqcharge->Draw("col");
                for (int i = 2; i < hefficiency_rl_jetpt_eqcharge->GetNbinsX(); ++i) {
                        for (int j = 3; j <= hefficiency_rl_jetpt_eqcharge->GetNbinsY(); ++j) {
                                double x = hefficiency_rl_jetpt_eqcharge->GetXaxis()->GetBinCenter(i);
                                double y = hefficiency_rl_jetpt_eqcharge->GetYaxis()->GetBinCenter(j);
                                double content = hefficiency_rl_jetpt_eqcharge->GetBinContent(i, j);
                                double error = hefficiency_rl_jetpt_eqcharge->GetBinError(i, j);

                                latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                        }
                }
                hefficiency_rl_jetpt_eqcharge->SetTitle("Efficiency Correction Eq. Charge;R_{L};p_{T,jet}(GeV)");
                hefficiency_rl_jetpt_eqcharge->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                hefficiency_rl_jetpt_eqcharge->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                c->Print("../src-analysis/plots/pair_efficiency_correction_rl_jetpt_eqcharge.pdf");

                // Op. charged EECs
                TH2F* hnum_pur_rl_jetpt_neqcharge    = new TH2F("hnum_pur_rl_jetpt_neqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hden_pur_rl_jetpt_neqcharge    = new TH2F("hden_pur_rl_jetpt_neqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hpurity_rl_jetpt_neqcharge     = new TH2F("hpurity_rl_jetpt_neqcharge"    , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hnum_eff_rl_jetpt_neqcharge    = new TH2F("hnum_eff_rl_jetpt_neqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hden_eff_rl_jetpt_neqcharge    = new TH2F("hden_eff_rl_jetpt_neqcharge"   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                TH2F* hefficiency_rl_jetpt_neqcharge = new TH2F("hefficiency_rl_jetpt_neqcharge", "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
                
                hnum_pur_rl_jetpt_neqcharge->Sumw2();
                hden_pur_rl_jetpt_neqcharge->Sumw2();
                hnum_eff_rl_jetpt_neqcharge->Sumw2();
                hden_eff_rl_jetpt_neqcharge->Sumw2();
                
                ntuple_purity->Project("hnum_pur_rl_jetpt_neqcharge", "jet_pt:R_L", pair_matching_cut + "eq_charge < 0");
                ntuple_purity->Project("hden_pur_rl_jetpt_neqcharge", "jet_pt:R_L", "eq_charge < 0");
                ntuple_efficiency_reco->Project("hnum_eff_rl_jetpt_neqcharge", "jet_pt_truth:R_L_truth", pair_matching_cut + "eq_charge < 0");
                ntuple_efficiency_mc->Project("hden_eff_rl_jetpt_neqcharge", "jet_pt:R_L", "eq_charge < 0");

                hpurity_rl_jetpt_neqcharge->Divide(hnum_pur_rl_jetpt_neqcharge, hden_pur_rl_jetpt_neqcharge, 1, 1, "B");
                hefficiency_rl_jetpt_neqcharge->Divide(hnum_eff_rl_jetpt_neqcharge, hden_eff_rl_jetpt_neqcharge, 1, 1, "B");

                hpurity_rl_jetpt_neqcharge->Draw("col");
                for (int i = 2; i < hpurity_rl_jetpt_neqcharge->GetNbinsX(); ++i) {
                        for (int j = 3; j <= hpurity_rl_jetpt_neqcharge->GetNbinsY(); ++j) {
                                double x = hpurity_rl_jetpt_neqcharge->GetXaxis()->GetBinCenter(i);
                                double y = hpurity_rl_jetpt_neqcharge->GetYaxis()->GetBinCenter(j);
                                double content = hpurity_rl_jetpt_neqcharge->GetBinContent(i, j);
                                double error = hpurity_rl_jetpt_neqcharge->GetBinError(i, j);

                                latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                        }
                }

                hpurity_rl_jetpt_neqcharge->SetTitle("Purity Correction Eq. Charge;R_{L};p_{T,jet}(GeV)");
                hpurity_rl_jetpt_neqcharge->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                hpurity_rl_jetpt_neqcharge->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                c->Print("../src-analysis/plots/pair_purity_correction_rl_jetpt_neqcharge.pdf");

                hefficiency_rl_jetpt_neqcharge->Draw("col");
                for (int i = 2; i < hefficiency_rl_jetpt_neqcharge->GetNbinsX(); ++i) {
                        for (int j = 3; j <= hefficiency_rl_jetpt_neqcharge->GetNbinsY(); ++j) {
                                double x = hefficiency_rl_jetpt_neqcharge->GetXaxis()->GetBinCenter(i);
                                double y = hefficiency_rl_jetpt_neqcharge->GetYaxis()->GetBinCenter(j);
                                double content = hefficiency_rl_jetpt_neqcharge->GetBinContent(i, j);
                                double error = hefficiency_rl_jetpt_neqcharge->GetBinError(i, j);

                                latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                        }
                }
                hefficiency_rl_jetpt_neqcharge->SetTitle("Efficiency Correction Eq. Charge;R_{L};p_{T,jet}(GeV)");
                hefficiency_rl_jetpt_neqcharge->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
                hefficiency_rl_jetpt_neqcharge->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                c->Print("../src-analysis/plots/pair_efficiency_correction_rl_jetpt_neqcharge.pdf");
        }

        // Create necessary 4vectors
        TLorentzVector* Jet_4vector = new TLorentzVector();
        TLorentzVector* Z0_4vector  = new TLorentzVector();
        TLorentzVector* mum_4vector = new TLorentzVector();
        TLorentzVector* mup_4vector = new TLorentzVector();
        TLorentzVector* h1_4vector  = new TLorentzVector();
        TLorentzVector* h2_4vector  = new TLorentzVector();
        
        TH3F* h_npair       = new TH3F("h_npair"      ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_npair_wmuon = new TH3F("h_npair_wmuon","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        h_npair->Sumw2();
        
        TH3F* h_eqchnpair        = new TH3F("h_eqchnpair"       ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_eqchnpair_wmuon  = new TH3F("h_eqchnpair_wmuon" ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_neqchnpair       = new TH3F("h_neqchnpair"      ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        TH3F* h_neqchnpair_wmuon = new TH3F("h_neqchnpair_wmuon","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning, nbin_jet_pt_unfolding, unfolding_jet_pt_binning, nbin_weight, weight_binning);
        h_eqchnpair->Sumw2();
        h_neqchnpair->Sumw2();
        
        TH1F* h_njet          = new TH1F("h_njet"         ,"",nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        TH1F* h_njet_wmuoneff = new TH1F("h_njet_wmuoneff","",nbin_jet_pt_unfolding, unfolding_jet_pt_binning);
        h_njet->Sumw2();

        TRandom3* rndm = new TRandom3(0);

        unsigned long long last_eventNum = 0;

        double printcounter = 0;
        
        std::cout<<"Working with 2016 data."<<std::endl;
        for (int evt = 0 ; evt < datatree_2016->fChain->GetEntries() ; evt++) {
                datatree_2016->GetEntry(evt);

                if (evt%10000 == 0) {
                        double percentage = 100.*evt/datatree_2016->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (evt != 0)
                        if (last_eventNum == datatree_2016->eventNumber) 
                                continue;

                if (datatree_2016->nPV != 1) 
                        continue;

                bool mum_trigger = (datatree_2016->mum_L0MuonEWDecision_TOS == 1 && 
                                    datatree_2016->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    datatree_2016->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

                bool mup_trigger = (datatree_2016->mup_L0MuonEWDecision_TOS == 1 && 
                                    datatree_2016->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    datatree_2016->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

                if (!mum_trigger && !mup_trigger) 
                        continue;
                
                Jet_4vector->SetPxPyPzE(datatree_2016->Jet_PX/1000.,
                                        datatree_2016->Jet_PY/1000.,
                                        datatree_2016->Jet_PZ/1000.,
                                        datatree_2016->Jet_PE/1000.);

                if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) 
                        continue;

                mum_4vector->SetPxPyPzE(datatree_2016->mum_PX/1000.,
                                        datatree_2016->mum_PY/1000.,
                                        datatree_2016->mum_PZ/1000.,
                                        datatree_2016->mum_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) 
                        continue;
                
                mup_4vector->SetPxPyPzE(datatree_2016->mup_PX/1000.,
                                        datatree_2016->mup_PY/1000.,
                                        datatree_2016->mup_PZ/1000.,
                                        datatree_2016->mup_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) 
                        continue;
                
                Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),
                                       mup_4vector->Py()+mum_4vector->Py(),
                                       mup_4vector->Pz()+mum_4vector->Pz(),
                                       mup_4vector->E() +mum_4vector->E());

                if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M()))
                        continue;

                double mup_pt  = (mup_4vector->Pt() > 70.) ? 69. : mup_4vector->Pt();
                double mum_pt  = (mum_4vector->Pt() > 70.) ? 69. : mum_4vector->Pt();
                double mup_eta = mup_4vector->Eta();
                double mum_eta = mum_4vector->Eta();

                double mup_eff_id  = h2_muon_2016_ideff_data->GetBinContent(h2_muon_2016_ideff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trk = h2_muon_2016_trkeff_data->GetBinContent(h2_muon_2016_trkeff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trg = h2_muon_2016_trgeff_data->GetBinContent(h2_muon_2016_trgeff_data->FindBin(mup_eta, mup_pt));

                double mum_eff_id  = h2_muon_2016_ideff_data->GetBinContent(h2_muon_2016_ideff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trk = h2_muon_2016_trkeff_data->GetBinContent(h2_muon_2016_trkeff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trg = h2_muon_2016_trgeff_data->GetBinContent(h2_muon_2016_trgeff_data->FindBin(mum_eta, mum_pt));

                double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(Jet_4vector->Pt()));
                
                double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(Jet_4vector->Pt()));

                if (get_muon) {
                        double mup_eff_id_err  = h2_muon_2016_ideff_data->GetBinError(h2_muon_2016_ideff_data->FindBin(mup_eta, mup_pt));
                        double mup_eff_trk_err = h2_muon_2016_trkeff_data->GetBinError(h2_muon_2016_trkeff_data->FindBin(mup_eta, mup_pt));
                        double mup_eff_trg_err = h2_muon_2016_trgeff_data->GetBinError(h2_muon_2016_trgeff_data->FindBin(mup_eta, mup_pt));
                        double mum_eff_id_err  = h2_muon_2016_ideff_data->GetBinError(h2_muon_2016_ideff_data->FindBin(mum_eta, mum_pt));
                        double mum_eff_trk_err = h2_muon_2016_trkeff_data->GetBinError(h2_muon_2016_trkeff_data->FindBin(mum_eta, mum_pt));
                        double mum_eff_trg_err = h2_muon_2016_trgeff_data->GetBinError(h2_muon_2016_trgeff_data->FindBin(mum_eta, mum_pt));

                        if (rndm->Integer(2)) {
                                mup_eff_id  += mup_eff_id_err;
                                mup_eff_trk  += mup_eff_trk_err;
                                mup_eff_trg  += mup_eff_trg_err;
                                mum_eff_id  += mum_eff_id_err;
                                mum_eff_trk  += mum_eff_trk_err;
                                mum_eff_trg  += mum_eff_trg_err;
                        } else {
                                mup_eff_id  -= mup_eff_id_err;
                                mup_eff_trk  -= mup_eff_trk_err;
                                mup_eff_trg  -= mup_eff_trg_err;
                                mum_eff_id  -= mum_eff_id_err;
                                mum_eff_trk  -= mum_eff_trk_err;
                                mum_eff_trg  -= mum_eff_trg_err;
                        }
                }
                
                double muon_weight = 1./(mum_eff_id*mup_eff_id*mum_eff_trk*mup_eff_trk*(mum_eff_trg+mup_eff_trg-mum_eff_trg*mup_eff_trg));

                for (int h1_index = 0 ; h1_index < datatree_2016->Jet_NDtr ; h1_index++) {
                        if (datatree_2016->Jet_Dtr_IsMeson[h1_index] != 1 && datatree_2016->Jet_Dtr_IsBaryon[h1_index] != 1)
                                continue;

                        h1_4vector->SetPxPyPzE(datatree_2016->Jet_Dtr_PX[h1_index]/1000.,
                                               datatree_2016->Jet_Dtr_PY[h1_index]/1000.,
                                               datatree_2016->Jet_Dtr_PZ[h1_index]/1000.,
                                               datatree_2016->Jet_Dtr_E[h1_index]/1000.);

                        if (!apply_chargedtrack_cuts(datatree_2016->Jet_Dtr_ThreeCharge[h1_index],
                                                     h1_4vector->P(),
                                                     h1_4vector->Pt(),
                                                     datatree_2016->Jet_Dtr_TrackChi2[h1_index]/datatree_2016->Jet_Dtr_TrackNDF[h1_index],
                                                     datatree_2016->Jet_Dtr_ProbNNghost[h1_index],
                                                     h1_4vector->Eta())) 
                                continue;
                
                        for (int h2_index = h1_index+1 ; h2_index < datatree_2016->Jet_NDtr ; h2_index++) {
                                if (datatree_2016->Jet_Dtr_IsMeson[h2_index] != 1 && datatree_2016->Jet_Dtr_IsBaryon[h2_index] != 1) 
                                        continue;

                                h2_4vector->SetPxPyPzE(datatree_2016->Jet_Dtr_PX[h2_index]/1000.,
                                                       datatree_2016->Jet_Dtr_PY[h2_index]/1000.,
                                                       datatree_2016->Jet_Dtr_PZ[h2_index]/1000.,
                                                       datatree_2016->Jet_Dtr_E[h2_index]/1000.);

                                if (!apply_chargedtrack_cuts(datatree_2016->Jet_Dtr_ThreeCharge[h2_index],
                                                             h2_4vector->P(),
                                                             h2_4vector->Pt(),
                                                             datatree_2016->Jet_Dtr_TrackChi2[h2_index]/datatree_2016->Jet_Dtr_TrackNDF[h2_index],
                                                             datatree_2016->Jet_Dtr_ProbNNghost[h2_index],
                                                             h2_4vector->Eta())) 
                                        continue;

                                double momentum_weight = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());

                                h_npair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                h_npair_wmuon->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight, muon_weight);

                                if (datatree_2016->Jet_Dtr_ThreeCharge[h1_index] * datatree_2016->Jet_Dtr_ThreeCharge[h2_index] > 0) {
                                        h_eqchnpair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                        h_eqchnpair_wmuon->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight, muon_weight);
                                }
                                else if (datatree_2016->Jet_Dtr_ThreeCharge[h1_index] * datatree_2016->Jet_Dtr_ThreeCharge[h2_index] < 0) {
                                        h_neqchnpair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                        h_neqchnpair_wmuon->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight, muon_weight);
                                }
                        }
                }
                h_njet->Fill(Jet_4vector->Pt());
                h_njet_wmuoneff->Fill(Jet_4vector->Pt(), muon_weight);

                last_eventNum = datatree_2016->eventNumber;
        }

        last_eventNum = 0;

        std::cout<<std::endl;
        std::cout<<"Working with 2017 data."<<std::endl;
        for (int evt = 0 ; evt < datatree_2017->fChain->GetEntries() ; evt++) {
                datatree_2017->GetEntry(evt);

                if (evt%10000 == 0) {
                        double percentage = 100.*evt/datatree_2017->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (evt != 0)
                        if (last_eventNum == datatree_2017->eventNumber) 
                                continue;

                if (datatree_2017->nPV != 1)
                        continue;

                bool mum_trigger = (datatree_2017->mum_L0MuonEWDecision_TOS == 1 && 
                                    datatree_2017->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    datatree_2017->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

                bool mup_trigger = (datatree_2017->mup_L0MuonEWDecision_TOS == 1 && 
                                    datatree_2017->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    datatree_2017->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);                                    

                if (!mum_trigger && !mup_trigger)
                        continue;
                
                Jet_4vector->SetPxPyPzE(datatree_2017->Jet_PX/1000.,
                                        datatree_2017->Jet_PY/1000.,
                                        datatree_2017->Jet_PZ/1000.,
                                        datatree_2017->Jet_PE/1000.);
                
                if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) 
                        continue;
                
                mum_4vector->SetPxPyPzE(datatree_2017->mum_PX/1000.,
                                        datatree_2017->mum_PY/1000.,
                                        datatree_2017->mum_PZ/1000.,
                                        datatree_2017->mum_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) 
                        continue;
                
                mup_4vector->SetPxPyPzE(datatree_2017->mup_PX/1000.,
                                        datatree_2017->mup_PY/1000.,
                                        datatree_2017->mup_PZ/1000.,
                                        datatree_2017->mup_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) 
                        continue;
                
                Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),
                                       mup_4vector->Py()+mum_4vector->Py(),
                                       mup_4vector->Pz()+mum_4vector->Pz(),
                                       mup_4vector->E() +mum_4vector->E());

                if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) 
                        continue;

                double mup_pt  = (mup_4vector->Pt() > 70.) ? 69. : mup_4vector->Pt();
                double mum_pt  = (mum_4vector->Pt() > 70.) ? 69. : mum_4vector->Pt();
                double mup_eta = mup_4vector->Eta();
                double mum_eta = mum_4vector->Eta();

                double mup_eff_id  = h2_muon_2017_ideff_data->GetBinContent(h2_muon_2017_ideff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trk = h2_muon_2017_trkeff_data->GetBinContent(h2_muon_2017_trkeff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trg = h2_muon_2017_trgeff_data->GetBinContent(h2_muon_2017_trgeff_data->FindBin(mup_eta, mup_pt));

                double mum_eff_id  = h2_muon_2017_ideff_data->GetBinContent(h2_muon_2017_ideff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trk = h2_muon_2017_trkeff_data->GetBinContent(h2_muon_2017_trkeff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trg = h2_muon_2017_trgeff_data->GetBinContent(h2_muon_2017_trgeff_data->FindBin(mum_eta, mum_pt));

                double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(Jet_4vector->Pt()));
                
                double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(Jet_4vector->Pt()));

                if (get_muon) {
                        double mup_eff_id_err  = h2_muon_2017_ideff_data->GetBinError(h2_muon_2017_ideff_data->FindBin(mup_eta, mup_pt));
                        double mup_eff_trk_err = h2_muon_2017_trkeff_data->GetBinError(h2_muon_2017_trkeff_data->FindBin(mup_eta, mup_pt));
                        double mup_eff_trg_err = h2_muon_2017_trgeff_data->GetBinError(h2_muon_2017_trgeff_data->FindBin(mup_eta, mup_pt));
                        double mum_eff_id_err  = h2_muon_2017_ideff_data->GetBinError(h2_muon_2017_ideff_data->FindBin(mum_eta, mum_pt));
                        double mum_eff_trk_err = h2_muon_2017_trkeff_data->GetBinError(h2_muon_2017_trkeff_data->FindBin(mum_eta, mum_pt));
                        double mum_eff_trg_err = h2_muon_2017_trgeff_data->GetBinError(h2_muon_2017_trgeff_data->FindBin(mum_eta, mum_pt));

                        if (rndm->Integer(2)) {
                                mup_eff_id  += mup_eff_id_err;
                                mup_eff_trk  += mup_eff_trk_err;
                                mup_eff_trg  += mup_eff_trg_err;
                                mum_eff_id  += mum_eff_id_err;
                                mum_eff_trk  += mum_eff_trk_err;
                                mum_eff_trg  += mum_eff_trg_err;
                        } else {
                                mup_eff_id  -= mup_eff_id_err;
                                mup_eff_trk  -= mup_eff_trk_err;
                                mup_eff_trg  -= mup_eff_trg_err;
                                mum_eff_id  -= mum_eff_id_err;
                                mum_eff_trk  -= mum_eff_trk_err;
                                mum_eff_trg  -= mum_eff_trg_err;
                        }
                }
                
                double muon_weight = 1./(mum_eff_id*mup_eff_id*mum_eff_trk*mup_eff_trk*(mum_eff_trg+mup_eff_trg-mum_eff_trg*mup_eff_trg));

                for (int h1_index = 0 ; h1_index < datatree_2017->Jet_NDtr ; h1_index++) {
                        if (datatree_2017->Jet_Dtr_IsMeson[h1_index] != 1 && datatree_2017->Jet_Dtr_IsBaryon[h1_index] != 1) 
                                continue;

                        h1_4vector->SetPxPyPzE(datatree_2017->Jet_Dtr_PX[h1_index]/1000.,
                                               datatree_2017->Jet_Dtr_PY[h1_index]/1000.,
                                               datatree_2017->Jet_Dtr_PZ[h1_index]/1000.,
                                               datatree_2017->Jet_Dtr_E[h1_index]/1000.);

                        if (!apply_chargedtrack_cuts(datatree_2017->Jet_Dtr_ThreeCharge[h1_index],
                                                     h1_4vector->P(),
                                                     h1_4vector->Pt(),
                                                     datatree_2017->Jet_Dtr_TrackChi2[h1_index]/datatree_2017->Jet_Dtr_TrackNDF[h1_index],
                                                     datatree_2017->Jet_Dtr_ProbNNghost[h1_index],
                                                     h1_4vector->Eta())) 
                                continue;
                        
                        for (int h2_index = h1_index+1 ; h2_index < datatree_2017->Jet_NDtr ; h2_index++) {
                                if (datatree_2017->Jet_Dtr_IsMeson[h2_index] != 1 && datatree_2017->Jet_Dtr_IsBaryon[h2_index] != 1) 
                                        continue;

                                h2_4vector->SetPxPyPzE(datatree_2017->Jet_Dtr_PX[h2_index]/1000.,
                                                       datatree_2017->Jet_Dtr_PY[h2_index]/1000.,
                                                       datatree_2017->Jet_Dtr_PZ[h2_index]/1000.,
                                                       datatree_2017->Jet_Dtr_E[h2_index]/1000.);

                                if (!apply_chargedtrack_cuts(datatree_2017->Jet_Dtr_ThreeCharge[h2_index],
                                                             h2_4vector->P(),
                                                             h2_4vector->Pt(),
                                                             datatree_2017->Jet_Dtr_TrackChi2[h2_index]/datatree_2017->Jet_Dtr_TrackNDF[h2_index],
                                                             datatree_2017->Jet_Dtr_ProbNNghost[h2_index],
                                                             h2_4vector->Eta())) 
                                        continue;

                                double momentum_weight = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());

                                h_npair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                h_npair_wmuon->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight, muon_weight);

                                if (datatree_2017->Jet_Dtr_ThreeCharge[h1_index] * datatree_2017->Jet_Dtr_ThreeCharge[h2_index] > 0) {
                                        h_eqchnpair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                        h_eqchnpair_wmuon->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight, muon_weight);
                                }
                                else if (datatree_2017->Jet_Dtr_ThreeCharge[h1_index] * datatree_2017->Jet_Dtr_ThreeCharge[h2_index] < 0) {
                                        h_neqchnpair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                        h_neqchnpair_wmuon->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight, muon_weight);
                                }
                        }
                }

                h_njet->Fill(Jet_4vector->Pt());
                h_njet_wmuoneff->Fill(Jet_4vector->Pt(), muon_weight);
                
                last_eventNum = datatree_2017->eventNumber;
        }

        last_eventNum = 0;

        std::cout<<std::endl;
        std::cout<<"Working with 2018 data."<<std::endl;
        for (int evt = 0 ; evt < datatree_2018->fChain->GetEntries() ; evt++) {
                datatree_2018->GetEntry(evt);

                if (evt%10000 == 0) {
                        double percentage = 100.*evt/datatree_2018->fChain->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (evt != 0)
                        if (last_eventNum == datatree_2018->eventNumber) 
                                continue;

                if (datatree_2018->nPV != 1) 
                        continue;

                bool mum_trigger = (datatree_2018->mum_L0MuonEWDecision_TOS == 1 && 
                                    datatree_2018->mum_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    datatree_2018->mum_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

                bool mup_trigger = (datatree_2018->mup_L0MuonEWDecision_TOS == 1 && 
                                    datatree_2018->mup_Hlt1SingleMuonHighPTDecision_TOS == 1 && 
                                    datatree_2018->mup_Hlt2EWSingleMuonVHighPtDecision_TOS == 1);

                if (!mum_trigger && !mup_trigger) 
                        continue;
                
                Jet_4vector->SetPxPyPzE(datatree_2018->Jet_PX/1000.,
                                        datatree_2018->Jet_PY/1000.,
                                        datatree_2018->Jet_PZ/1000.,
                                        datatree_2018->Jet_PE/1000.);

                if (!apply_jet_cuts(Jet_4vector->Eta(), Jet_4vector->Pt())) 
                        continue;
                
                mum_4vector->SetPxPyPzE(datatree_2018->mum_PX/1000.,
                                        datatree_2018->mum_PY/1000.,
                                        datatree_2018->mum_PZ/1000.,
                                        datatree_2018->mum_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mum_4vector), mum_4vector->Pt(), mum_4vector->Eta())) 
                        continue;
                
                mup_4vector->SetPxPyPzE(datatree_2018->mup_PX/1000.,
                                        datatree_2018->mup_PY/1000.,
                                        datatree_2018->mup_PZ/1000.,
                                        datatree_2018->mup_PE/1000.);

                if (!apply_muon_cuts(Jet_4vector->DeltaR(*mup_4vector), mup_4vector->Pt(), mup_4vector->Eta())) 
                        continue;
                
                Z0_4vector->SetPxPyPzE(mup_4vector->Px()+mum_4vector->Px(),
                                       mup_4vector->Py()+mum_4vector->Py(),
                                       mup_4vector->Pz()+mum_4vector->Pz(),
                                       mup_4vector->E() +mum_4vector->E());

                if (!apply_zboson_cuts(TMath::Abs(Jet_4vector->DeltaPhi(*Z0_4vector)), Z0_4vector->M())) 
                        continue;

                double mup_pt  = (mup_4vector->Pt() > 70.) ? 69. : mup_4vector->Pt();
                double mum_pt  = (mum_4vector->Pt() > 70.) ? 69. : mum_4vector->Pt();
                double mup_eta = mup_4vector->Eta();
                double mum_eta = mum_4vector->Eta();

                double mup_eff_id  = h2_muon_2018_ideff_data->GetBinContent(h2_muon_2018_ideff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trk = h2_muon_2018_trkeff_data->GetBinContent(h2_muon_2018_trkeff_data->FindBin(mup_eta, mup_pt));
                double mup_eff_trg = h2_muon_2018_trgeff_data->GetBinContent(h2_muon_2018_trgeff_data->FindBin(mup_eta, mup_pt));

                double mum_eff_id  = h2_muon_2018_ideff_data->GetBinContent(h2_muon_2018_ideff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trk = h2_muon_2018_trkeff_data->GetBinContent(h2_muon_2018_trkeff_data->FindBin(mum_eta, mum_pt));
                double mum_eff_trg = h2_muon_2018_trgeff_data->GetBinContent(h2_muon_2018_trgeff_data->FindBin(mum_eta, mum_pt));

                double jet_efficiency = hefficiency_jet->GetBinContent(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity     = hpurity_jet->GetBinContent(hpurity_jet->FindBin(Jet_4vector->Pt()));
                
                double jet_efficiency_error = hefficiency_jet->GetBinError(hefficiency_jet->FindBin(Jet_4vector->Pt()));
                double jet_purity_error     = hpurity_jet->GetBinError(hpurity_jet->FindBin(Jet_4vector->Pt()));

                if (get_muon) {
                        double mup_eff_id_err  = h2_muon_2018_ideff_data->GetBinError(h2_muon_2018_ideff_data->FindBin(mup_eta, mup_pt));
                        double mup_eff_trk_err = h2_muon_2018_trkeff_data->GetBinError(h2_muon_2018_trkeff_data->FindBin(mup_eta, mup_pt));
                        double mup_eff_trg_err = h2_muon_2018_trgeff_data->GetBinError(h2_muon_2018_trgeff_data->FindBin(mup_eta, mup_pt));
                        double mum_eff_id_err  = h2_muon_2018_ideff_data->GetBinError(h2_muon_2018_ideff_data->FindBin(mum_eta, mum_pt));
                        double mum_eff_trk_err = h2_muon_2018_trkeff_data->GetBinError(h2_muon_2018_trkeff_data->FindBin(mum_eta, mum_pt));
                        double mum_eff_trg_err = h2_muon_2018_trgeff_data->GetBinError(h2_muon_2018_trgeff_data->FindBin(mum_eta, mum_pt));

                        if (rndm->Integer(2)) {
                                mup_eff_id  += mup_eff_id_err;
                                mup_eff_trk  += mup_eff_trk_err;
                                mup_eff_trg  += mup_eff_trg_err;
                                mum_eff_id  += mum_eff_id_err;
                                mum_eff_trk  += mum_eff_trk_err;
                                mum_eff_trg  += mum_eff_trg_err;
                        } else {
                                mup_eff_id  -= mup_eff_id_err;
                                mup_eff_trk  -= mup_eff_trk_err;
                                mup_eff_trg  -= mup_eff_trg_err;
                                mum_eff_id  -= mum_eff_id_err;
                                mum_eff_trk  -= mum_eff_trk_err;
                                mum_eff_trg  -= mum_eff_trg_err;
                        }
                }
                
                double muon_weight = 1./(mum_eff_id*mup_eff_id*mum_eff_trk*mup_eff_trk*(mum_eff_trg+mup_eff_trg-mum_eff_trg*mup_eff_trg));

                for (int h1_index = 0 ; h1_index < datatree_2018->Jet_NDtr ; h1_index++) {
                        if (datatree_2018->Jet_Dtr_IsMeson[h1_index] != 1 && datatree_2018->Jet_Dtr_IsBaryon[h1_index] != 1) 
                                continue;

                        h1_4vector->SetPxPyPzE(datatree_2018->Jet_Dtr_PX[h1_index]/1000.,
                                               datatree_2018->Jet_Dtr_PY[h1_index]/1000.,
                                               datatree_2018->Jet_Dtr_PZ[h1_index]/1000.,
                                               datatree_2018->Jet_Dtr_E[h1_index]/1000.);

                        if (!apply_chargedtrack_cuts(datatree_2018->Jet_Dtr_ThreeCharge[h1_index],
                                                     h1_4vector->P(),
                                                     h1_4vector->Pt(),
                                                     datatree_2018->Jet_Dtr_TrackChi2[h1_index]/datatree_2018->Jet_Dtr_TrackNDF[h1_index],
                                                     datatree_2018->Jet_Dtr_ProbNNghost[h1_index],
                                                     h1_4vector->Eta()))
                                continue;

                        for (int h2_index = h1_index+1 ; h2_index < datatree_2018->Jet_NDtr ; h2_index++) {
                                if (datatree_2018->Jet_Dtr_IsMeson[h2_index] != 1 && datatree_2018->Jet_Dtr_IsBaryon[h2_index] != 1) 
                                        continue;

                                h2_4vector->SetPxPyPzE(datatree_2018->Jet_Dtr_PX[h2_index]/1000.,
                                                       datatree_2018->Jet_Dtr_PY[h2_index]/1000.,
                                                       datatree_2018->Jet_Dtr_PZ[h2_index]/1000.,
                                                       datatree_2018->Jet_Dtr_E[h2_index]/1000.);

                                if (!apply_chargedtrack_cuts(datatree_2018->Jet_Dtr_ThreeCharge[h2_index],
                                                             h2_4vector->P(),
                                                             h2_4vector->Pt(),
                                                             datatree_2018->Jet_Dtr_TrackChi2[h2_index]/datatree_2018->Jet_Dtr_TrackNDF[h2_index],
                                                             datatree_2018->Jet_Dtr_ProbNNghost[h2_index],
                                                             h2_4vector->Eta()))
                                        continue;

                                double momentum_weight = weight(h1_4vector->Pt(), h2_4vector->Pt(), Jet_4vector->Pt());

                                h_npair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                h_npair_wmuon->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight, muon_weight);

                                if (datatree_2018->Jet_Dtr_ThreeCharge[h1_index] * datatree_2018->Jet_Dtr_ThreeCharge[h2_index] > 0) {
                                        h_eqchnpair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                        h_eqchnpair_wmuon->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight, muon_weight);
                                }
                                else if (datatree_2018->Jet_Dtr_ThreeCharge[h1_index] * datatree_2018->Jet_Dtr_ThreeCharge[h2_index] < 0) {
                                        h_neqchnpair->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight);
                                        h_neqchnpair_wmuon->Fill(h1_4vector->DeltaR(*h2_4vector), Jet_4vector->Pt(), momentum_weight, muon_weight);
                                }
                        }
                }

                h_njet->Fill(Jet_4vector->Pt());
                h_njet_wmuoneff->Fill(Jet_4vector->Pt(), muon_weight);
                
                last_eventNum = datatree_2018->eventNumber;
        }

        fout->cd();
        hpurity_jet->Write();
        hefficiency_jet->Write();
        hpurity->Write();
        hefficiency->Write();
        h_njet->Write();
        h_njet_wmuoneff->Write();
        h_npair->Write();
        h_npair_wmuon->Write();
        h_eqchnpair->Write();
        h_eqchnpair_wmuon->Write();
        h_neqchnpair->Write();
        h_neqchnpair_wmuon->Write();
        hpurity_eqcharge->Write();
        hefficiency_eqcharge->Write();
        hpurity_neqcharge->Write();
        hefficiency_neqcharge->Write();
        fout->Close();
        
        std::cout<<std::endl;

        return 0;
}
