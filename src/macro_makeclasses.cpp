#include <iostream>
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils.cpp"
#include "../include/utils.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

void macro_makeclasses()
{
        TChain* sim_mc     = new TChain("mcjettuple/MCJetTree");  
        TChain* sim_mcreco = new TChain("StdHltZJets/DecayTree");
        TChain* data_2016  = new TChain("StdHltZJets/DecayTree");
        TChain* data_2017  = new TChain("StdHltZJets/DecayTree");
        TChain* data_2018  = new TChain("StdHltZJets/DecayTree");

        TChain* data_all  = new TChain("StdHltZJets/DecayTree");
        
        TChain* sim_mcreco_ct = new TChain("StdHltZJets/DecayTree");
        TChain* pseudodata_ct = new TChain("StdHltZJets/DecayTree");
        TChain* sim_mc_ct     = new TChain("mcjettuple/MCJetTree");
        TChain* truth_ct      = new TChain("mcjettuple/MCJetTree");

        TChain* sim_28r1 = new TChain("StdHltZJets/DecayTree");
        TChain* sim_28r2 = new TChain("StdHltZJets/DecayTree");
        
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MD_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MU_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MD_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MU_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MD_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MU_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MD_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MU_2016_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2017_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2017_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2018_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2018_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MD_2018_02222025.root").c_str());
        sim_mc->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MU_2018_02222025.root").c_str());
        sim_mc->MakeClass("TZJetsMC");
        
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MD_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MU_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MD_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MU_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MD_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MU_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MD_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MU_2016_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2017_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2017_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2018_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2018_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MD_2018_02222025.root").c_str());
        sim_mcreco->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MU_2018_02222025.root").c_str());
        sim_mcreco->MakeClass("TZJetsMCReco");

        data_2016->Add((input_folder + "Zjet_Data_2016_MU_04212025.root").c_str());
        data_2016->Add((input_folder + "Zjet_Data_2016_MD_04212025.root").c_str());
        data_2017->Add((input_folder + "Zjet_Data_2017_MU_04222025.root").c_str());
        data_2017->Add((input_folder + "Zjet_Data_2017_MD_04222025.root").c_str());
        data_2018->Add((input_folder + "Zjet_Data_2018_MU_04232025.root").c_str());
        data_2018->Add((input_folder + "Zjet_Data_2018_MD_04242025.root").c_str());
        data_2016->MakeClass("TZJets2016Data");
        data_2017->MakeClass("TZJets2017Data");
        data_2018->MakeClass("TZJets2018Data");
        
        data_all->Add((input_folder + "Zjet_Data_2016_MU_04212025.root").c_str());
        data_all->Add((input_folder + "Zjet_Data_2016_MD_04212025.root").c_str());
        data_all->Add((input_folder + "Zjet_Data_2017_MU_04222025.root").c_str());
        data_all->Add((input_folder + "Zjet_Data_2017_MD_04222025.root").c_str());
        data_all->Add((input_folder + "Zjet_Data_2018_MU_04232025.root").c_str());
        data_all->Add((input_folder + "Zjet_Data_2018_MD_04242025.root").c_str());
        data_all->MakeClass("TZJetsData");
        
        // CT Correction Files
        sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MD_2016_02222025.root").c_str());
        // sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MU_2016_02222025.root").c_str());
        sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MD_2016_02222025.root").c_str());
        // sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MU_2016_02222025.root").c_str());
        sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MD_2016_02222025.root").c_str());
        // sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MU_2016_02222025.root").c_str());
        sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2016_02222025.root").c_str());
        // sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2016_02222025.root").c_str());
        sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MD_2016_02222025.root").c_str());
        // sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MU_2016_02222025.root").c_str());
        sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2017_02222025.root").c_str());
        // sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2017_02222025.root").c_str());
        sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2018_02222025.root").c_str());
        // sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2018_02222025.root").c_str());
        sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MD_2018_02222025.root").c_str());
        // sim_mc_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MU_2018_02222025.root").c_str());
        sim_mc_ct->MakeClass("TZJetsMCCTCorr");

        sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MD_2016_02222025.root").c_str());
        // sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MU_2016_02222025.root").c_str());
        sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MD_2016_02222025.root").c_str());
        // sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MU_2016_02222025.root").c_str());
        sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MD_2016_02222025.root").c_str());
        // sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MU_2016_02222025.root").c_str());
        sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2016_02222025.root").c_str());
        // sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2016_02222025.root").c_str());
        sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MD_2016_02222025.root").c_str());
        // sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MU_2016_02222025.root").c_str());
        sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2017_02222025.root").c_str());
        // sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2017_02222025.root").c_str());
        sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2018_02222025.root").c_str());
        // sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2018_02222025.root").c_str());
        sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MD_2018_02222025.root").c_str());
        // sim_mcreco_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MU_2018_02222025.root").c_str());
        sim_mcreco_ct->MakeClass("TZJetsMCRecoCTCorr");
        
        // PseudoData
        // pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MD_2016_02222025.root").c_str());
        pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MU_2016_02222025.root").c_str());
        // pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MD_2016_02222025.root").c_str());
        pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MU_2016_02222025.root").c_str());
        // pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MD_2016_02222025.root").c_str());
        pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MU_2016_02222025.root").c_str());
        // pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2016_02222025.root").c_str());
        pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2016_02222025.root").c_str());
        // pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MD_2016_02222025.root").c_str());
        pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MU_2016_02222025.root").c_str());
        // pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2017_02222025.root").c_str());
        pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2017_02222025.root").c_str());
        // pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2018_02222025.root").c_str());
        pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2018_02222025.root").c_str());
        // pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MD_2018_02222025.root").c_str());
        pseudodata_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MU_2018_02222025.root").c_str());
        pseudodata_ct->MakeClass("TZJetsPseudoData");
        
        //Truth 
        // truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MD_2016_02222025.root").c_str());
        truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09b_MU_2016_02222025.root").c_str());
        // truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MD_2016_02222025.root").c_str());
        truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MU_2016_02222025.root").c_str());
        // truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MD_2016_02222025.root").c_str());
        truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MU_2016_02222025.root").c_str());
        // truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2016_02222025.root").c_str());
        truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2016_02222025.root").c_str());
        // truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MD_2016_02222025.root").c_str());
        truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MU_2016_02222025.root").c_str());
        // truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2017_02222025.root").c_str());
        truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2017_02222025.root").c_str());
        // truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2018_02222025.root").c_str());
        truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2018_02222025.root").c_str());
        // truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MD_2018_02222025.root").c_str());
        truth_ct->Add((input_folder + "ic-files/Zhadron_MC_Sim10b_MU_2018_02222025.root").c_str());
        truth_ct->MakeClass("TZJetsTruth");

        sim_28r1->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MD_2016_02222025.root").c_str());
        sim_28r1->Add((input_folder + "ic-files/Zhadron_MC_Sim09c_MU_2016_02222025.root").c_str());
        sim_28r1->MakeClass("TZJetsMCReco28r1");

        sim_28r2->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MD_2016_02222025.root").c_str());
        sim_28r2->Add((input_folder + "ic-files/Zhadron_MC_Sim09j_MU_2016_02222025.root").c_str());
        sim_28r2->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MD_2016_02222025.root").c_str());
        sim_28r2->Add((input_folder + "ic-files/Zhadron_MC_Sim09l_MU_2016_02222025.root").c_str());
        sim_28r2->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MD_2016_02222025.root").c_str());
        sim_28r2->Add((input_folder + "ic-files/Zhadron_MC_Sim10a_MU_2016_02222025.root").c_str());
        sim_28r2->MakeClass("TZJetsMCReco28r2");
}