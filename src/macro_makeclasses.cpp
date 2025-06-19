#include <iostream>
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/directories.h"
#include "../include/names.h"

void macro_makeclasses()
{
    TChain* sim_mc     = new TChain("mcjettuple/MCJetTree");  
    TChain* sim_mcreco = new TChain("StdHltZJets/DecayTree");
    TChain* data_2016  = new TChain("StdHltZJets/DecayTree");
    TChain* data_2017  = new TChain("StdHltZJets/DecayTree");
    TChain* data_2018  = new TChain("StdHltZJets/DecayTree");
    
    TChain* sim_mcreco_ct = new TChain("StdHltZJets/DecayTree");
    TChain* pseudodata_ct = new TChain("StdHltZJets/DecayTree");
    TChain* sim_mc_ct     = new TChain("mcjettuple/MCJetTree");
    TChain* truth_ct      = new TChain("mcjettuple/MCJetTree");
    
    TChain* data_2016_new  = new TChain("StdHltZJets/DecayTree");
    
    sim_mc->Add((input_folder+"Zjet_MC_Sim09j_2016_MD_04112025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim09j_2016_MU_04112025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim09l_2016_MD_04112025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim09l_2016_MU_04112025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim10a_2016_MD_04142025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim10a_2016_MU_04142025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim09l_2017_MD_04142025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim09l_2017_MU_04142025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim09l_2018_MD_04142025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim09l_2018_MU_04142025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim10a_2018_MD_04142025.root").c_str());
    sim_mc->Add((input_folder+"Zjet_MC_Sim10a_2018_MU_04142025.root").c_str());
    sim_mc->MakeClass("TZJetsMC");
    
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim09j_2016_MD_04112025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim09j_2016_MU_04112025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim09l_2016_MD_04112025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim09l_2016_MU_04112025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim10a_2016_MD_04142025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim10a_2016_MU_04142025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim09l_2017_MD_04142025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim09l_2017_MU_04142025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim09l_2018_MD_04142025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim09l_2018_MU_04142025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim10a_2018_MD_04142025.root").c_str());
    sim_mcreco->Add((input_folder+"Zjet_MC_Sim10a_2018_MU_04142025.root").c_str());
    sim_mcreco->MakeClass("TZJetsMCReco");

    data_2016->Add((input_folder+"Zjet_Data_2016_MU_04062024.root").c_str());
    data_2016->Add((input_folder+"Zjet_Data_2016_MD_04062024.root").c_str());
    data_2017->Add((input_folder+"Zjet_Data_2017_MU_04062024.root").c_str());
    data_2017->Add((input_folder+"Zjet_Data_2017_MD_04062024.root").c_str());
    data_2018->Add((input_folder+"Zjet_Data_2018_MU_04062024.root").c_str());
    data_2018->Add((input_folder+"Zjet_Data_2018_MD_04062024.root").c_str());
    data_2016->MakeClass("TZJets2016Data");
    data_2017->MakeClass("TZJets2017Data");
    data_2018->MakeClass("TZJets2018Data");
    
    // CT Correction Files
    sim_mc_ct->Add((input_folder+"Zjet_MC_Sim09j_2016_MD_04112025.root").c_str());
    sim_mc_ct->Add((input_folder+"Zjet_MC_Sim09j_2016_MU_04112025.root").c_str());
    sim_mc_ct->Add((input_folder+"Zjet_MC_Sim09l_2016_MD_04112025.root").c_str());
    sim_mc_ct->Add((input_folder+"Zjet_MC_Sim09l_2016_MU_04112025.root").c_str());
    sim_mc_ct->Add((input_folder+"Zjet_MC_Sim10a_2016_MD_04142025.root").c_str());
    sim_mc_ct->Add((input_folder+"Zjet_MC_Sim10a_2016_MU_04142025.root").c_str());
    sim_mc_ct->MakeClass("TZJetsMCCTCorr");
    
    sim_mcreco_ct->Add((input_folder+"Zjet_MC_Sim09j_2016_MD_04112025.root").c_str());
    sim_mcreco_ct->Add((input_folder+"Zjet_MC_Sim09j_2016_MU_04112025.root").c_str());
    sim_mcreco_ct->Add((input_folder+"Zjet_MC_Sim09l_2016_MD_04112025.root").c_str());
    sim_mcreco_ct->Add((input_folder+"Zjet_MC_Sim09l_2016_MU_04112025.root").c_str());
    sim_mcreco_ct->Add((input_folder+"Zjet_MC_Sim10a_2016_MD_04142025.root").c_str());
    sim_mcreco_ct->Add((input_folder+"Zjet_MC_Sim10a_2016_MU_04142025.root").c_str());
    sim_mcreco_ct->MakeClass("TZJetsMCRecoCTCorr");
    
    // PseudoData
    pseudodata_ct->Add((input_folder+"Zjet_MC_Sim09l_2017_MD_04142025.root").c_str());
    pseudodata_ct->Add((input_folder+"Zjet_MC_Sim09l_2017_MU_04142025.root").c_str());
    pseudodata_ct->MakeClass("TZJetsPseudoData");
    
    //Truth 
    truth_ct->Add((input_folder+"Zjet_MC_Sim09l_2017_MD_04142025.root").c_str());
    truth_ct->Add((input_folder+"Zjet_MC_Sim09l_2017_MU_04142025.root").c_str());
    truth_ct->MakeClass("TZJetsTruth");
    
    // New format data for the JEC Syst
    data_2016_new->Add((input_folder+"Zjet_Data_2016_MU_04212025.root").c_str());
    data_2016_new->Add((input_folder+"Zjet_Data_2016_MD_04212025.root").c_str());
    data_2016_new->MakeClass("TZJets2016NewData");

}
