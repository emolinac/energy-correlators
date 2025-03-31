#include <iostream>
#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"

void macro_makeclasses()
{
    // Add the files located in the input folder
    TChain* chain1 = new TChain("mcjettuple/MCJetTree");  
    TChain* chain2 = new TChain("StdHltZJets/DecayTree");
    TChain* chain3 = new TChain("StdHltZJets/DecayTree");
    TChain* chain4 = new TChain("StdHltZJets/DecayTree");
    TChain* chain5 = new TChain("StdHltZJets/DecayTree");
    TChain* chain6 = new TChain("mcjettuple/MCJetTree");  

    chain1->Add((input_folder+"Zjet_MC_Sim09j_2016_MD_03112025_full.root").c_str());
    chain1->Add((input_folder+"Zjet_MC_Sim09j_2016_MU_03112025_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim09j_2016_MD_03112025_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim09j_2016_MU_03112025_full.root").c_str());
    chain1->Add((input_folder+"Zjet_MC_Sim09l_2016_MD_03112025_full.root").c_str());
    chain1->Add((input_folder+"Zjet_MC_Sim09l_2016_MU_03112025_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim09l_2016_MD_03112025_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim09l_2016_MU_03112025_full.root").c_str());

    chain3->Add((input_folder+"Zjet_Data_2016_MU_04062024.root").c_str());
    chain3->Add((input_folder+"Zjet_Data_2016_MD_04062024.root").c_str());
    chain4->Add((input_folder+"Zjet_Data_2017_MU_04062024.root").c_str());
    chain4->Add((input_folder+"Zjet_Data_2017_MD_04062024.root").c_str());
    chain5->Add((input_folder+"Zjet_Data_2018_MU_04062024.root").c_str());
    chain5->Add((input_folder+"Zjet_Data_2018_MD_04062024.root").c_str());

    chain6->Add((input_folder+"Zjet_MC_Sim09j_2016_MD_03252025_0_full.root").c_str());
    chain6->Add((input_folder+"Zjet_MC_Sim09j_2016_MU_03252025_0_full.root").c_str());
    chain6->Add((input_folder+"Zjet_MC_Sim09l_2016_MD_03272025.root").c_str());
    chain6->Add((input_folder+"Zjet_MC_Sim09l_2016_MU_03272025.root").c_str());
    chain6->Add((input_folder+"Zjet_MC_Sim10a_2016_MD_03272025.root").c_str());
    chain6->Add((input_folder+"Zjet_MC_Sim10a_2016_MU_03272025.root").c_str());
    
    chain1->MakeClass("TZJetsMC");
    chain2->MakeClass("TZJetsMCReco");
    chain3->MakeClass("TZJets2016Data");
    chain4->MakeClass("TZJets2017Data");
    chain5->MakeClass("TZJets2018Data");
    chain6->MakeClass("TZJetsMCAT");
}

