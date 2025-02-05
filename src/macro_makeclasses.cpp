#include <iostream>
#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"

void macro_makeclasses()
{
    // Add the files located in the input folder
    TChain* chain1 = new TChain("mcjettuple/MCJetTree");  
    TChain* chain2 = new TChain("StdHltZJets/DecayTree");
    TChain* chain3 = new TChain("StdHltZJets/DecayTree");    // Used for data

    // Matching analysis
    TChain* chain4 = new TChain("StdHltZJets/DecayTree");
    TChain* chain5 = new TChain("StdHltZJets/DecayTree");

    // IBRAHIMS FILES
    //chain1->Add((input_folder+"Zhadron_MCReco_Sim10a_MU_2016_08162023.root").c_str());
    //chain1->Add((input_folder+"Zhadron_MCReco_Sim10a_MD_2016_08162023.root").c_str());
    //chain1->Add((input_folder+"Zhadron_MCReco_Sim09j_MU_2016_08162023.root").c_str());
    //chain1->Add((input_folder+"Zhadron_MCReco_Sim09j_MD_2016_08162023.root").c_str());
    //chain1->Add((input_folder+"Zhadron_MCReco_Sim09l_MU_2016_08162023.root").c_str());
    //chain1->Add((input_folder+"Zhadron_MCReco_Sim09l_MD_2016_08162023.root").c_str());

    //chain2->Add((input_folder+"Zhadron_MCReco_Sim10a_MU_2016_08162023.root").c_str());
    //chain2->Add((input_folder+"Zhadron_MCReco_Sim10a_MD_2016_08162023.root").c_str());
    //chain2->Add((input_folder+"Zhadron_MCReco_Sim09j_MU_2016_08162023.root").c_str());
    //chain2->Add((input_folder+"Zhadron_MCReco_Sim09j_MD_2016_08162023.root").c_str());
    //chain2->Add((input_folder+"Zhadron_MCReco_Sim09l_MU_2016_08162023.root").c_str());
    //chain2->Add((input_folder+"Zhadron_MCReco_Sim09l_MD_2016_08162023.root").c_str());

//    chain1->Add((input_folder+"Zjet_MC_Sim09_2016_MD_10052024_full.root").c_str());
//    chain1->Add((input_folder+"Zjet_MC_Sim09_2016_MU_10052024_full.root").c_str());
//    chain2->Add((input_folder+"Zjet_MC_Sim09_2016_MD_10052024_full.root").c_str());
//    chain2->Add((input_folder+"Zjet_MC_Sim09_2016_MU_10052024_full.root").c_str());

    // My files with matched jet dtrs in matching procedure
    chain1->Add((input_folder+"Zjet_MC_Sim09_2016_MD_10122024_full.root").c_str());
    chain1->Add((input_folder+"Zjet_MC_Sim09_2016_MU_10122024_full.root").c_str());
    chain1->Add((input_folder+"Zjet_MC_Sim09l_2016_MD_10142024_full.root").c_str());
    chain1->Add((input_folder+"Zjet_MC_Sim09l_2016_MU_10142024_full.root").c_str());
    chain1->Add((input_folder+"Zjet_MC_Sim10_2016_MD_10142024_full.root").c_str());
    chain1->Add((input_folder+"Zjet_MC_Sim10_2016_MU_10142024_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim09_2016_MD_10122024_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim09_2016_MU_10122024_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim09l_2016_MD_10142024_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim09l_2016_MU_10142024_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim10_2016_MD_10142024_full.root").c_str());
    chain2->Add((input_folder+"Zjet_MC_Sim10_2016_MU_10142024_full.root").c_str());
//


    // My files without matched jet dtrs in matching procedure
    //chain1->Add((input_folder+"Zjet_MC_Sim09_2016_MU_12162024_full.root").c_str());
    //chain2->Add((input_folder+"Zjet_MC_Sim09_2016_MU_12162024_full.root").c_str());
    
    chain3->Add((input_folder+"Zjet_Data_2016_MU_04062024.root").c_str());
    chain3->Add((input_folder+"Zjet_Data_2016_MD_04062024.root").c_str());
    chain3->Add((input_folder+"Zjet_Data_2017_MU_04062024.root").c_str());
    chain3->Add((input_folder+"Zjet_Data_2017_MD_04062024.root").c_str());
    chain3->Add((input_folder+"Zjet_Data_2018_MU_04062024.root").c_str());
    chain3->Add((input_folder+"Zjet_Data_2018_MD_04062024.root").c_str());
    
    chain4->Add((input_folder+"match-test/Zjet_MC_Sim09_2016_MU_12132024_full.root").c_str());
    chain5->Add((input_folder+"match-test/Zjet_MC_Sim09_2016_MU_12132024_full_matchedjetinc.root").c_str());

    chain1->MakeClass("TZJetsMC");
    chain2->MakeClass("TZJetsMCReco");
    chain3->MakeClass("TZJetsData");
    chain4->MakeClass("TZJetsMCRecoNoMatchedJetDtrs");
    chain5->MakeClass("TZJetsMCRecoMatchedJetDtrs");
}

