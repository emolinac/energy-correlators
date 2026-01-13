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
#include "../include/TZJetsData.h"
#include "../include/TZJetsData.C"
#include "../include/TZJetsMCReco.h"
#include "../include/TZJetsMCReco.C"

void macro_print_dimuon_mass_all()
{
        gStyle->SetOptStat(1110);
        
        TZJetsData*   datatree   = new TZJetsData();
        TZJetsMCReco* mcrecotree = new TZJetsMCReco();

        TH1F* h_data   = new TH1F("h_data"  ,"DiMuon Invariant Mass",100,55,130);
        TH1F* h_mcreco = new TH1F("h_mcreco","DiMuon Invariant Mass",100,55,130);
        set_histogram_style(h_data,   868, std_line_width, std_marker_style, std_marker_size);
        set_histogram_style(h_mcreco, 797, std_line_width, std_marker_style, std_marker_size);

        for (int evt = 0 ; evt < datatree->fChain->GetEntries() ; evt++) {
                datatree->GetEntry(evt);

                if(datatree->Z0_M < 0 || TMath::IsNaN(datatree->Z0_M))
                        continue;

                h_data->Fill(datatree->Z0_M/1000.);
        }

        for (int evt = 0 ; evt < mcrecotree->fChain->GetEntries() ; evt++) {
                mcrecotree->GetEntry(evt);

                if(mcrecotree->Z0_M < 0 || TMath::IsNaN(mcrecotree->Z0_M))
                        continue;

                h_mcreco->Fill(mcrecotree->Z0_M/1000.);
        }

        h_data->Scale(1./h_data->Integral());
        h_mcreco->Scale(1./h_mcreco->Integral());

        std::cout<<h_data->GetMean()<<std::endl;
        std::cout<<h_mcreco->GetMean()<<std::endl;

        THStack* hs = new THStack();
        hs->Add(h_data);
        hs->Add(h_mcreco);

        
        TCanvas* c = new TCanvas("c","",1920,1080);
        hs->Draw("NOSTACK");
        hs->SetTitle(";M_{#mu^{+}#mu^{-}} (GeV);");

        TLegend* l = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());
        l->AddEntry(h_data,"Data","lpf");
        l->AddEntry(h_mcreco,"MCReco","lpf");
        l->Draw("SAME");

        c->Print("./plots/dimuon_mass_all.pdf");
}               