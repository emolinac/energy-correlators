#include "../include/analysis-constants.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"
#include "../include/TZJetsUnfold.h"
#include "../include/TZJetsUnfold.C"


void macro_print_unfold_rl()
{
    TZJetsUnfold* tree = new TZJetsUnfold();

    gStyle->SetOptStat();

    // Create histograms with the respective true and matched reco 
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    TH1F* hmeas = new TH1F("hmeas","",Nbin_R_L,R_L_min,R_L_max);
    TH1F* htrue = new TH1F("htrue","",Nbin_R_L,R_L_min,R_L_max);
    TH2F* hresp = new TH2F("hresp","",Nbin_R_L,R_L_min,R_L_max,Nbin_R_L,R_L_min,R_L_max);
    hmeas->Sumw2();

    // Create response matrix object
    RooUnfoldResponse response(hmeas, htrue, hresp);

    for(int evt = 0 ; evt < tree->fChain->GetEntries() ; evt++)
    {
        // Access entry of tree
        tree->GetEntry(evt);

        // Jet cuts
        bool jet_cuts = (tree->jet_eta>2.5&&tree->jet_eta<4&&tree->jet_pt>jet_pt_min&&tree->jet_pt<jet_pt_max) ? 1 : 0;

        // Track cuts
        bool chi2ndf_cut  = (tree->h1_chi2/tree->h1_ndf<3&&tree->h2_chi2/tree->h2_ndf<3) ? 1 : 0;
        bool p_cut        = (tree->h1_p>4&&tree->h1_p<1000&&tree->h2_p>4&&tree->h2_p<1000) ? 1 : 0;
        bool pt_cut       = (tree->h1_pt>0.250&&tree->h2_pt>0.250) ? 1 : 0;
        bool pnnghost_cut = (tree->h1_probnnghost<0.5&&tree->h2_probnnghost<0.5) ? 1 : 0;

        // Topological cuts
        bool phi_zjet_cut     = (TMath::Abs(tree->deltaphi_z_jet)>7*TMath::Pi()/8.)? 1 : 0;
        bool phi_mumjet_cut   = (TMath::Abs(tree->deltaphi_mum_jet)>0.4)? 1 : 0;
        bool phi_mupjet_cut   = (TMath::Abs(tree->deltaphi_mup_jet)>0.4)? 1 : 0;

        // Muon cuts
        bool mu_pt_cut        = (tree->mum_pt>20.&&tree->mup_pt>20.)? 1 : 0;
        bool mu_eta_cut       = (tree->mum_eta>2&&tree->mum_eta<4.5&&tree->mup_eta>2&&tree->mup_eta<4.5)? 1 : 0;
        bool mum_mup_mass_cut = (sqrt(tree->mum_m*tree->mum_m + tree->mup_m*tree->mup_m + 2*(tree->mum_pe*tree->mup_pe - tree->mum_px*tree->mup_px - tree->mum_py*tree->mup_py - tree->mum_pz*tree->mup_pz))>60 &&
                                 sqrt(tree->mum_m*tree->mum_m + tree->mup_m*tree->mup_m + 2*(tree->mum_pe*tree->mup_pe - tree->mum_px*tree->mup_px - tree->mum_py*tree->mup_py - tree->mum_pz*tree->mup_pz))<120)? 1 : 0;
        bool mu_trackprob_cut = (tree->mum_probchi2>0.001&&tree->mup_probchi2>0.001)? 1 : 0;

        bool h1_noneutrals_cut = (tree->h1_charge!=0)? 1 : 0;
        bool h2_noneutrals_cut = (tree->h2_charge!=0)? 1 : 0;

        bool total_req = (jet_cuts&&chi2ndf_cut&&p_cut&&pt_cut&&pnnghost_cut&&phi_zjet_cut&&phi_mumjet_cut&&phi_mupjet_cut&&mu_pt_cut&&mu_eta_cut&&mum_mup_mass_cut&&mu_trackprob_cut&&h1_noneutrals_cut&&h2_noneutrals_cut);

        //if(total_req) hmeas->Fill(tree->R_L);

        if(total_req&&tree->R_L_truth!=-999) {response.Fill(tree->R_L, tree->R_L_truth); htrue->Fill(tree->R_L_truth);hmeas->Fill(tree->R_L);}
        else if(total_req&&tree->R_L_truth==-999) response.Fake(tree->R_L);
    }

    // Draw response matrix
    TH2F* hresponse = (TH2F*) response.HresponseNoOverflow();
    RooUnfoldBayes unfold_bayes(&response, hmeas, 4);
    auto* hreco_bayes = unfold_bayes.Hunfold();

    set_histogram_style(hmeas      , kViolet+2, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hreco_bayes, kRed+2, std_line_width, std_marker_style, std_marker_size);


    //THStack* s = new THStack();
    //s->Add(hreco_bayes);
    //s->Add(hmeas);
    //s->Draw("NOSTACK");

    //gPad->SetLogx(1);
    //gPad->SetLogy(1);

    TLegend* l = new TLegend();
    l->AddEntry(hreco_bayes,"Unfolded","lpf");
    l->AddEntry(htrue      ,"True"    ,"lpf");
    l->Draw("SAME");

    TRatioPlot* rp = new TRatioPlot(htrue, hreco_bayes);
    rp->Draw();
}