#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_closuretest_jetpt(int Niter = 1, double jet_pt_min_local = unfolding_jetpt_binning[0], double jet_pt_max_local = unfolding_jetpt_binning[5])
{
    TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

    float jet_pt, jet_pt_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    
    // Create histograms with the respective true and matched reco 
    TH1F* hmeas = new TH1F("hmeas","",nbin_jet_pt+2,unfolding_jetpt_binning);
    TH1F* htrue = new TH1F("htrue","",nbin_jet_pt+2,unfolding_jetpt_binning);
    TH2F* hresp = new TH2F("hresp","",nbin_jet_pt+2,unfolding_jetpt_binning,nbin_jet_pt+2,unfolding_jetpt_binning);

    TH1F* htrue_ref = new TH1F("htrue_ref","",nbin_jet_pt+2,unfolding_jetpt_binning);
    TH1F* h_ct      = new TH1F("h_ct"     ,"",nbin_jet_pt+2,unfolding_jetpt_binning);

    TRandom3* rndm = new TRandom3();
    for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        // if (jet_pt<jet_pt_min_local||jet_pt>jet_pt_max_local) continue;
        if (jet_pt_truth==-999) continue;
        if (rndm->Uniform(1)<=0.5) 
        {
            htrue_ref->Fill(jet_pt_truth);
            continue;
        }

        hresp->Fill(jet_pt, jet_pt_truth);
        hmeas->Fill(jet_pt);
        htrue->Fill(jet_pt_truth);
    }

    // Create response matrix object
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, hresp);

    // // Draw response matrix
    // TH2F* hresponse = (TH2F*) response->Hresponse();
    // hresponse->Draw("col text");

    TCanvas* c = new TCanvas("c", "", 1920, 1080);
    c->Draw();

    RooUnfoldBayes unfold(response, hmeas, Niter);
    auto* hunfolded_bayes = unfold.Hunfold();

    htrue_ref->Scale(1./htrue_ref->Integral());
    hmeas->Scale(1./hmeas->Integral());
    htrue->Scale(1./htrue->Integral());
    hunfolded_bayes->Scale(1./hunfolded_bayes->Integral());

    h_ct->Divide(hunfolded_bayes,htrue_ref,1,1,"B");
    hmeas->Divide(htrue);

    set_histogram_style(h_ct, kViolet, std_line_width, std_marker_style, std_marker_size+1);
    set_histogram_style(hmeas, kCyan, std_line_width, std_marker_style, std_marker_size+1);

    h_ct->GetYaxis()->SetNdivisions(505);
    h_ct->Draw("E0");
    h_ct->SetTitle(";p^{jet}_{t}(GeV);Unfolded/True");
    // h_ct->GetYaxis()->SetRangeUser(0.9,1.1);

    // gPad->SetLogx(1);
    gPad->SetGridy(1);
    // c->Print(Form("./plots/closuretest_roounfold_jetpt_bayes_iter%i_jetptfrom%.0fto%.0f.pdf",Niter,jet_pt_min_local,jet_pt_max_local));    
}