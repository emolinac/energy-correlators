#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_closuretest_puritycorrunfolded3d_corre2c(int niter = 10)
{
    // Open the necessary files
    TFile* fcorr = new TFile((output_folder+namef_ntuple_e2c_corr).c_str()); 
    if(fcorr->IsZombie()) return;
    
    TNtuple* ntuple_data = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_jet  = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
    
    // Set the branches of data
    float R_L, jet_pt, weight_pt;
    float efficiency, purity, efficiency_relerror, purity_relerror;
    ntuple_data->SetBranchAddress("R_L",&R_L);
    ntuple_data->SetBranchAddress("jet_pt",&jet_pt);
    ntuple_data->SetBranchAddress("weight_pt",&weight_pt);
    ntuple_data->SetBranchAddress("efficiency",&efficiency);
    ntuple_data->SetBranchAddress("efficiency_relerror",&efficiency_relerror);
    ntuple_data->SetBranchAddress("purity",&purity);
    ntuple_data->SetBranchAddress("purity_relerror",&purity_relerror);

    // UNFOLDING FIRST
    TFile* f = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

    float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt_reco);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    ntuple->SetBranchAddress("R_L",&R_L_reco);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    ntuple->SetBranchAddress("weight_pt",&weight_pt_reco);
    ntuple->SetBranchAddress("weight_pt_truth",&weight_pt_truth);
    
    // Create histograms with the respective true and matched reco 
    TH3D* hpurcorr = new TH3D("hpurcorr","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* hmeas    = new TH3D("hmeas"   ,"",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* htrue    = new TH3D("htrue"   ,"",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");

    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        if(jet_pt_reco<unfolding_jetpt_binning[0]||jet_pt_reco>unfolding_jetpt_binning[4]) continue;
        if(R_L_truth==-999) continue;
    
        response->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
    }

    // Fill the purity corrected distributions
    TH2D* hunfolded_ratio   = new TH2D("hunfolded_ratio"  ,"",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    TH3D* hpuritycorrected  = new TH3D("hpuritycorrected" ,"",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* hpuritycorrected2 = new TH3D("hpuritycorrected2","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    ntuple_data->Project("hpuritycorrected" ,"weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    ntuple_data->Project("hpuritycorrected2","weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    
    // Unfold the purity corrected pairs
    RooUnfoldBayes unfold(response, hpuritycorrected, niter);
    TH3D* hunfolded_bayes = (TH3D*) unfold.Hunfold();

    TH2D* hunfolded_bayes_rl_jetpt = (TH2D*) hunfolded_bayes->Project3D("yx");
    TH2D* hpuritycorrected2_rl_jetpt = (TH2D*) hpuritycorrected2->Project3D("yx");
    
    hunfolded_ratio->Divide(hunfolded_bayes_rl_jetpt,hpuritycorrected2_rl_jetpt,1,1);

    TH1F* hcorr_data[Nbin_jet_pt]; 
    TH1F* hcorr_jet[Nbin_jet_pt]; 
    TH1F* hcorrref_jet[Nbin_jet_pt];
   
    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    gStyle->SetPaintTextFormat("4.2f");
    hunfolded_ratio->Draw("col text");
    hunfolded_ratio->SetTitle("Purity Corrected Unfolded/Purity Corrected;R_{L};p^{jet}_{T}GeV");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    c->Print(Form("./plots/unfolded3d_%initer_ratio.pdf",niter));

    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend();

    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        hcorrref_jet[bin] = new TH1F(Form("hcorrref_jet[%i]" ,bin) ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]); 
        hcorr_data[bin]   = new TH1F(Form("hcorr_data[%i]",bin),"",Nbin_R_L,rl_binning);
        set_histogram_style(hcorr_data[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
    
        for(int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++)
        {
            ntuple_data->GetEntry(entry);

            if(jet_pt<jet_pt_binning[bin]||jet_pt>jet_pt_binning[bin+1]) continue;
            if(efficiency_relerror>corr_rel_error) continue;
            if(purity_relerror>corr_rel_error) continue;
            if(efficiency<=0||efficiency>1) continue;
            if(purity<=0||purity>1) continue;

            double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt,weight_pt));
            hcorr_data[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
        }
        ntuple_jet->Project(Form("hcorrref_jet[%i]" ,bin),"jet_pt",jet_full_corr[bin]);
        hcorr_data[bin]->Scale(1./hcorrref_jet[bin]->Integral());
        
        s_data->Add(hcorr_data[bin],"E");
        l_data->AddEntry(hcorr_data[bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]),"lf");
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print(Form("./plots/corr_e2c_unf3d_%initer_finnerlogbinning_logscale.pdf",niter));
}
