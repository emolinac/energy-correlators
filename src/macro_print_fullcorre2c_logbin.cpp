#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorre2c_logbin(int niter = 20)
{
    // Open the necessary files
    TFile* fout = new TFile((output_folder+namef_histos_corr_e2c_logbin).c_str(),"RECREATE");
    gROOT->cd();

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
    TH3D* hpurcorr = new TH3D("hpurcorr","",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* hmeas    = new TH3D("hmeas"   ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* htrue    = new TH3D("htrue"   ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");

    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        // if(jet_pt_reco<unfolding_jetpt_binning[0]||jet_pt_reco>unfolding_jetpt_binning[4]) continue;
        if(R_L_truth==-999) continue;
    
        response->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
    }

    // Fill the purity corrected distributions
    TH2D* hunfolded_ratio   = new TH2D("hunfolded_ratio"  ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH3D* hpuritycorrected  = new TH3D("hpuritycorrected" ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* hpuritycorrected2 = new TH3D("hpuritycorrected2","",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    ntuple_data->Project("hpuritycorrected" ,"weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    ntuple_data->Project("hpuritycorrected2","weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    
    // Unfold the purity corrected pairs
    RooUnfoldBayes unfold(response, hpuritycorrected, niter);
    TH3D* hunfolded_bayes = (TH3D*) unfold.Hunfold();

    TH2D* hunfolded_bayes_rl_jetpt = (TH2D*) hunfolded_bayes->Project3D("yx");
    TH2D* hpuritycorrected2_rl_jetpt = (TH2D*) hpuritycorrected2->Project3D("yx");
    
    hunfolded_ratio->Divide(hunfolded_bayes_rl_jetpt,hpuritycorrected2_rl_jetpt,1,1);

    TH1F* hcorr_e2c[Nbin_jet_pt]; 
    TH1F* hcorr_npair[Nbin_jet_pt]; 
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
    hunfolded_ratio->GetXaxis()->SetRangeUser(R_L_min,R_L_max);
    hunfolded_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0],jet_pt_binning[3]);
    gPad->SetLogx(0);
    gPad->SetLogy(1);
    // c->Print(Form("./plots/unfolded3d_%initer_ratio_sepyears_linearbinning_unfbinvarv2.pdf",niter));

    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);

    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        hcorrref_jet[bin] = new TH1F(Form("hcorrref_jet[%i]" ,bin) ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]); 
        hcorr_e2c[bin]    = new TH1F(Form("hcorr_e2c[%i]",bin)     ,"",Nbin_R_L_logbin,rl_logbinning);
        hcorr_npair[bin]  = new TH1F(Form("hcorr_npair[%i]",bin)   ,"",Nbin_R_L_logbin,rl_logbinning);
        set_histogram_style(hcorr_e2c[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
    
        for(int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++)
        {
            ntuple_data->GetEntry(entry);

            if(jet_pt<jet_pt_binning[bin]||jet_pt>jet_pt_binning[bin+1]) continue;
            if(efficiency_relerror>corr_rel_error) continue;
            if(purity_relerror>corr_rel_error) continue;
            if(efficiency<=0||efficiency>1) continue;
            if(purity<=0||purity>1) continue;
            if((jet_pt>20&&jet_pt<30)&&weight_pt>0.1) continue;
            if((jet_pt>30&&jet_pt<50)&&weight_pt>0.07) continue;
            if((jet_pt>50&&jet_pt<100)&&weight_pt>0.04) continue;

            double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt,weight_pt));
            // double unfolding_weight = 1.;

            hcorr_e2c[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
            hcorr_npair[bin]->Fill(R_L,purity*unfolding_weight/efficiency);

            // hcorr_e2c[bin]->Fill(R_L,weight_pt*purity/efficiency);
            // hcorr_npair[bin]->Fill(R_L,purity/efficiency);
        }
        ntuple_jet->Project(Form("hcorrref_jet[%i]" ,bin),"jet_pt",jet_full_corr[bin]);
        hcorr_e2c[bin]->Scale(1./hcorrref_jet[bin]->Integral());
        hcorr_npair[bin]->Scale(1./hcorrref_jet[bin]->Integral());
        
        s_data->Add(hcorr_e2c[bin],"E");
        l_data->AddEntry(hcorr_e2c[bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]),"lf");

        fout->cd();
        hcorr_e2c[bin]->Write();
        hcorr_npair[bin]->Write();
        gROOT->cd();
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print(Form("./plots/corr_e2c_unf3d_%initer_sepyears_logbin_rlnbin%i.pdf",niter,Nbin_R_L_logbin));
}
