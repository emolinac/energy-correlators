#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorre2c(int niter = 20)
{
    // Open the necessary files
    TFile* fout        = new TFile((output_folder+namef_histos_corr_e2c_logbin).c_str(),"RECREATE");
    TFile* fout_linear = new TFile((output_folder+namef_histos_corr_e2c).c_str(),"RECREATE");
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
    
    // Create histograms with different types of binning
    TH3D* hpurcorr = new TH3D("hpurcorr","",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* hmeas    = new TH3D("hmeas"   ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* htrue    = new TH3D("htrue"   ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");
    
    TH3D* hpurcorr_l = new TH3D("hpurcorr_l","",Nbin_R_L_unfolding,unfolding_rl_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* hmeas_l    = new TH3D("hmeas_l"   ,"",Nbin_R_L_unfolding,unfolding_rl_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* htrue_l    = new TH3D("htrue_l"   ,"",Nbin_R_L_unfolding,unfolding_rl_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    RooUnfoldResponse* response_l = new RooUnfoldResponse(hmeas_l, htrue_l, "response_l");

    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        // if(jet_pt_reco<unfolding_jetpt_binning[0]||jet_pt_reco>unfolding_jetpt_binning[4]) continue;
        if(R_L_truth==-999) continue;
    
        response->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
        response_l->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
    }

    // Fill the purity corrected distributions
    TH2D* hunfolded_ratio     = new TH2D("hunfolded_ratio"  ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH3D* hpuritycorrected    = new TH3D("hpuritycorrected" ,"",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* hpuritycorrected2   = new TH3D("hpuritycorrected2","",Nbin_R_L_logbin_unfolding,unfolding_rl_logbinning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH2D* hunfolded_ratio_l   = new TH2D("hunfolded_ratio_l"  ,"",Nbin_R_L_unfolding,unfolding_rl_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning);
    TH3D* hpuritycorrected_l  = new TH3D("hpuritycorrected_l" ,"",Nbin_R_L_unfolding,unfolding_rl_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    TH3D* hpuritycorrected2_l = new TH3D("hpuritycorrected2_l","",Nbin_R_L_unfolding,unfolding_rl_binning,Nbin_jet_pt_unfolding,unfolding_jetpt_binning,Nbin_weight,weight_binning);
    
    ntuple_data->Project("hpuritycorrected" ,"weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    ntuple_data->Project("hpuritycorrected2","weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    ntuple_data->Project("hpuritycorrected_l" ,"weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    ntuple_data->Project("hpuritycorrected2_l","weight_pt:jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    
    // Unfold the purity corrected pairs
    RooUnfoldBayes unfold(response, hpuritycorrected, niter);
    TH3D* hunfolded_bayes            = (TH3D*) unfold.Hunfold();
    TH2D* hunfolded_bayes_rl_jetpt   = (TH2D*) hunfolded_bayes->Project3D("yx");
    TH2D* hpuritycorrected2_rl_jetpt = (TH2D*) hpuritycorrected2->Project3D("yx");
    hunfolded_ratio->Divide(hunfolded_bayes_rl_jetpt,hpuritycorrected2_rl_jetpt,1,1);

    RooUnfoldBayes unfold_l(response_l, hpuritycorrected_l, niter);
    TH3D* hunfolded_bayes_l            = (TH3D*) unfold_l.Hunfold();
    TH2D* hunfolded_bayes_rl_jetpt_l   = (TH2D*) hunfolded_bayes_l->Project3D("yx");
    TH2D* hpuritycorrected2_rl_jetpt_l = (TH2D*) hpuritycorrected2_l->Project3D("yx");
    hunfolded_ratio_l->Divide(hunfolded_bayes_rl_jetpt_l,hpuritycorrected2_rl_jetpt_l,1,1);

    TH1F* hcorr_jet[Nbin_jet_pt];
    TH1F* hcorr_e2c[Nbin_jet_pt]; 
    TH1F* hcorr_npair[Nbin_jet_pt]; 
    TH1F* hcorr_e2c_l[Nbin_jet_pt]; 
    TH1F* hcorr_npair_l[Nbin_jet_pt]; 
    
    TCanvas* c = new TCanvas("c","",1920,1080);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    // gStyle->SetPaintTextFormat("4.2f");
    // hunfolded_ratio->Draw("col text");
    // hunfolded_ratio->SetTitle("Purity Corrected Unfolded/Purity Corrected;R_{L};p^{jet}_{T}GeV");
    // hunfolded_ratio->GetXaxis()->SetRangeUser(R_L_min,R_L_max);
    // hunfolded_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0],jet_pt_binning[3]);
    // gPad->SetLogx(0);
    // gPad->SetLogy(1);
    // c->Print(Form("./plots/unfolded3d_initer%_ratio_sepyears_linearbinning_unfbinvarv2.pdf",niter));

    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);

    // Fill the histograms
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        hcorr_jet[bin]      = new TH1F(Form("hcorr_jet[%i]" ,bin)  ,"", 1,jet_pt_binning[bin],jet_pt_binning[bin+1]); 
        
        hcorr_e2c[bin]      = new TH1F(Form("hcorr_e2c[%i]",bin)    ,"", Nbin_R_L_logbin,rl_logbinning);
        hcorr_npair[bin]    = new TH1F(Form("hcorr_npair[%i]",bin)  ,"", Nbin_R_L_logbin,rl_logbinning);
        hcorr_e2c_l[bin]    = new TH1F(Form("hcorr_e2c_l[%i]",bin)  ,"", Nbin_R_L       ,rl_binning   );
        hcorr_npair_l[bin]  = new TH1F(Form("hcorr_npair_l[%i]",bin),"", Nbin_R_L       ,rl_binning   );
        set_histogram_style(hcorr_e2c[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        set_histogram_style(hcorr_e2c_l[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
    
        for(int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++)
        {
            ntuple_data->GetEntry(entry);

            if(jet_pt<jet_pt_binning[bin]||jet_pt>jet_pt_binning[bin+1]) continue;
            if(efficiency_relerror>corr_rel_error) continue;
            if(purity_relerror>corr_rel_error) continue;
            if(efficiency<=0||efficiency>1) efficiency = 1;//continue;
            if(purity<=0||purity>1) purity = 1;//continue;
            if(weight_pt>weight_pt_cut[bin]) continue;
            
            double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt,weight_pt));
            if(unfolding_weight<=0) unfolding_weight = 1;

            double unfolding_weight_l = hunfolded_ratio_l->GetBinContent(hunfolded_ratio_l->FindBin(R_L,jet_pt,weight_pt));
            if(unfolding_weight_l<=0) unfolding_weight_l = 1;

            hcorr_e2c[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
            hcorr_npair[bin]->Fill(R_L,purity*unfolding_weight/efficiency);
            hcorr_e2c_l[bin]->Fill(R_L,purity*unfolding_weight_l*weight_pt/efficiency);
            hcorr_npair_l[bin]->Fill(R_L,purity*unfolding_weight_l/efficiency);
        }

        ntuple_jet->Project(Form("hcorr_jet[%i]" ,bin),"jet_pt",jet_full_corr[bin]);
    }

    // Normalize the distributions
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        std::cout<<"Jet Pt Bin = "<<bin<<std::endl;
        // Log binning
        for(int i = 1 ; i <= hcorr_e2c[bin]->GetNbinsX() ; i++)
        {
            hcorr_e2c[bin]->SetBinContent(i, hcorr_e2c[bin]->GetBinContent(i)/(rl_logbinning[i]-rl_logbinning[i-1]));
            hcorr_npair[bin]->SetBinContent(i, hcorr_npair[bin]->GetBinContent(i)/(rl_logbinning[i]-rl_logbinning[i-1]));

            std::cout<<hcorr_e2c[bin]->GetBinContent(i)<<std::endl;
        }
        hcorr_e2c[bin]->Scale(1./hcorr_jet[bin]->Integral());
        hcorr_npair[bin]->Scale(1./hcorr_jet[bin]->Integral());
        
        // Linear binning
        hcorr_e2c_l[bin]->Scale(1./hcorr_jet[bin]->Integral());
        hcorr_npair_l[bin]->Scale(1./hcorr_jet[bin]->Integral());

        fout->cd();
        hcorr_e2c[bin]->Write();
        hcorr_npair[bin]->Write();
        fout_linear->cd();
        hcorr_e2c_l[bin]->Write();
        hcorr_npair_l[bin]->Write();
        gROOT->cd();
    }

    // Draw Log binning distributions
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        s_data->Add(hcorr_e2c[bin],"E");
        l_data->AddEntry(hcorr_e2c[bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]),"lf");
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
    gPad->SetLogx(1);
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print(Form("./plots/fullcorre2c_initer%_relerrorleq%.2f_logbinning.pdf",niter,corr_rel_error));
    
    // Draw Linear binning distribution
    s_data = new THStack();
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        s_data->Add(hcorr_e2c_l[bin],"E");
    }
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    c->Print(Form("./plots/fullcorre2c_initer%_relerrorleq%.2f_linbinning.pdf",niter,corr_rel_error));
}
