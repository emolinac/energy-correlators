#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorre2c_mc_comp(int niter = 20)
{
    // Open the necessary files
    TFile* fcorr = new TFile((output_folder+namef_ntuple_e2c_corr).c_str()); 
    TFile* fmc   = new TFile((output_folder+namef_ntuple_mc_e2c).c_str());
    if(fcorr->IsZombie()) return;
    
    TNtuple* ntuple_data       = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_jet        = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
    TNtuple* ntuple_mc         = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
    TNtuple* ntuple_mc_jet     = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());
    TNtuple* ntuple_mcreco     = (TNtuple*) fmc->Get((name_ntuple_mcreco).c_str());
    TNtuple* ntuple_mcreco_jet = (TNtuple*) fmc->Get((name_ntuple_mcreco_jet).c_str());

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

    // Set the branches of data
    float R_L_mc, jet_pt_mc, weight_pt_mc;
    ntuple_mc->SetBranchAddress("R_L",&R_L_mc);
    ntuple_mc->SetBranchAddress("jet_pt",&jet_pt_mc);
    ntuple_mc->SetBranchAddress("weight_pt",&weight_pt_mc);
    
    // Set the branches of data
    float R_L_mcreco, jet_pt_mcreco, weight_pt_mcreco;
    ntuple_mcreco->SetBranchAddress("R_L",&R_L_mcreco);
    ntuple_mcreco->SetBranchAddress("jet_pt",&jet_pt_mcreco);
    ntuple_mcreco->SetBranchAddress("weight_pt",&weight_pt_mcreco);
    
    // UNFOLDING FIRST
    TFile* f = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

    float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_reco, weight_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt_reco);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    ntuple->SetBranchAddress("R_L",&R_L_reco);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    ntuple->SetBranchAddress("weight_pt",&weight_reco);
    ntuple->SetBranchAddress("weight_pt_truth",&weight_truth);
    
    // Create histograms with the respective true and matched reco 
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
    
        response->Fill(R_L_reco, jet_pt_reco, weight_reco, R_L_truth, jet_pt_truth, weight_truth);
        response_l->Fill(R_L_reco, jet_pt_reco, weight_reco, R_L_truth, jet_pt_truth, weight_truth);
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

    TH1F* hcorr_e2c[Nbin_jet_pt]; 
    TH1F* hmc[Nbin_jet_pt]; 
    TH1F* hmcreco[Nbin_jet_pt]; 
    TH1F* hcorr_e2c_l[Nbin_jet_pt]; 
    TH1F* hmc_l[Nbin_jet_pt]; 
    TH1F* hmcreco_l[Nbin_jet_pt]; 

    TH1F* hcorr_jet[Nbin_jet_pt];
    TH1F* hmc_jet[Nbin_jet_pt]; 
    TH1F* hmcreco_jet[Nbin_jet_pt]; 
    
    TCanvas* c = new TCanvas("c","",1800,600);
    c->Draw();
    c->Divide(3,1);
    
    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    THStack* s_data[3];
    TLegend* l_data[3];

    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        hcorr_e2c[bin]   = new TH1F(Form("hcorr_e2c[%i]",bin) ,"",Nbin_R_L_logbin,rl_logbinning);
        hmc[bin]         = new TH1F(Form("hmc[%i]" ,bin)      ,"",Nbin_R_L_logbin,rl_logbinning);
        hmcreco[bin]     = new TH1F(Form("hmcreco[%i]" ,bin)  ,"",Nbin_R_L_logbin,rl_logbinning);
        hcorr_e2c_l[bin] = new TH1F(Form("hcorr_e2c_l[%i]",bin) ,"",Nbin_R_L,rl_binning);
        hmc_l[bin]       = new TH1F(Form("hmc_l[%i]" ,bin)      ,"",Nbin_R_L,rl_binning);
        hmcreco_l[bin]   = new TH1F(Form("hmcreco_l[%i]" ,bin)  ,"",Nbin_R_L,rl_binning);
        
        hcorr_jet[bin]   = new TH1F(Form("hcorr_jet[%i]" ,bin)   ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]); 
        hmc_jet[bin]     = new TH1F(Form("hmc_jet[%i]" ,bin)     ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]);
        hmcreco_jet[bin] = new TH1F(Form("hmcreco_jet[%i]" ,bin) ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]);
        
        set_histogram_style(hcorr_e2c[bin]  , corr_marker_color_jet_pt[0], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size);
        set_histogram_style(hmc[bin]        , corr_marker_color_jet_pt[1], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
        set_histogram_style(hmcreco[bin]    , corr_marker_color_jet_pt[2], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
        set_histogram_style(hcorr_e2c_l[bin], corr_marker_color_jet_pt[0], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size);
        set_histogram_style(hmc_l[bin]      , corr_marker_color_jet_pt[1], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
        set_histogram_style(hmcreco_l[bin]  , corr_marker_color_jet_pt[2], std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);
    
        // Fill and normalize data
        for(int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++)
        {
            ntuple_data->GetEntry(entry);

            if(jet_pt<jet_pt_binning[bin]||jet_pt>jet_pt_binning[bin+1]) continue;
            if(efficiency_relerror>corr_rel_error) continue;
            if(purity_relerror>corr_rel_error) continue;
            if(efficiency<=0||efficiency>1) efficiency = 1 ; //continue;
            if(purity<=0||purity>1) purity = 1 ;// continue;
            if(weight_pt>weight_pt_cut[bin]) continue;
            
            double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt,weight_pt));
            if(unfolding_weight<=0) unfolding_weight = 1;

            double unfolding_weight_l = hunfolded_ratio_l->GetBinContent(hunfolded_ratio_l->FindBin(R_L,jet_pt,weight_pt));
            if(unfolding_weight_l<=0) unfolding_weight_l = 1;

            hcorr_e2c[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
            hcorr_e2c_l[bin]->Fill(R_L,purity*unfolding_weight_l*weight_pt/efficiency);
        }
        
        ntuple_jet->Project(Form("hcorr_jet[%i]" ,bin),"jet_pt",jet_full_corr[bin]);
        
        hcorr_e2c[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
        hcorr_e2c_l[bin]->Scale(1./hcorr_jet[bin]->Integral());

        // Fill and normalize MC        
        for(int entry = 0 ; entry < ntuple_mc->GetEntries() ; entry++)
        {
            ntuple_mc->GetEntry(entry);

            if(jet_pt_mc<jet_pt_binning[bin]||jet_pt_mc>jet_pt_binning[bin+1]) continue;
            if(weight_pt_mc>weight_pt_cut[bin]) continue;
            
            hmc[bin]->Fill(R_L_mc,weight_pt_mc);
            hmc_l[bin]->Fill(R_L_mc,weight_pt_mc);
        }
        
        ntuple_mc_jet->Project(Form("hmc_jet[%i]" ,bin),"jet_pt",pair_jetpt_cut[bin]);
        hmc[bin]->Scale(1./hmc_jet[bin]->Integral(),"width");
        hmc_l[bin]->Scale(1./hmc_jet[bin]->Integral());

        // Fill and normalize MCReco
        for(int entry = 0 ; entry < ntuple_mcreco->GetEntries() ; entry++)
        {
            ntuple_mcreco->GetEntry(entry);

            if(jet_pt_mcreco<jet_pt_binning[bin]||jet_pt_mcreco>jet_pt_binning[bin+1]) continue;
            if(weight_pt_mcreco>weight_pt_cut[bin]) continue;

            hmcreco[bin]->Fill(R_L_mcreco,weight_pt_mcreco);
            hmcreco_l[bin]->Fill(R_L_mcreco,weight_pt_mcreco);
        }

        ntuple_mcreco_jet->Project(Form("hmcreco_jet[%i]" ,bin),"jet_pt",pair_jetpt_cut[bin]);
        hmcreco[bin]->Scale(1./hmcreco_jet[bin]->Integral(),"width");
        hmcreco_l[bin]->Scale(1./hmcreco_jet[bin]->Integral());
    }

    // Draw the log binning histos
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin ++)
    {
        c->cd(bin+1);
        s_data[bin] = new THStack();
        l_data[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,1-gPad->GetTopMargin()-0.15,1-gPad->GetRightMargin()-0.01,1-gPad->GetTopMargin()-0.01);

        s_data[bin]->Add(hmc[bin],"E");
        s_data[bin]->Add(hmcreco[bin],"E");
        s_data[bin]->Add(hcorr_e2c[bin],"E");
        s_data[bin]->Draw("NOSTACK");
        s_data[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};#Sigma_{EEC}(R_{L})",jet_pt_binning[bin],jet_pt_binning[bin+1]));
        s_data[bin]->SetMaximum(1.7);
        s_data[bin]->SetMinimum(40E-03);

        l_data[bin]->AddEntry(hmc[bin]      ,"mc"    ,"p");
        l_data[bin]->AddEntry(hmcreco[bin]  ,"mcreco","p");
        l_data[bin]->AddEntry(hcorr_e2c[bin],"data"  ,"p");
        gPad->SetLogx(1);
        l_data[bin]->Draw("SAME");    
    }
    c->Print(Form("./plots/fullcorre2c_niter%i_relerrorleq%.2f_mccomp_logbinning.pdf",niter,corr_rel_error));
    
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin ++)
    {
        c->cd(bin+1);
        s_data[bin]->SetMaximum(1.7);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
    }
    c->Print(Form("./plots/fullcorre2c_niter%i_relerrorleq%.2f_mccomp_logbinning_logscale.pdf",niter,corr_rel_error));
    
    // Draw the log binning histos
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin ++)
    {
        c->cd(bin+1);
        s_data[bin] = new THStack();
        l_data[bin] = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);

        s_data[bin]->Add(hmc_l[bin],"E");
        s_data[bin]->Add(hmcreco_l[bin],"E");
        s_data[bin]->Add(hcorr_e2c_l[bin],"E");
        s_data[bin]->Draw("NOSTACK");
        s_data[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};#Sigma_{EEC}(R_{L})",jet_pt_binning[bin],jet_pt_binning[bin+1]));
        s_data[bin]->SetMaximum(9E-02);
        s_data[bin]->SetMinimum(1E-04);

        l_data[bin]->AddEntry(hmc_l[bin]      ,"mc"    ,"p");
        l_data[bin]->AddEntry(hmcreco_l[bin]  ,"mcreco","p");
        l_data[bin]->AddEntry(hcorr_e2c_l[bin],"data"  ,"p");
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        l_data[bin]->Draw("SAME");    
    }
    
    c->Print(Form("./plots/fullcorre2c_niter%i_relerrorleq%.2f_mccomp_linbinning.pdf",niter,corr_rel_error));
}