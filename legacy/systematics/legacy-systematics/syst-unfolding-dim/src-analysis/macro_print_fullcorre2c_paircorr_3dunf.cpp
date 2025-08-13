#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorre2c_paircorr_3dunf(int niter = nominal_niter, bool do_print = true, bool integrate_weight = true)
{
    // Open the necessary files
    TFile* fout        = new TFile((output_folder + namef_histos_paircorr_e2c).c_str(),"RECREATE");
    TFile* fout_linear = new TFile((output_folder + namef_histos_paircorr_e2c).c_str(),"RECREATE");
    gROOT->cd();

    TFile* fcorr = new TFile((output_folder + namef_ntuple_e2c_paircorr).c_str()); 
    if (fcorr->IsZombie()) return;
    
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

    // Unfold the purity corrected data
    TFile* f = new TFile((output_folder + namef_ntuple_e2c_paircorrections).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_correction_reco.c_str());

    float R_L_reco, R_L_truth, jet_pt_reco, jet_pt_truth, weight_pt_reco, weight_pt_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt_reco);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    ntuple->SetBranchAddress("R_L",&R_L_reco);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    ntuple->SetBranchAddress("weight_pt",&weight_pt_reco);
    ntuple->SetBranchAddress("weight_pt_truth",&weight_pt_truth);
    
    // Create histograms with different types of binning
    TH3D* hpurcorr = new TH3D("hpurcorr","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    TH3D* hmeas    = new TH3D("hmeas"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    TH3D* htrue    = new TH3D("htrue"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");
    
    TH3D* hpurcorr_l = new TH3D("hpurcorr_l","",nbin_rl_nominal_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    TH3D* hmeas_l    = new TH3D("hmeas_l"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    TH3D* htrue_l    = new TH3D("htrue_l"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    RooUnfoldResponse* response_l = new RooUnfoldResponse(hmeas_l, htrue_l, "response_l");

    for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        if (abs(R_L_truth-R_L_reco)>rl_resolution) continue;
    
        response->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
        response_l->Fill(R_L_reco, jet_pt_reco, weight_pt_reco, R_L_truth, jet_pt_truth, weight_pt_truth);
    }

    // Fill the purity corrected distributions
    TH2D* hunfolded_ratio_2d   = new TH2D("hunfolded_ratio_2d"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH3D* hunfolded_ratio_3d   = new TH3D("hunfolded_ratio_3d"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    TH3D* hpuritycorrected     = new TH3D("hpuritycorrected"    ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    TH3D* hpuritycorrected_ref    = new TH3D("hpuritycorrected_ref"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    
    TH2D* hunfolded_ratio_2d_l = new TH2D("hunfolded_ratio_2d_l","",nbin_rl_nominal_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH3D* hunfolded_ratio_3d_l = new TH3D("hunfolded_ratio_3d_l","",nbin_rl_nominal_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    TH3D* hpuritycorrected_l   = new TH3D("hpuritycorrected_l"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    TH3D* hpuritycorrected_ref_l  = new TH3D("hpuritycorrected_ref_l" ,"",nbin_rl_nominal_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning,nbin_weight_unfolding,weight_unfoldingbinning);
    
    ntuple_data->Project("hpuritycorrected" ,"weight_pt:jet_pt:R_L","purity");
    ntuple_data->Project("hpuritycorrected_ref","weight_pt:jet_pt:R_L","purity");
    ntuple_data->Project("hpuritycorrected_l" ,"weight_pt:jet_pt:R_L","purity");
    ntuple_data->Project("hpuritycorrected_ref_l","weight_pt:jet_pt:R_L","purity");
    
    // Unfold the purity corrected pairs
    RooUnfoldBayes unfold(response, hpuritycorrected, niter);
    TH3D* hunfolded_bayes = (TH3D*) unfold.Hunfold();
    hunfolded_ratio_3d->Divide(hunfolded_bayes,hpuritycorrected_ref,1,1);
    
    TH2D* hunfolded_bayes_rl_jet_pt   = (TH2D*) hunfolded_bayes->Project3D("yx");
    TH2D* hpuritycorrected_ref_rl_jet_pt = (TH2D*) hpuritycorrected_ref->Project3D("yx");
    hunfolded_ratio_2d->Divide(hunfolded_bayes_rl_jet_pt,hpuritycorrected_ref_rl_jet_pt,1,1);
    
    RooUnfoldBayes unfold_l(response_l, hpuritycorrected_l, niter);
    TH3D* hunfolded_bayes_l = (TH3D*) unfold_l.Hunfold();
    hunfolded_ratio_3d_l->Divide(hunfolded_bayes_l,hpuritycorrected_ref_l,1,1);
    
    TH2D* hunfolded_bayes_rl_jet_pt_l   = (TH2D*) hunfolded_bayes_l->Project3D("yx");
    TH2D* hpuritycorrected_ref_rl_jet_pt_l = (TH2D*) hpuritycorrected_ref_l->Project3D("yx");
    hunfolded_ratio_2d_l->Divide(hunfolded_bayes_rl_jet_pt_l,hpuritycorrected_ref_rl_jet_pt_l,1,1);
    

    TH1F* hcorr_jet[nbin_jet_pt];
    TH1F* hcorr_e2c[nbin_jet_pt]; 
    TH1F* hcorr_e2c_nounf[nbin_jet_pt]; 
    TH1F* hcorr_e2c_l[nbin_jet_pt]; 
    TH1F* hcorr_e2c_nounf_l[nbin_jet_pt]; 
    
    TCanvas* c = new TCanvas("c", "", 1920, 1080);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    gStyle->SetPaintTextFormat("4.2f");
    hunfolded_ratio_2d->Draw("col text");
    hunfolded_ratio_2d->SetTitle("Purity Corrected Unfolded/Purity Corrected;R_{L};p_{T,jet} (GeV)");
    // hunfolded_ratio_2d->GetXaxis()->SetRangeUser(rl_min, rl_max);
    hunfolded_ratio_2d->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    if (do_print) c->Print(Form("./plots/unfolded3d_niter%i_ratio.pdf",niter));

    hunfolded_ratio_2d_l->Draw("col text");
    hunfolded_ratio_2d_l->SetTitle("Purity Corrected Unfolded/Purity Corrected;R_{L};p_{T,jet} (GeV)");
    // hunfolded_ratio_2d_l->GetXaxis()->SetRangeUser(rl_min, rl_max);
    hunfolded_ratio_2d_l->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
    gPad->SetLogx(0);
    gPad->SetLogy(1);
    if (do_print) c->Print(Form("./plots/unfolded3d_niter%i_ratio_linbinning.pdf",niter));
    
    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);

    // Fill the histograms
    for (int bin = 0 ; bin < nbin_jet_pt ; bin++)
    {
        hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)  ,"", 1,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
        
        hcorr_e2c[bin]          = new TH1F(Form("hcorr_e2c%i",bin)         ,"", nbin_rl_nominal,rl_nominal_binning);
        hcorr_e2c_nounf[bin]    = new TH1F(Form("hcorr_e2c_nounf%i",bin)   ,"", nbin_rl_nominal,rl_nominal_binning);
        hcorr_e2c_l[bin]        = new TH1F(Form("hcorr_e2c_l%i",bin)       ,"", nbin_rl_nominal       ,rl_binning   );
        hcorr_e2c_nounf_l[bin]  = new TH1F(Form("hcorr_e2c_nounf_l%i",bin) ,"", nbin_rl_nominal       ,rl_binning   );
        set_histogram_style(hcorr_e2c[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        set_histogram_style(hcorr_e2c_l[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
    
        for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++)
        {
            ntuple_data->GetEntry(entry);

            if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) continue;
            if (efficiency <= 0 || efficiency > 1) efficiency = 1;//continue;
            if (purity <= 0 || purity > 1) purity = 1;//continue;
            
            double unfolding_weight   = (integrate_weight) ?  hunfolded_ratio_2d->GetBinContent(hunfolded_ratio_2d->FindBin(R_L,jet_pt)): 
                                                              hunfolded_ratio_3d->GetBinContent(hunfolded_ratio_3d->FindBin(R_L,jet_pt,weight_pt));

            double unfolding_weight_l = (integrate_weight) ?  hunfolded_ratio_2d_l->GetBinContent(hunfolded_ratio_2d_l->FindBin(R_L,jet_pt)):
                                                              hunfolded_ratio_3d_l->GetBinContent(hunfolded_ratio_3d_l->FindBin(R_L,jet_pt,weight_pt));
            
            if (unfolding_weight <= 0) unfolding_weight = 1;
            if (unfolding_weight_l <= 0) unfolding_weight_l = 1;

            hcorr_e2c[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
            hcorr_e2c_nounf[bin]->Fill(R_L,purity*weight_pt/efficiency);
            hcorr_e2c_l[bin]->Fill(R_L,purity*unfolding_weight_l*weight_pt/efficiency);
            hcorr_e2c_nounf_l[bin]->Fill(R_L,purity*weight_pt/efficiency);
        }

        // Normalize the distributions
        ntuple_jet->Project(Form("hcorr_jet%i" ,bin), "jet_pt",jet_full_corr[bin]);

        // Log binning
        hcorr_e2c[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
        hcorr_e2c_nounf[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
        
        // Linear binning
        hcorr_e2c_l[bin]->Scale(1./hcorr_jet[bin]->Integral());
        hcorr_e2c_nounf_l[bin]->Scale(1./hcorr_jet[bin]->Integral());

        fout->cd();
        hcorr_e2c[bin]->Write();
        hcorr_e2c_nounf[bin]->Write();
        fout_linear->cd();
        hcorr_e2c_l[bin]->Write();
        hcorr_e2c_nounf_l[bin]->Write();
        gROOT->cd();
    }

    // Draw Log binning distributions
    for (int bin = 0 ; bin < nbin_jet_pt ; bin++)
    {
        s_data->Add(hcorr_e2c[bin],"E");
        l_data->AddEntry(hcorr_e2c[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lf");
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
    l_data->Draw("SAME");
    gPad->SetLogx(1);
    gPad->SetLogy(0);
    
    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    if (do_print) c->Print(Form("./plots/paircorre2c_niter%i_3dunf_wintegrated%s.pdf",niter,integrate_weight?"true":"false"));
    
    // Draw Linear binning distribution
    s_data = new THStack();
    for (int bin = 0 ; bin < nbin_jet_pt ; bin++)
    {
        s_data->Add(hcorr_e2c_l[bin],"E");
    }
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";R_{L};#Sigma_{EEC}(R_{L})");
    l_data->Draw("SAME");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    
    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    if (do_print) c->Print(Form("./plots/paircorre2c_niter%i_linbinning_3dunf_wintegrated%s.pdf",niter,integrate_weight?"true":"false"));
}
