#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_fullcorrchargede2c_paircorr_2dunf(int niter = nominal_niter, bool do_print = true)
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
    float R_L, jet_pt, weight_pt, eq_charge;
    float efficiency, purity, efficiency_relerror, purity_relerror;
    ntuple_data->SetBranchAddress("R_L",&R_L);
    ntuple_data->SetBranchAddress("jet_pt",&jet_pt);
    ntuple_data->SetBranchAddress("weight_pt",&weight_pt);
    ntuple_data->SetBranchAddress("eq_charge",&eq_charge);
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
    TH2D* hpurcorr = new TH2D("hpurcorr","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH2D* hmeas    = new TH2D("hmeas"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH2D* htrue    = new TH2D("htrue"   ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");
    
    TH2D* hpurcorr_l = new TH2D("hpurcorr_l","",nbin_rl_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH2D* hmeas_l    = new TH2D("hmeas_l"   ,"",nbin_rl_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH2D* htrue_l    = new TH2D("htrue_l"   ,"",nbin_rl_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    RooUnfoldResponse* response_l = new RooUnfoldResponse(hmeas_l, htrue_l, "response_l");

    for (int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        if (abs(R_L_truth-R_L_reco)>rl_resolution) continue;
    
        response->Fill(R_L_reco, jet_pt_reco, R_L_truth, jet_pt_truth);
        response_l->Fill(R_L_reco, jet_pt_reco, R_L_truth, jet_pt_truth);
    }

    // Fill the purity corrected distributions
    TH2D* hunfolded_ratio     = new TH2D("hunfolded_ratio"  ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH2D* hpuritycorrected    = new TH2D("hpuritycorrected" ,"",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH2D* hpuritycorrected2   = new TH2D("hpuritycorrected2","",nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH2D* hunfolded_ratio_l   = new TH2D("hunfolded_ratio_l"  ,"",nbin_rl_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH2D* hpuritycorrected_l  = new TH2D("hpuritycorrected_l" ,"",nbin_rl_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    TH2D* hpuritycorrected2_l = new TH2D("hpuritycorrected2_l","",nbin_rl_unfolding,unfolding_rl_binning,nbin_jet_pt_unfolding,unfolding_jet_pt_binning);
    
    ntuple_data->Project("hpuritycorrected" , "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    ntuple_data->Project("hpuritycorrected2", "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    ntuple_data->Project("hpuritycorrected_l" , "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    ntuple_data->Project("hpuritycorrected2_l", "jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    
    // Unfold the purity corrected pairs
    RooUnfoldBayes unfold(response, hpuritycorrected, niter);
    TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();
    hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected2,1,1);

    RooUnfoldBayes unfold_l(response_l, hpuritycorrected_l, niter);
    TH2D* hunfolded_bayes_l = (TH2D*) unfold_l.Hunfold();
    hunfolded_ratio_l->Divide(hunfolded_bayes_l,hpuritycorrected2_l,1,1);

    TH1F* hcorr_e2c_eqcharge[nbin_jet_pt]; 
    TH1F* hcorr_e2c_neqcharge[nbin_jet_pt]; 

    TH1F* hcorr_jet[nbin_jet_pt];
    TH1F* hcorr_jet_centroid[nbin_jet_pt];
    TH1F* hcorr_e2c[nbin_jet_pt]; 
    TH1F* hcorr_e2c_nounf[nbin_jet_pt]; 
    TH1F* hcorr_e2c_l[nbin_jet_pt]; 
    TH1F* hcorr_e2c_nounf_l[nbin_jet_pt]; 
    TH1F* hcorr_tau[nbin_jet_pt]; 
    TH1F* hcorr_tau_nounf[nbin_jet_pt]; 
    
    TCanvas* c = new TCanvas("c", "", 1920, 1080);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);

    gStyle->SetPaintTextFormat("4.2f");
    hunfolded_ratio->Draw("col text");
    hunfolded_ratio->SetTitle("Purity Corrected Unfolded/Purity Corrected;R_{L};p_{T,jet} (GeV)");
    // hunfolded_ratio->GetXaxis()->SetRangeUser(rl_min, rl_max);
    hunfolded_ratio->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    if (do_print) c->Print(Form("./plots/unfolded2d_niter%i_ratio.pdf",niter));

    hunfolded_ratio_l->Draw("col text");
    hunfolded_ratio_l->SetTitle("Purity Corrected Unfolded/Purity Corrected;R_{L};p_{T,jet} (GeV)");
    hunfolded_ratio_l->GetXaxis()->SetRangeUser(rl_min, rl_max);
    hunfolded_ratio_l->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
    gPad->SetLogx(0);
    gPad->SetLogy(1);
    if (do_print) c->Print(Form("./plots/unfolded2d_niter%i_ratio_linbinning.pdf",niter));

    THStack* s_data     = new THStack();
    TLegend* l_data     = new TLegend(0.4,gPad->GetBottomMargin()+0.01,0.6,0.2+gPad->GetBottomMargin()+0.01);
    TLegend* l_data_tau = new TLegend(gPad->GetLeftMargin()+0.01,1-gPad->GetTopMargin()-0.01,gPad->GetLeftMargin()+0.21,1-gPad->GetTopMargin()-0.21);
    
    TLegend* l_data_chargede2c = new TLegend(1-gPad->GetRightMargin()-0.31,gPad->GetBottomMargin()+0.01,1-gPad->GetRightMargin()-0.01,gPad->GetBottomMargin()+0.21);

    // Fill the histograms
    for (int bin = 0 ; bin < nbin_jet_pt ; bin++)
    {
        hcorr_jet[bin]          = new TH1F(Form("hcorr_jet%i" ,bin)         ,"", 1  ,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 
        hcorr_jet_centroid[bin] = new TH1F(Form("hcorr_jet_centroid%i" ,bin),"", 200,jet_pt_binning[bin],jet_pt_binning[bin + 1]); 

        hcorr_e2c_eqcharge[bin]  = new TH1F(Form("hcorr_e2c_eqcharge%i",bin) ,"", nbin_rl_nominal,rl_nominal_binning );
        hcorr_e2c_neqcharge[bin] = new TH1F(Form("hcorr_e2c_neqcharge%i",bin),"", nbin_rl_nominal,rl_nominal_binning );

        hcorr_e2c[bin]          = new TH1F(Form("hcorr_e2c%i",bin)         ,"", nbin_rl_nominal,rl_nominal_binning );
        hcorr_e2c_nounf[bin]    = new TH1F(Form("hcorr_e2c_nounf%i",bin)   ,"", nbin_rl_nominal,rl_nominal_binning );
        hcorr_tau[bin]          = new TH1F(Form("hcorr_tau%i",bin)         ,"", nbin_rl_nominal,tau_nominal_binning);
        hcorr_tau_nounf[bin]    = new TH1F(Form("hcorr_tau_nounf%i",bin)   ,"", nbin_rl_nominal,tau_nominal_binning);
        hcorr_e2c_l[bin]        = new TH1F(Form("hcorr_e2c_l%i",bin)       ,"", nbin_rl       ,rl_binning    );
        hcorr_e2c_nounf_l[bin]  = new TH1F(Form("hcorr_e2c_nounf_l%i",bin) ,"", nbin_rl       ,rl_binning    );
        
        set_histogram_style(hcorr_e2c_eqcharge[bin]  , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
        set_histogram_style(hcorr_e2c_neqcharge[bin] , corr_marker_color_jet_pt[bin], std_line_width-1, std_marker_style_jet_pt[bin] , std_marker_size+1);
        
        set_histogram_style(hcorr_e2c[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        set_histogram_style(hcorr_tau[bin]  , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
        set_histogram_style(hcorr_e2c_l[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
    
        ntuple_jet->Project(Form("hcorr_jet%i" ,bin)         , "jet_pt",jet_full_corr[bin]);
        ntuple_jet->Project(Form("hcorr_jet_centroid%i" ,bin), "jet_pt",jet_full_corr[bin]);

        double jet_pt_centroid = hcorr_jet_centroid[bin]->GetMean();
        for (int entry = 0 ; entry < ntuple_data->GetEntries() ; entry++)
        {
            ntuple_data->GetEntry(entry);

            if (jet_pt < jet_pt_binning[bin] || jet_pt > jet_pt_binning[bin + 1]) continue;
            if (efficiency <= 0 || efficiency > 1) efficiency = 1;//continue;
            if (purity <= 0 || purity > 1) purity = 1;//continue;
            
            double unfolding_weight = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(R_L,jet_pt));
            if (unfolding_weight <= 0) unfolding_weight = 1;

            double unfolding_weight_l = hunfolded_ratio_l->GetBinContent(hunfolded_ratio_l->FindBin(R_L,jet_pt));
            if (unfolding_weight_l <= 0) unfolding_weight_l = 1;

            hcorr_e2c[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
            hcorr_e2c_nounf[bin]->Fill(R_L,purity*weight_pt/efficiency);
            hcorr_tau[bin]->Fill(R_L*jet_pt_centroid,purity*unfolding_weight*weight_pt/efficiency);
            hcorr_tau_nounf[bin]->Fill(R_L*jet_pt_centroid,purity*weight_pt/efficiency);
            hcorr_e2c_l[bin]->Fill(R_L,purity*unfolding_weight_l*weight_pt/efficiency);
            hcorr_e2c_nounf_l[bin]->Fill(R_L,purity*weight_pt/efficiency);

            // Filling the charged e2cs
            if (eq_charge>0)  hcorr_e2c_eqcharge[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
            else if (eq_charge<0) hcorr_e2c_neqcharge[bin]->Fill(R_L,purity*unfolding_weight*weight_pt/efficiency);
        }

        // Normalize charged e2cs
        hcorr_e2c_eqcharge[bin]->Divide(hcorr_e2c[bin]);
        hcorr_e2c_neqcharge[bin]->Divide(hcorr_e2c[bin]);
        
        // Normalize the distributions
        // Log binning
        hcorr_e2c[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
        hcorr_e2c_nounf[bin]->Scale(1./hcorr_jet[bin]->Integral(),"width");
        
        // Linear binning
        hcorr_e2c_l[bin]->Scale(1./hcorr_jet[bin]->Integral());
        hcorr_e2c_nounf_l[bin]->Scale(1./hcorr_jet[bin]->Integral());

        fout->cd();
        hcorr_e2c[bin]->Write();
        hcorr_e2c_nounf[bin]->Write();
        hcorr_tau[bin]->Write();
        hcorr_tau_nounf[bin]->Write();
        hcorr_e2c_eqcharge[bin]->Write();
        hcorr_e2c_neqcharge[bin]->Write();
        fout_linear->cd();
        hcorr_e2c_l[bin]->Write();
        hcorr_e2c_nounf_l[bin]->Write();
        gROOT->cd();
    }

    // Draw Linear binning distribution
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

    if (do_print) c->Print(Form("./plots/paircorre2c_niter%i_linbinning_2dunf.pdf",niter));

    // Draw Log binning distributions
    s_data = new THStack();
    for (int bin = 0 ; bin < nbin_jet_pt ; bin++)
    {
        s_data->Add(hcorr_tau[bin],"E");
        l_data_tau->AddEntry(hcorr_tau[bin],Form("%.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lf");
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetTitle(";R_{L} #LT p_{T,jet} #GT(GeV);#Sigma_{EEC}(R_{L})");
    l_data_tau->Draw("SAME");
    gPad->SetLogx(1);
    gPad->SetLogy(0);
    
    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    if (do_print) c->Print(Form("./plots/paircorrtau_niter%i_2dunf.pdf",niter));

    s_data = new THStack();
    for (int bin = 0 ; bin < nbin_jet_pt ; bin++)
    {
        s_data->Add(hcorr_e2c_eqcharge[bin] ,"E1 X0 ");
        s_data->Add(hcorr_e2c_neqcharge[bin],"E1 X0 ");
        l_data_chargede2c->AddEntry(hcorr_e2c_eqcharge[bin],Form("eq. charge %.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lfp");
        l_data_chargede2c->AddEntry(hcorr_e2c_neqcharge[bin],Form("op. charge %.1f<p_{T,jet}<%.1f (GeV)",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"lfp");
    }
    
    s_data->Draw("NOSTACK");
    s_data->SetMaximum(1);
    s_data->SetTitle(";R_{L};Charged #Sigma_{EEC}(R_{L})");
    l_data_chargede2c->Draw("SAME");
    gPad->SetLogx(1);
    gPad->SetLogy(0);
    
    tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    if (do_print) c->Print(Form("./plots/paircorrchargede2c_niter%i_2dunf.pdf",niter));

    s_data = new THStack();
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

    if (do_print) c->Print(Form("./plots/paircorre2c_niter%i_2dunf.pdf",niter));
}

