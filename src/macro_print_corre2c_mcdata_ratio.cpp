#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_corre2c_mcdata_ratio()
{
    // Open the necessary files
    TFile* fcorr = new TFile((output_folder+namef_ntuple_e2c_corr).c_str()); 
    TFile* fmc   = new TFile((output_folder+namef_ntuple_mc_e2c).c_str());
    if(fcorr->IsZombie()) return;
    
    TFile* ftest = new TFile("./test.root","RECREATE");
    gROOT->cd();

    TNtuple* ntuple_data   = (TNtuple*) fcorr->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_jet    = (TNtuple*) fcorr->Get((name_ntuple_corrjet).c_str());
    TNtuple* ntuple_mc     = (TNtuple*) fmc->Get((name_ntuple_mc).c_str());
    TNtuple* ntuple_mc_jet = (TNtuple*) fmc->Get((name_ntuple_mc_jet).c_str());

    // UNFOLDING FIRST
    TFile* f = new TFile((output_folder+namef_ntuple_e2c_pairpurity).c_str());
    TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_purity.c_str());

    float R_L, R_L_truth, jet_pt, jet_pt_truth, weight, weight_truth;
    ntuple->SetBranchAddress("jet_pt",&jet_pt);
    ntuple->SetBranchAddress("jet_pt_truth",&jet_pt_truth);
    ntuple->SetBranchAddress("R_L",&R_L);
    ntuple->SetBranchAddress("R_L_truth",&R_L_truth);
    ntuple->SetBranchAddress("weight",&weight);
    ntuple->SetBranchAddress("weight_truth",&weight_truth);
    
    // Create histograms with the respective true and matched reco 
    TH2D* hpurcorr = new TH2D("hpurcorr","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    TH2D* hmeas    = new TH2D("hmeas","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    TH2D* htrue    = new TH2D("htrue","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    RooUnfoldResponse* response = new RooUnfoldResponse(hmeas, htrue, "response");

    for(int evt = 0 ; evt < ntuple->GetEntries() ; evt++)
    {
        // Access entry of ntuple
        ntuple->GetEntry(evt);
        if(jet_pt<unfolding_jetpt_binning[0]||jet_pt>unfolding_jetpt_binning[4]) continue;
        if(R_L_truth==-999) continue;
    
        response->Fill(R_L, jet_pt, R_L_truth, jet_pt_truth);
    }

    // Fill the purity corrected distributions
    TH2D* hunfolded_ratio   = new TH2D("hunfolded_ratio","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    TH2D* hpuritycorrected  = new TH2D("hpuritycorrected","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    TH2D* hpuritycorrected2 = new TH2D("hpuritycorrected2","",Nbin_R_L+2,unfolding_rl_binning,Nbin_jet_pt+2,unfolding_jetpt_binning);
    ntuple_data->Project("hpuritycorrected" ,"jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    ntuple_data->Project("hpuritycorrected2","jet_pt:R_L",pair_purity_corr_singletrack_weightpt);
    
    RooUnfoldBayes unfold(response, hpuritycorrected, 10);
    TH2D* hunfolded_bayes = (TH2D*) unfold.Hunfold();

    hunfolded_ratio->Divide(hunfolded_bayes,hpuritycorrected2,1,1);
    
    TH1F* hcorr_data[Nbin_jet_pt]; 
    TH1F* hall_data[Nbin_jet_pt];  
    TH1F* hcorr_jet[Nbin_jet_pt]; 
    TH1F* hall_jet[Nbin_jet_pt];
    TH1F* hallref_jet[Nbin_jet_pt];
    TH1F* hcorrref_jet[Nbin_jet_pt];

    TH1F* hmc[Nbin_jet_pt]; 
    TH1F* hmc_jet[Nbin_jet_pt]; 
    TH1F* hmcref_jet[Nbin_jet_pt];

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
        c->cd(bin+1);

        s_data[bin] = new THStack();
        l_data[bin] = new TLegend();

        hallref_jet[bin]  = new TH1F(Form("hallref_jet[%i]" ,bin) ,"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]); 
        hcorrref_jet[bin] = new TH1F(Form("hcorrref_jet[%i]" ,bin),"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]); 

        hall_jet[bin]     = new TH1F(Form("hall_jet[%i]" ,bin) ,"",Nbin_R_L,rl_binning);
        hcorr_jet[bin]    = new TH1F(Form("hcorr_jet[%i]" ,bin),"",Nbin_R_L,rl_binning);
        hcorr_data[bin]   = new TH1F(Form("hcorr_data[%i]",bin),"",Nbin_R_L,rl_binning);
        hall_data[bin]    = new TH1F(Form("hall_data[%i]" ,bin),"",Nbin_R_L,rl_binning);

        hmcref_jet[bin] = new TH1F(Form("hmcref_jet[%i]" ,bin),"",1,jet_pt_binning[bin],jet_pt_binning[bin+1]); 
        hmc[bin]        = new TH1F(Form("hmc[%i]" ,bin)       ,"",Nbin_R_L,rl_binning);
        hmc_jet[bin]    = new TH1F(Form("hmc_jet[%i]" ,bin)   ,"",Nbin_R_L,rl_binning);
        
        set_histogram_style(hcorr_data[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size);
        set_histogram_style(hmc[bin] , std_marker_color_jet_pt[bin] , std_line_width, std_marker_style_jet_pt[bin] , std_marker_size);

        // Project into the histograms
        ntuple_data->Project(Form("hcorr_data[%i]",bin),"R_L",e2c_jetpt_full_corr_singletrack_weightpt[bin]);
        ntuple_data->Project(Form("hall_data[%i]" ,bin),"R_L",e2c_jetpt_cut_weightpt[bin]);
        ntuple_jet->Project(Form("hcorrref_jet[%i]" ,bin),"jet_pt",jet_full_corr[bin]);
        ntuple_jet->Project(Form("hallref_jet[%i]" ,bin),"jet_pt",pair_jetpt_cut[bin]);

        ntuple_mc->Project(Form("hmc[%i]" ,bin),"R_L",e2c_jetpt_cut_weightpt[bin]);
        ntuple_mc_jet->Project(Form("hmcref_jet[%i]" ,bin),"jet_pt",pair_jetpt_cut[bin]);
        for(int rl_bin = 0 ; rl_bin < Nbin_R_L ; rl_bin++)
        {
            hall_jet[bin]->SetBinContent(rl_bin+1,  hallref_jet[bin]->GetBinContent(1)); 
            hcorr_jet[bin]->SetBinContent(rl_bin+1, hcorrref_jet[bin]->GetBinContent(1));
            hmc_jet[bin]->SetBinContent(rl_bin+1, hmcref_jet[bin]->GetBinContent(1));
            hall_jet[bin]->SetBinError(rl_bin+1,  hallref_jet[bin]->GetBinError(1)); 
            hcorr_jet[bin]->SetBinError(rl_bin+1, hcorrref_jet[bin]->GetBinError(1));
            hmc_jet[bin]->SetBinError(rl_bin+1, hmcref_jet[bin]->GetBinError(1));
        }

        hcorr_data[bin]->Divide(hcorr_jet[bin]);
        hall_data[bin]->Divide(hall_jet[bin]);
        hmc[bin]->Divide(hmc_jet[bin]);

        // Apply unfolding results
        for(int bin_rl = 0 ; bin_rl < hcorr_data[bin]->GetNbinsX() ; bin_rl++)
        {
            double content_rl  = hcorr_data[bin]->GetBinContent(bin_rl+1);
            double value_rl    = hcorr_data[bin]->GetBinCenter(bin_rl+1);
            double value_jetpt = jet_pt_binning[bin] + (jet_pt_binning[bin+1]-jet_pt_binning[bin])/2.;

            double unfolding_ratio_content = hunfolded_ratio->GetBinContent(hunfolded_ratio->FindBin(value_rl,value_jetpt));
            double corr_content = content_rl*unfolding_ratio_content;

            hcorr_data[bin]->SetBinContent(bin_rl+1,corr_content);
        }

        hcorr_data[bin]->Scale(1./hcorr_data[bin]->Integral());
        hall_data[bin]->Scale(1./hall_data[bin]->Integral());
        hmc[bin]->Scale(1./hmc[bin]->Integral());

        hcorr_data[bin]->Divide(hmc[bin]);

        s_data[bin]->Add(hcorr_data[bin],"E");
        // s_data[bin]->Add(hmc[bin],"E");
        l_data[bin]->AddEntry(hcorr_data[bin],Form("%.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]),"lf");
        s_data[bin]->Draw("NOSTACK");
        s_data[bin]->SetTitle(";R_{L};Data/MC");
        gPad->SetLogx(1);
        s_data[bin]->SetMinimum(0.8);
        s_data[bin]->SetMaximum(1.2);
        l_data[bin]->Draw("SAME");
    }
    
    // tex->DrawLatexNDC(0.25,0.25,"LHCb Internal");

    // c->Print(Form("./plots/corr_norme2c_jetpt_relerrorleq%.2f_weightpt_mccompratio.pdf",corr_rel_error));
}
