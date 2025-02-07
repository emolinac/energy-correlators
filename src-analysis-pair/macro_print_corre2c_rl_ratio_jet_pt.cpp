#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_corre2c_rl_ratio_jet_pt()
{
    // Open the necessary files
    TFile* fdata       = new TFile((output_folder+namef_ntuple_e2c).c_str());
    TFile* fefficiency = new TFile((output_folder+namef_ntuple_e2c_efficiency).c_str());
    TFile* fpurity     = new TFile((output_folder+namef_ntuple_e2c_purity).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_data            = (TNtuple*) fdata->Get((name_ntuple_data).c_str());
    TNtuple* ntuple_mc              = (TNtuple*) fdata->Get((name_ntuple_mc).c_str());
    TNtuple* ntuple_efficiency_reco = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_reco).c_str());
    TNtuple* ntuple_efficiency_mc   = (TNtuple*) fefficiency->Get((name_ntuple_efficiency_mc).c_str());
    TNtuple* ntuple_purity          = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Determine log binnning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    // Define the necessary histograms to calculate efficiency
    TH1F* hsig_eff[Nbin_jet_pt]; TH1F* hall_eff[Nbin_jet_pt]; TH1F* hefficiency[Nbin_jet_pt];
    TH1F* hsig_pur[Nbin_jet_pt];
    TH1F* hall_pur[Nbin_jet_pt];
    TH1F* hpurity[Nbin_jet_pt];

    // Define MC plots
    TH1F* hmc[Nbin_jet_pt];

    // Define the necessary histograms to show data and corrected data
    TH1F* hcorr_data[Nbin_jet_pt];
    TH1F* hall_data[Nbin_jet_pt];

    // DATA PLOTS
    THStack* s_data = new THStack();
    TLegend* l_data = new TLegend();

    for(int jet_pt_bin = 0 ; jet_pt_bin < Nbin_jet_pt ; jet_pt_bin++)
    {
        hsig_eff[jet_pt_bin]    = new TH1F(Form("hsig_eff[%i]",jet_pt_bin)   ,"",Nbin_R_L,R_L_min,R_L_max);
        hall_eff[jet_pt_bin]    = new TH1F(Form("hall_eff[%i]",jet_pt_bin)   ,"",Nbin_R_L,R_L_min,R_L_max);
        hsig_pur[jet_pt_bin]    = new TH1F(Form("hsig_pur[%i]",jet_pt_bin)   ,"",Nbin_R_L,R_L_min,R_L_max);
        hall_pur[jet_pt_bin]    = new TH1F(Form("hall_pur[%i]",jet_pt_bin)   ,"",Nbin_R_L,R_L_min,R_L_max);
        hefficiency[jet_pt_bin] = new TH1F(Form("hefficiency[%i]",jet_pt_bin),"",Nbin_R_L,R_L_min,R_L_max);
        hpurity[jet_pt_bin]     = new TH1F(Form("hpurity[%i]",jet_pt_bin)    ,"",Nbin_R_L,R_L_min,R_L_max);

        hsig_eff[jet_pt_bin]->Sumw2();
        hall_eff[jet_pt_bin]->Sumw2();
        hsig_pur[jet_pt_bin]->Sumw2();
        hall_pur[jet_pt_bin]->Sumw2();
    
        hcorr_data[jet_pt_bin] = new TH1F(Form("hcorr_data[%i]",jet_pt_bin),"",Nbin_R_L,R_L_min,R_L_max);
        hall_data[jet_pt_bin] = new TH1F(Form("hall_data[%i]",jet_pt_bin),"",Nbin_R_L,R_L_min,R_L_max);

        hmc[jet_pt_bin] = new TH1F(Form("hmc[%i]",jet_pt_bin),"",Nbin_R_L,R_L_min,R_L_max);

        hcorr_data[jet_pt_bin]->Sumw2();
        hall_data[jet_pt_bin]->Sumw2();

        ntuple_efficiency_reco->Project(Form("hsig_eff[%i]",jet_pt_bin),"R_L",pair_jetpt_signal_cut[jet_pt_bin]);
        ntuple_efficiency_mc->Project(Form("hall_eff[%i]",jet_pt_bin),"R_L",pair_jetpt_cut[jet_pt_bin]);
        ntuple_purity->Project(Form("hsig_pur[%i]",jet_pt_bin),"R_L",pair_jetpt_signal_cut[jet_pt_bin]);
        ntuple_purity->Project(Form("hall_pur[%i]",jet_pt_bin),"R_L",pair_jetpt_cut[jet_pt_bin]);
        ntuple_data->Project(Form("hcorr_data[%i]",jet_pt_bin),"R_L",e2c_jetpt_cut[jet_pt_bin]);
        ntuple_data->Project(Form("hall_data[%i]",jet_pt_bin),"R_L", e2c_jetpt_cut[jet_pt_bin]);
        ntuple_mc->Project(Form("hmc[%i]",jet_pt_bin),"R_L", e2c_jetpt_cut[jet_pt_bin]);

        // Normalize reference
        hmc[jet_pt_bin]->Scale(1./hmc[jet_pt_bin]->Integral());

        // Get correction factors
        hefficiency[jet_pt_bin]->Divide(hsig_eff[jet_pt_bin],hall_eff[jet_pt_bin],1,1,"B");
        hpurity[jet_pt_bin]->Divide(hsig_pur[jet_pt_bin],hall_pur[jet_pt_bin],1,1,"B");
    
        // Apply corrections
        hcorr_data[jet_pt_bin]->Divide(hefficiency[jet_pt_bin]);
        hcorr_data[jet_pt_bin]->Multiply(hpurity[jet_pt_bin]);
        hcorr_data[jet_pt_bin]->Scale(1./hcorr_data[jet_pt_bin]->Integral());

        hcorr_data[jet_pt_bin]->Divide(hmc[jet_pt_bin]);

        set_histogram_style(hcorr_data[jet_pt_bin], corr_marker_color_jet_pt[jet_pt_bin], std_line_width, corr_marker_style_jet_pt[jet_pt_bin], std_marker_size);
        s_data->Add(hcorr_data[jet_pt_bin]);
        //s_data->Add(hall_data[jet_pt_bin]);

        l_data->AddEntry(hcorr_data[jet_pt_bin],Form("%.1f<Jet P_{T}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1]),"lpf");
        //l_data->AddEntry(hall_data[jet_pt_bin],Form("Data : %.1f<Jet P_{T}(GeV)<%.1f",jet_pt_binning[jet_pt_bin],jet_pt_binning[jet_pt_bin+1])      ,"lpf");
    }
    
    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    tex->SetTextColorAlpha(16,0.3);
    tex->SetTextSize(0.1991525);
    tex->SetTextAngle(26.15998);
    tex->SetLineWidth(2);

    // MCRECO PLOTS
    gPad->SetLogx(1);
   
    s_data->Draw("NOSTACK");
    s_data->GetXaxis()->SetRangeUser(R_L_min,1);
    s_data->SetMaximum(1.4);
    s_data->SetMinimum(0.6);    
    s_data->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};Norm. E2C (Meas/MC)",R_L_res));
    l_data->Draw("SAME");

    tex->DrawLatexNDC(0.2,0.2,"LHCb Internal");
    
    c->Print(Form("../plots/corr_e2c_ratio_jetpt_deltarleq%.3f.pdf",R_L_res));
}