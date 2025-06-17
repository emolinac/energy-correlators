#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_jes_chisquare(const int nbin = 50, double ptratio_min = 0.4 , double ptratio_max = 1.6, bool do_print = true)
{
    // Open the necessary files
    TFile* f = new TFile((output_folder+namef_ntuple_jes_jer).c_str());
    
    // Get the corresponding Ntuples
    TNtuple* ntuple_jes_data = (TNtuple*) f->Get((name_ntuple_jes_data).c_str());
    TNtuple* ntuple_jes_reco = (TNtuple*) f->Get((name_ntuple_jes_reco).c_str());
    
    // Define the necessary histograms to calculate purity
    TH1F* hdata_nojec[Nbin_jet_pt]; 
    TH1F* hreco_nojec[Nbin_jet_pt];  
    TH1F* hreco_newjec[Nbin_jet_pt];  
    TH1F* hbetastar_balance[Nbin_jet_pt];  
    TH1F* hbetastar_chisquare[Nbin_jet_pt];

    TGraph* chisquare_graph[Nbin_jet_pt];

    THStack* hs[Nbin_jet_pt];
    TLegend* l[Nbin_jet_pt];
    TCanvas* c = new TCanvas("c","",1500,500);
    c->Draw();
    c->Divide(3,1);
    
    const double beta_star_init = 0.9;
    const double beta_star_end  = 1.1;
    const int    beta_star_bins = 60;
    double beta_star_step = (beta_star_end - beta_star_init)/beta_star_bins;
    
    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        c->cd(bin+1);
        
        hs[bin] = new THStack();
        l[bin]  = new TLegend(gPad->GetLeftMargin()+0.01,0.8,gPad->GetLeftMargin()+0.26,0.9,Form(" %.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]));
        hdata_nojec[bin]         = new TH1F(Form("hdata_nojec[%i]",bin)        ,"", nbin           , ptratio_min    , ptratio_max); 
        hreco_nojec[bin]         = new TH1F(Form("hreco_nojec[%i]",bin)        ,"", nbin           , ptratio_min    , ptratio_max); 
        hreco_newjec[bin]        = new TH1F(Form("hreco_newjec[%i]",bin)       ,"", nbin           , ptratio_min    , ptratio_max); 
        hbetastar_balance[bin]   = new TH1F(Form("hbetastar_balance[%i]",bin)  ,"", beta_star_bins , beta_star_init , beta_star_end); 
        hbetastar_chisquare[bin] = new TH1F(Form("hbetastar_chisquare[%i]",bin),"", beta_star_bins , beta_star_init , beta_star_end); 
    
        set_histogram_style(hbetastar_balance[bin]  , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size);
        set_histogram_style(hbetastar_chisquare[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size);
        
        // Undo the JEC
        ntuple_jes_data->Project(Form("hdata_nojec[%i]",bin),"(jet_pt/z_pt)/jet_jec_cor",pair_jetpt_cut[bin]);
        ntuple_jes_reco->Project(Form("hreco_nojec[%i]",bin),"(jet_pt/z_pt)/jet_jec_cor",pair_jetpt_cut[bin]);

        for(int beta_star_bin = 0 ; beta_star_bin < beta_star_bins ; beta_star_bin++)
        {
            double beta_star = beta_star_init + beta_star_bin*beta_star_step;
            ntuple_jes_reco->Project(Form("hreco_newjec[%i]",bin),Form("%f*(jet_pt/z_pt)/jet_jec_cor",beta_star),pair_jetpt_cut[bin]);
            
            double delta_mean = abs(hreco_newjec[bin]->GetMean() - hdata_nojec[bin]->GetMean());
            double chisquare  = hreco_newjec[bin]->Chi2Test(hdata_nojec[bin],"CHI2");

            hbetastar_balance[bin]->SetBinContent(beta_star_bin+1, delta_mean);
            hbetastar_balance[bin]->SetBinError(beta_star_bin+1,sqrt(hreco_newjec[bin]->GetMeanError()*hreco_newjec[bin]->GetMeanError() + hdata_nojec[bin]->GetMeanError()*hdata_nojec[bin]->GetMeanError()));
            hbetastar_chisquare[bin]->SetBinContent(beta_star_bin+1, chisquare);
        }

        double beta_star_min = hbetastar_balance[bin]->GetBinCenter(hbetastar_balance[bin]->GetMinimumBin());

        hbetastar_balance[bin]->Draw();
        hbetastar_balance[bin]->SetMinimum(0);
        hbetastar_balance[bin]->SetMaximum(0.1);
        if(bin==0) hbetastar_balance[bin]->SetTitle(";#beta;");
        else hbetastar_balance[bin]->SetTitle(";#beta;");
        l[bin]->Clear();
        l[bin]->SetHeader(Form(" %.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]));
        l[bin]->AddEntry(hbetastar_balance[bin],Form("#Delta Mean(p^{jet}_{t}/p^{Z}_{t}), #beta^{*}=%.4f",beta_star_min),"p");
        l[bin]->Draw("SAME");
    }

    if(do_print) c->Print(Form("./plots/jes_beta_balance.pdf"));

    for(int bin = 0 ; bin < Nbin_jet_pt ; bin++)
    {
        c->cd(bin+1);
        double beta_star_min = hbetastar_chisquare[bin]->GetBinCenter(hbetastar_chisquare[bin]->GetMinimumBin());

        hbetastar_chisquare[bin]->Draw("E");
        // hbetastar_chisquare[bin]->SetMinimum(0);
        // hbetastar_chisquare[bin]->SetMaximum(0.1);
        if(bin==0) hbetastar_chisquare[bin]->SetTitle(";#beta;");
        else hbetastar_chisquare[bin]->SetTitle(";#beta;");
        l[bin]->Clear();
        l[bin]->SetHeader(Form(" %.1f<p^{jet}_{t}<%.1f GeV",jet_pt_binning[bin],jet_pt_binning[bin+1]));
        l[bin]->AddEntry(hbetastar_chisquare[bin],Form("#chi^{2}(p^{jet}_{t}/p^{Z}_{t}), #beta^{*}=%.4f",beta_star_min),"p");
        l[bin]->Draw("SAME");
    }

    if(do_print) c->Print(Form("./plots/jes_beta_chisquare.pdf"));
}