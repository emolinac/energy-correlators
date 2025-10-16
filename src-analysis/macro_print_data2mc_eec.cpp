#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

void macro_print_data2mc_eec()
{
        TFile* fmc   = new TFile((output_folder + "histos_eec_3dcorr_rl_jetpt_weightpt_niter8_niterjet4_ct.root").c_str());
        TFile* fdata = new TFile((output_folder + "histos_eec_3dcorr_rl_jetpt_weightpt_niter8_niterjet4.root").c_str());
        
        TCanvas* c = new TCanvas("c","",1800,600);
        c->Draw();
        c->Divide(3,1);
        
        TLatex* tex = new TLatex();
        tex->SetTextColorAlpha(16,0.3);
        tex->SetTextSize(0.1991525);
        tex->SetTextAngle(26.15998);
        tex->SetLineWidth(2);

        THStack* s[3];
        TLegend* l[3];

        TH1F* heec_data[nbin_jet_pt];
        TH1F* heec_mc[nbin_jet_pt];

        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                heec_data[bin] = (TH1F*) fdata->Get(Form("hcorr_eec%i",bin));
                set_histogram_style(heec_data[bin], corr_marker_color_jet_pt[bin], std_line_width-2, std_marker_style_jet_pt[bin] , std_marker_size);

                heec_mc[bin] = (TH1F*) fmc->Get(Form("hcorr_eec_truth%i",bin));
                set_histogram_style(heec_mc[bin], corr_marker_color_jet_pt[bin], std_line_width-2, corr_marker_style_jet_pt[bin] , std_marker_size);                
                
                c->cd(bin+1);
                s[bin] = new THStack();
                l[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,gPad->GetBottomMargin()+0.01,1-gPad->GetRightMargin()-0.01,gPad->GetBottomMargin()+0.12);

                s[bin]->Add(heec_data[bin],"E");
                s[bin]->Add(heec_mc[bin],"E1");
                
                s[bin]->Draw("NOSTACK");

                s[bin]->SetMaximum(1.2);
                
                s[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};N_{pair}(R_{L})",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s[bin]->GetXaxis()->SetRangeUser(unfolding_rl_nominal_binning[1],unfolding_rl_nominal_binning[nbin_rl_nominal_unfolding-1]);

                l[bin]->AddEntry(heec_data[bin],"Data" ,"p");
                l[bin]->AddEntry(heec_mc[bin]  ,"Truth","p");
                
                gPad->SetLogx(1);
                
                l[bin]->Draw("SAME");    
        }

        c->Print("./plots/data2mc_eec.pdf");
}
