#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

std::string out_plot = "./plots/closure-test-npairs-eecs.pdf";
void macro_print_npair_eec_ct_plots()
{
        TFile* f = new TFile((output_folder + "histos_eec_3dcorr_rl_jetpt_weightpt_niter8_niterjet4_ct.root").c_str());
        
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

        TH1F* hnpair_ct[nbin_jet_pt];
        TH1F* heec_ct[nbin_jet_pt];

        for(int bin = 0 ; bin < nbin_jet_pt ; bin ++) {
                hnpair_ct[bin] = (TH1F*) f->Get(Form("pseudodata_to_truth_npair%i",bin));
                set_histogram_style(hnpair_ct[bin], corr_marker_color_jet_pt[bin], std_line_width-1, std_marker_style_jet_pt[bin] , std_marker_size);

                heec_ct[bin] = (TH1F*) f->Get(Form("pseudodata_to_truth_eec%i",bin));
                set_histogram_style(heec_ct[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin] , std_marker_size);                
                
                c->cd(bin+1);
                s[bin] = new THStack();
                l[bin] = new TLegend(1-gPad->GetRightMargin()-0.21,gPad->GetBottomMargin()+0.01,1-gPad->GetRightMargin()-0.01,gPad->GetBottomMargin()+0.12);

                s[bin]->Add(hnpair_ct[bin],"E");
                s[bin]->Add(heec_ct[bin],"E1");
                
                s[bin]->Draw("NOSTACK");
                
                s[bin]->SetTitle(Form("%.1f<p^{jet}_{t}(GeV)<%.1f;R_{L};N_{pair}(R_{L})",jet_pt_binning[bin],jet_pt_binning[bin+1]));
                s[bin]->SetMaximum(1.4);
                s[bin]->SetMinimum(0.6);
                s[bin]->GetXaxis()->SetRangeUser(unfolding_rl_nominal_binning[1],unfolding_rl_nominal_binning[nbin_rl_nominal_unfolding-1]);

                l[bin]->AddEntry(hnpair_ct[bin],"N_{pair}","p");
                l[bin]->AddEntry(heec_ct[bin]  ,"EEC"     ,"p");
                
                gPad->SetLogx(1);
                
                l[bin]->Draw("SAME");    
        }

        c->Print(out_plot.c_str());
}
