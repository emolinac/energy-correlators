#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_matching_fraction_jetpt()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder+namef_ntuple_e2c_hadroncorrections).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Determine log binnning
    double binning[Nbin_R_L+1];
    determine_log10binning(Nbin_R_L, R_L_min, R_L_max, binning);

    TH1F* hall[Nbin_jet_pt];           
    TH1F* hmatched[Nbin_jet_pt];       
    TH1F* hunmatched[Nbin_jet_pt];     
    TH1F* hhalfunmatched[Nbin_jet_pt]; 
    
    TH1F* hratio_matched[Nbin_jet_pt];       
    TH1F* hratio_unmatched[Nbin_jet_pt];     
    TH1F* hratio_halfunmatched[Nbin_jet_pt]; 

    // MCRECO PLOTS
    THStack* s = new THStack();
    
    // Define the necessary histograms to calculate purity
    for (int i = 0 ; i < Nbin_jet_pt ; i++)
    {
        hall[i]           = new TH1F(Form("hall[%i]",i)          ,"",Nbin_R_L,R_L_min,R_L_max);
        hmatched[i]       = new TH1F(Form("hmatched[%i]",i)      ,"",Nbin_R_L,R_L_min,R_L_max);
        hunmatched[i]     = new TH1F(Form("hunmatched[%i]",i)    ,"",Nbin_R_L,R_L_min,R_L_max);
        hhalfunmatched[i] = new TH1F(Form("hhalfunmatched[%i]",i),"",Nbin_R_L,R_L_min,R_L_max);
        hall[i]->Sumw2();
        hmatched[i]->Sumw2();
        hunmatched[i]->Sumw2();
        hhalfunmatched[i]->Sumw2();
    
        // Project into the histograms
        ntuple_dtrmatch->Project(Form("hall[%i]",i)          ,"R_L",pair_jetpt_cut[i]         );
        ntuple_dtrmatch->Project(Form("hmatched[%i]",i)      ,"R_L",pair_jetpt_signal_cut[i]  );
        ntuple_dtrmatch->Project(Form("hunmatched[%i]",i)    ,"R_L",pair_jetpt_pairbg_cut[i]  );
        ntuple_dtrmatch->Project(Form("hhalfunmatched[%i]",i),"R_L",pair_jetpt_singlebg_cut[i]);

        hratio_matched[i]       = new TH1F(Form("hratio_matched[%i]",i)  ,"",Nbin_R_L,R_L_min,R_L_max);
        hratio_unmatched[i]     = new TH1F(Form("hratio_unmatched[%i]",i),"",Nbin_R_L,R_L_min,R_L_max);
        hratio_halfunmatched[i] = new TH1F(Form("hratio_halfunmatched[%i]",i),"",Nbin_R_L,R_L_min,R_L_max);

        hratio_matched[i]->Divide(hmatched[i],hall[i],1,1,"B");
        hratio_unmatched[i]->Divide(hunmatched[i],hall[i],1,1,"B");
        hratio_halfunmatched[i]->Divide(hhalfunmatched[i],hall[i],1,1,"B");

        set_histogram_style(hratio_matched[i]      , kViolet, std_line_width, std_marker_style_jet_pt[i], std_marker_size);
        set_histogram_style(hratio_unmatched[i]    , kCyan  , std_line_width, std_marker_style_jet_pt[i], std_marker_size);
        set_histogram_style(hratio_halfunmatched[i], kCyan+4, std_line_width, std_marker_style_jet_pt[i], std_marker_size); 

        s->Add(hratio_matched[i]);
        s->Add(hratio_unmatched[i]);
        s->Add(hratio_halfunmatched[i]);
    }

    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);
    
    s->Draw("NOSTACK");
    s->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};",R_L_res));

    gPad->SetLogx(1);
    
    s->SetMaximum(1);

    TLegend* l = new TLegend();
    l->AddEntry(hratio_matched[0]      ,"\% N_{pair} Matched && #Delta R<0.02"  ,"lf");
    l->AddEntry(hratio_halfunmatched[0],"\% N_{pair} One unmatched","lf");
    l->AddEntry(hratio_unmatched[0]    ,"\% N_{pair} Both unmatched","lf");
    l->Draw("SAME");

    TLegend* l2 = new TLegend();
    l2->AddEntry(hratio_matched[0],Form("%.1f<p^{jet}_{t}(GeV)<%.1f",jet_pt_binning[0],jet_pt_binning[1]),"p");
    l2->AddEntry(hratio_matched[1],Form("%.1f<p^{jet}_{t}(GeV)<%.1f",jet_pt_binning[1],jet_pt_binning[2]),"p");
    l2->AddEntry(hratio_matched[2],Form("%.1f<p^{jet}_{t}(GeV)<%.1f",jet_pt_binning[2],jet_pt_binning[3]),"p");
    l2->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/matched_unmatched_jetpt_npair_deltarleq%.3f.pdf",R_L_res));    
}