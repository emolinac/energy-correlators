#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/utils-algorithms.h"
#include "../../include/utils-visual.h"

void macro_print_matching_fraction_e2c()
{
    // Open the necessary files
    TFile* fpurity = new TFile((output_folder + namef_ntuple_e2c_hadroncorrections).c_str());

    // Get the corresponding Ntuples
    TNtuple* ntuple_dtrmatch = (TNtuple*) fpurity->Get((name_ntuple_purity).c_str());

    // Determine log binnning
    double binning[nbin_rl+1];
    determine_log10binning(nbin_rl, rl_min, rl_max, binning);

    // Define the necessary histograms to calculate purity
    TH1F* hall           = new TH1F("hall"      ,"",nbin_rl,rl_min, rl_max);
    TH1F* hmatched       = new TH1F("hmatched"  ,"",nbin_rl,rl_min, rl_max);
    TH1F* hunmatched     = new TH1F("hunmatched","",nbin_rl,rl_min, rl_max);
    TH1F* hhalfunmatched = new TH1F("hhalfunmatched","",nbin_rl,rl_min, rl_max);
    hall->Sumw2();
    hmatched->Sumw2();
    hunmatched->Sumw2();
    hhalfunmatched->Sumw2();
    
    // Project into the histograms
    ntuple_dtrmatch->Project("hall"          ,"R_L",e2c_cut);
    ntuple_dtrmatch->Project("hmatched"      ,"R_L",e2c_signal_cut);
    ntuple_dtrmatch->Project("hunmatched"    ,"R_L",e2c_pairbg_cut);
    ntuple_dtrmatch->Project("hhalfunmatched","R_L",e2c_singlebg_cut);

    TH1F* hratio_matched       = new TH1F("hratio_matched"  ,"",nbin_rl,rl_min, rl_max);
    TH1F* hratio_unmatched     = new TH1F("hratio_unmatched","",nbin_rl,rl_min, rl_max);
    TH1F* hratio_halfunmatched = new TH1F("hratio_halfunmatched","",nbin_rl,rl_min, rl_max);
    
    hratio_matched->Divide(hmatched,hall,1,1,"B");
    hratio_unmatched->Divide(hunmatched,hall,1,1,"B");
    hratio_halfunmatched->Divide(hhalfunmatched,hall,1,1,"B");


    set_histogram_style(hratio_matched, kViolet, std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hratio_unmatched, kCyan  , std_line_width, std_marker_style, std_marker_size);
    set_histogram_style(hratio_halfunmatched, kCyan+4  , std_line_width, std_marker_style, std_marker_size);

    TCanvas* c = new TCanvas("c","",800,600);
    c->Draw();

    TLatex* tex = new TLatex();
    set_lhcb_watermark_properties(tex);
    
    // MCRECO PLOTS
    THStack* s = new THStack();
    s->Add(hratio_matched);
    s->Add(hratio_unmatched);
    s->Add(hratio_halfunmatched);
    s->Draw("NOSTACK");
    s->SetTitle(Form("#Delta R_{L}(truth-reco)<%.3f;R_{L};",rl_resolution));

    gPad->SetLogx(1);
    
    s->SetMaximum(1);

    TLegend* l = new TLegend();
    l->AddEntry(hratio_matched      ,"\% E2C Matched && #Delta R<0.02"  ,"lpf");
    l->AddEntry(hratio_halfunmatched,"\% E2C One unmatched","lpf");
    l->AddEntry(hratio_unmatched    ,"\% E2C Both unmatched","lpf");
    l->Draw("SAME");

    tex->DrawLatexNDC(0.3,0.3,"simulations");

    c->Print(Form("./plots/matched_unmatched_e2c_deltarleq%.3f.pdf",rl_resolution));
    
}