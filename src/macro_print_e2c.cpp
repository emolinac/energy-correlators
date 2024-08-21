#include "../include/analysis-constants.h"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils-algorithms.h"
#include "../include/utils-visual.h"

TCut e2c_cut = "weight*(jet_eta>2.5&&jet_eta<4.&&jet_pt>20&&jet_pt<30&&h1_chi2/h1_ndf<3&&h2_chi2/h2_ndf<3&&h1_p>4&&h1_p<1000&&h2_p>4&&h2_p<1000&&h1_pt>0.250&&h2_pt>0.250&&\
                h1_probnnghost<0.5&&h2_probnnghost<0.5&&TMath::Abs(z0_phi-jet_phi)>7*TMath::Pi()/8.&&TMath::Abs(jet_phi-mum_phi)>0.4&&TMath::Abs(jet_phi-mup_phi)>0.4&&\
                mum_pt>20.&&mup_pt>20.&&mum_eta>2&&mum_eta<4.5&&mup_eta>2&&mup_eta<4.5&&sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))>60&&\
                sqrt(mum_m*mum_m + mup_m*mup_m + 2*(mum_pe*mup_pe - mum_px*mup_px - mum_py*mup_py - mum_pz*mup_pz))<120&&mum_probchi2>0.001&&mup_probchi2>0.001)";

void macro_print_e2c()
{
    // Open ROOT file with ntuple
    TFile* f = new TFile((output_folder+namef_ntuple_e2c).c_str());

    // Get the ntuple
    TNtuple* ntuple = (TNtuple*) f->Get((name_ntuple_data).c_str());

    // Determine binning
    double binning[Nbin_X_L+1];
    determine_log10binning(Nbin_X_L, X_L_min, X_L_max, binning);

    TH1F* h = new TH1F("h","",Nbin_X_L, binning);
    h->Sumw2();

    ntuple->Draw("X_L>>h",e2c_cut);

    set_histogram_style(h, kViolet+2, std_line_width, std_marker_style, std_marker_size);
    h->SetTitle(";X_{L};E2C");
    gPad->SetLogx(1);
    gPad->SetLogy(1);
}