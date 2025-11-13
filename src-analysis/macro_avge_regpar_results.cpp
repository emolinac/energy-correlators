#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"
#include "../include/utils.cpp"
#include "../include/utils.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

void macro_avge_regpar_results(int niter = 4, int niter_jet = 4)
{
        const int nvariations = 3;

        const int effective_nvariations = nvariations - 1;

        TFile* f[nvariations];

        for (int i = 0 ; i < nvariations ; i++) {
                const int niter_syst = niter - 1 + i;

                f[i] = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niter%i_niterjet%i--get-regpar.root", niter_syst, niter_jet)).c_str());
        }

        TFile* fout = new TFile((output_folder + Form("histos_eec_3dcorr_rl_jetpt_weightpt_niterjet%i--get-regpar.root", niter_jet)).c_str(),"RECREATE");
        gROOT->cd();

        TH1F* hcorr_eec[nbin_jet_pt][nvariations]; 
        TH1F* hcorr_tau[nbin_jet_pt][nvariations]; 
        TH1F* hcorr_eec_eqcharge[nbin_jet_pt][nvariations]; 
        TH1F* hcorr_eec_neqcharge[nbin_jet_pt][nvariations]; 
        
        TH1F* hcorr_eec_avge[nbin_jet_pt]; 
        TH1F* hcorr_tau_avge[nbin_jet_pt]; 
        TH1F* hcorr_eec_eqcharge_avge[nbin_jet_pt]; 
        TH1F* hcorr_eec_neqcharge_avge[nbin_jet_pt]; 
        
        // Get the different variations
        double tau_binning[nbin_jet_pt][nbin_rl_nominal + 1];

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                double avge_pt2_jet = (jet_pt_binning[bin + 1] + jet_pt_binning[bin])/2.;
                
                get_tau_binning_from_eec_binning(tau_binning[bin], rl_nominal_binning, avge_pt2_jet);
                
                hcorr_eec_avge[bin] = new TH1F(Form("relerror_eec%i",bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hcorr_tau_avge[bin] = new TH1F(Form("relerror_tau%i",bin), "", nbin_rl_nominal, tau_binning[bin]);
                hcorr_eec_eqcharge_avge[bin] = new TH1F(Form("relerror_eqcheec%i",bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hcorr_eec_neqcharge_avge[bin] = new TH1F(Form("relerror_neqcheec%i",bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);

                for (int var = 0 ; var < nvariations ; var++) {
                        hcorr_eec[bin][var]           = (TH1F*) f[var]->Get(Form("relerror_eec%i",bin));
                        hcorr_tau[bin][var]           = (TH1F*) f[var]->Get(Form("relerror_tau%i",bin));
                        hcorr_eec_eqcharge[bin][var]  = (TH1F*) f[var]->Get(Form("relerror_eqcheec%i",bin));
                        hcorr_eec_neqcharge[bin][var] = (TH1F*) f[var]->Get(Form("relerror_neqcheec%i",bin));

                        hcorr_eec[bin][var]->Multiply(hcorr_eec[bin][var]);
                        hcorr_tau[bin][var]->Multiply(hcorr_tau[bin][var]);
                        hcorr_eec_eqcharge[bin][var]->Multiply(hcorr_eec_eqcharge[bin][var]);
                        hcorr_eec_neqcharge[bin][var]->Multiply(hcorr_eec_neqcharge[bin][var]);

                        hcorr_eec_avge[bin]->Add(hcorr_eec[bin][var]);
                        hcorr_tau_avge[bin]->Add(hcorr_tau[bin][var]);
                        hcorr_eec_eqcharge_avge[bin]->Add(hcorr_eec_eqcharge[bin][var]);
                        hcorr_eec_neqcharge_avge[bin]->Add(hcorr_eec_neqcharge[bin][var]);
                }
                
                hcorr_eec_avge[bin]->Scale(1./(double) effective_nvariations);
                hcorr_tau_avge[bin]->Scale(1./(double) effective_nvariations);
                hcorr_eec_eqcharge_avge[bin]->Scale(1./(double) effective_nvariations);
                hcorr_eec_neqcharge_avge[bin]->Scale(1./(double) effective_nvariations);

                square_root_bins(hcorr_eec_avge[bin]);
                square_root_bins(hcorr_tau_avge[bin]);
                square_root_bins(hcorr_eec_eqcharge_avge[bin]);
                square_root_bins(hcorr_eec_neqcharge_avge[bin]);

                fout->cd();
                hcorr_eec_avge[bin]->Write();
                hcorr_tau_avge[bin]->Write();
                hcorr_eec_eqcharge_avge[bin]->Write();
                hcorr_eec_neqcharge_avge[bin]->Write();
                gROOT->cd();
        }
}
