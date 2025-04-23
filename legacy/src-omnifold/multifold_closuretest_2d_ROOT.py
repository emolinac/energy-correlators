from omnifold import DataLoader
from omnifold import MultiFold
from omnifold import MLP

import numpy as np
import os
import uproot
import ROOT

from matplotlib import pyplot as plt
# from keras.layers import Dense, Input
# from keras.models import Model
from array import array
from sklearn.model_selection import train_test_split

ROOT.gStyle.SetPaintTextFormat(".3f")

file      = uproot.open("../output-files/ntuple_e2c_unfolding.root")
# file_data = uproot.open("../output-files/ntuple_corre2c.root")
root_tree      = file["ntuple_unfolding"]
# root_tree_data = file_data["ntuple_hadron"]

# mc_set_rl    = np.array(root_tree["R_L_truth"])
# mc_set_jetpt = np.array(root_tree["jet_pt_truth"])
# recodata_set_rl    = np.array(root_tree_data["R_L"])
# recodata_set_jetpt = np.array(root_tree_data["jet_pt"])

# n_to_remove = np.size(recodata_set_rl) - np.size(mc_set_rl)
# for index in range(0,n_to_remove):
#     recodata_set_rl = np.delete(recodata_set_rl, index)
#     recodata_set_jetpt = np.delete(recodata_set_jetpt, index)

# Synthetic
mc_set        = np.column_stack((root_tree["R_L_truth"],root_tree["jet_pt_truth"]))
mcreco_set    = np.column_stack((root_tree["R_L"],root_tree["jet_pt"]))

# Data
# recodata_set = np.column_stack((recodata_set_rl,recodata_set_jetpt))

# Split sets to perform CT
# CT Steps : 
# - Divide the simulation set in two groups
# - Use one group to perform the Unfolding procedure.
#     - In this case you need to reweight not against data but to mc_set reco!
# - Then with the unfolded results of the first group ratio it to the mc_set of the second group

# mc_set_ofinput , mcreco_set_ofinput, mcreco_set_validate : Used to calculate the weights
# mc_set_validate : unfolded result should be compared to this set of values!

mc_set_ofinput , mc_set_validate, mcreco_set_ofinput, recodata_set_ofinput = train_test_split(mc_set, mcreco_set, random_state=42, test_size=0.5, train_size=0.5)

# New procedure for Multifold
sim_dataloader  = DataLoader(reco = mcreco_set_ofinput, gen = mc_set_ofinput)
data_dataloader = DataLoader(reco = recodata_set_ofinput)

nfeatures  = 2
reco_model = MLP(nfeatures)
gen_model  = MLP(nfeatures)

omnifold = MultiFold("EEC Unfolding", reco_model, gen_model, data = data_dataloader, mc = sim_dataloader, niter = 5)
omnifold.Unfold()

of_weights = omnifold.reweight(mc_set_validate, omnifold.model2, batch_size = 1000)

# Visualize
unfolding_rl_nbins = 17
R_L_min = 0.0099
R_L_max = 0.49

unfolding_rl_binning = array('d',[R_L_min-0.005,R_L_min, 0.0419067, 0.0739133, 0.10592, 0.137927, 0.169933,
                                          0.20194, 0.233947, 0.265953, 0.29796, 0.329967, 0.361973, 
                                          0.39398, 0.425987, 0.457993, R_L_max, R_L_max + 0.04])
unfolding_jetpt_binning = array('d',[15,20,30,50,100,150])

# Visualize
htruth_rl     = ROOT.TH1F("htruth_rl"    ,"",unfolding_rl_nbins,unfolding_rl_binning)
hunfol_rl     = ROOT.TH1F("hunfol_rl"    ,"",unfolding_rl_nbins,unfolding_rl_binning)
hct_rl        = ROOT.TH1F("hct_rl"       ,"",unfolding_rl_nbins,unfolding_rl_binning)
htruth_jetpt  = ROOT.TH1F("htruth_jetpt" ,"",5,unfolding_jetpt_binning)
hunfol_jetpt  = ROOT.TH1F("hunfol_jetpt" ,"",5,unfolding_jetpt_binning)
hct_jetpt     = ROOT.TH1F("hct_jetpt"    ,"",5,unfolding_jetpt_binning)

for entry in range(len(of_weights)-1):
    htruth_rl.Fill(mc_set_ofinput[entry,0])
    hunfol_rl.Fill(mc_set_validate[entry,0],of_weights[entry])
    
    htruth_jetpt.Fill(mc_set_ofinput[entry,1])
    hunfol_jetpt.Fill(mc_set_validate[entry,1],of_weights[entry])

hct_rl.Divide(htruth_rl,hunfol_rl,1,1)
hct_jetpt.Divide(htruth_jetpt,hunfol_jetpt,1,1)

fout = ROOT.TFile("../output-files/multifold_closuretest_2d_newof.root","RECREATE")
fout.cd()
hct_rl.Write()
hct_jetpt.Write()
fout.Close()

# # plt.savefig("./plots/closuretest-2d-multifold-{}iterations-{}nodes-{}epochs-{}nbatchsize.pdf".format(niterations, nnodes, nepochs, nbatch_size), format="pdf", bbox_inches="tight")
# # plt.show()
