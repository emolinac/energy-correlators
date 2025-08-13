import numpy as np
import omnifold as of
import os
import tensorflow as tf
import uproot
import ROOT

from matplotlib import pyplot as plt
from keras.layers import Dense, Input
from keras.models import Model
from array import array
from sklearn.model_selection import train_test_split

ROOT.gStyle.SetPaintTextFormat(".3f")

file      = uproot.open("../output-files/ntuple_eec_unfolding.root")
file_data = uproot.open("../output-files/ntuple_correec.root")
root_tree      = file["ntuple_unfolding"]
root_tree_data = file_data["ntuple_data"]

# Regularize input 
theta0_G_rl_np = np.array(root_tree["R_L_truth"])
theta0_G_jp_np = np.array(root_tree["jet_pt_truth"])
theta_unknown_S_rl_np = np.array(root_tree_data["R_L"])
theta_unknown_S_jp_np = np.array(root_tree_data["jet_pt"])

n_to_remove = np.size(theta_unknown_S_rl_np) - np.size(theta0_G_rl_np)
for index in range(0,n_to_remove):
    theta_unknown_S_rl_np = np.delete(theta_unknown_S_rl_np, index)
    theta_unknown_S_jp_np = np.delete(theta_unknown_S_jp_np, index)

# Synthetic
theta0_G = np.column_stack((theta0_G_rl_np,theta0_G_jp_np))
theta0_S = np.column_stack((root_tree["R_L"],root_tree["jet_pt"]))

# Data
theta_unknown_S = np.column_stack((theta_unknown_S_rl_np,theta_unknown_S_jp_np))

# Split sets to perform CT
# CT Steps : 
# - Divide the simulation set in two groups
# - Use one group to perform the Unfolding procedure.
#     - In this case you need to reweight not against data but to MC reco!
# - Then with the unfolded results of the first group ratio it to the MC of the second group

# theta0_G_train , theta0_S_train, theta0_S_validate : Used to calculate the weights
# theta0_G_validate : unfolded result should be compared to this set of values!

theta0_G_train , theta0_G_validate, theta0_S_train, theta_unknown_S_train = train_test_split(theta0_G, theta0_S, random_state=42, test_size=0.5, train_size=0.5)

# Calling layers of NN
nnodes      = 10
nepochs     = 25
nbatch_size = 5000
niterations = 20

inputs = Input((np.shape(theta0_G_train)[-1], ))
hidden_layer_1 = Dense(nnodes, activation='relu')(inputs) # relu = rectified linear unit activation function
hidden_layer_2 = Dense(nnodes, activation='relu')(hidden_layer_1)
hidden_layer_3 = Dense(nnodes, activation='relu')(hidden_layer_2)
outputs = Dense(1, activation='sigmoid')(hidden_layer_3) # Must have one output!
model = Model(inputs=inputs, outputs=outputs)

myweights = of.omnifold(theta0_G_train, theta0_S_train, theta_unknown_S_train, niterations, model, 0,nepochs, nbatch_size)

# Visualize
# rl_nbins = 16
unfolding_rl_nbins = 17

# rl_binning     = array('d',np.linspace(0.00999,0.5,rl_nbins))
# jet_pt_binnning = array('d',[20,30,50,100])

rl_min = 0.0099
rl_max = 0.49

unfolding_rl_binning = array('d',[rl_min-0.005,rl_min, 0.0419067, 0.0739133, 0.10592, 0.137927, 0.169933,
                                          0.20194, 0.233947, 0.265953, 0.29796, 0.329967, 0.361973, 
                                          0.39398, 0.425987, 0.457993, rl_max, rl_max + 0.04])
unfolding_jet_pt_binning = array('d',[15,20,30,50,100,150])

# Visualize
htruth_rl     = ROOT.TH1F("htruth_rl"    ,"",unfolding_rl_nbins,unfolding_rl_binning)
hunfol_rl     = ROOT.TH1F("hunfol_rl"    ,"",unfolding_rl_nbins,unfolding_rl_binning)
hct_rl        = ROOT.TH1F("hct_rl"       ,"",unfolding_rl_nbins,unfolding_rl_binning)
htruth_jet_pt  = ROOT.TH1F("htruth_jet_pt" ,"",5,unfolding_jet_pt_binning)
hunfol_jet_pt  = ROOT.TH1F("hunfol_jet_pt" ,"",5,unfolding_jet_pt_binning)
hct_jet_pt     = ROOT.TH1F("hct_jet_pt"    ,"",5,unfolding_jet_pt_binning)

for entry in range(len(myweights[niterations-1, 1, :])-1):
    htruth_rl.Fill(theta0_G_validate[entry,0])
    hunfol_rl.Fill(theta0_G_train[entry,0],myweights[niterations-1, 1,entry])
    
    htruth_jet_pt.Fill(theta0_G_validate[entry,1])
    hunfol_jet_pt.Fill(theta0_G_train[entry,1],myweights[niterations-1, 1,entry])

hct_rl.Divide(htruth_rl,hunfol_rl,1,1)
hct_jet_pt.Divide(htruth_jet_pt,hunfol_jet_pt,1,1)

fout = ROOT.TFile("../output-files/multifold_closuretest_2d.root","RECREATE")
fout.cd()
hct_rl.Write()
hct_jet_pt.Write()
fout.Close()

# plt.savefig("./plots/closuretest-2d-multifold-{}iterations-{}nodes-{}epochs-{}nbatchsize.pdf".format(niterations, nnodes, nepochs, nbatch_size), format="pdf", bbox_inches="tight")
# plt.show()
