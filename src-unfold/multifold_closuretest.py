import numpy as np
from matplotlib import pyplot as plt

from keras.layers import Dense, Input
from keras.models import Model

import omnifold as of
import os
import tensorflow as tf

import uproot
from sklearn.model_selection import train_test_split

file      = uproot.open("../output-files/ntuple_e2c_unfolding.root")
file_data = uproot.open("../output-files/ntuple_corre2c.root")
root_tree      = file["ntuple_unfolding"]
root_tree_data = file_data["ntuple_data"]

# Regularize input 
theta0_G_rl_np = np.array(root_tree["R_L_truth"])
theta0_G_jp_np = np.array(root_tree["jet_pt_truth"])
theta0_G_w_np  = np.array(root_tree["weight_truth"])
theta_unknown_S_rl_np = np.array(root_tree_data["R_L"])
theta_unknown_S_jp_np = np.array(root_tree_data["jet_pt"])
theta_unknown_S_w_np  = np.array(root_tree_data["weight"])

n_to_remove = np.size(theta_unknown_S_rl_np) - np.size(theta0_G_rl_np)
for index in range(0,n_to_remove):
    theta_unknown_S_rl_np = np.delete(theta_unknown_S_rl_np, index)
    theta_unknown_S_jp_np = np.delete(theta_unknown_S_jp_np, index)
    theta_unknown_S_w_np  = np.delete(theta_unknown_S_w_np, index)

# Synthetic
theta0_G = np.column_stack((theta0_G_rl_np,theta0_G_jp_np,theta0_G_w_np))
theta0_S = np.column_stack((root_tree["R_L"],root_tree["jet_pt"],root_tree["weight"]))

# Data
theta_unknown_S = np.column_stack((theta_unknown_S_rl_np,theta_unknown_S_jp_np,theta_unknown_S_w_np))

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
niterations = 6

inputs = Input((np.shape(theta0_G_train)[-1], ))
hidden_layer_1 = Dense(nnodes, activation='relu')(inputs) # relu = rectified linear unit activation function
hidden_layer_2 = Dense(nnodes, activation='relu')(hidden_layer_1)
hidden_layer_3 = Dense(nnodes, activation='relu')(hidden_layer_2)
outputs = Dense(1, activation='sigmoid')(hidden_layer_3) # Must have one output!
model = Model(inputs=inputs, outputs=outputs)

myweights = of.omnifold(theta0_G_train, theta0_S_train, theta_unknown_S_train, niterations, model, 0,nepochs, nbatch_size)

# Visualize
figtest, ((ax11test),(ax22test),(ax33test)) = plt.subplots(1, 3, layout='tight')
fig, ((ax11),(ax22),(ax33)) = plt.subplots(1, 3, layout='tight', sharey=True)
fig.set_figwidth(15)
fig.set_figheight(5)

h_rl_validate    ,bins_rl_validate    ,_ = ax11test.hist(theta0_G_validate[:,0],bins=np.linspace(0.00999,0.5,15),color='orange',histtype="step",label="Data, reco",lw="2")
h_jetpt_validate ,bins_jetpt_validate ,_ = ax22test.hist(theta0_G_validate[:,1],bins=[20,30,50,100]     ,color='orange',histtype="step",label="Data, reco",lw="2")
h_weight_validate,bins_weight_validate,_ = ax33test.hist(theta0_G_validate[:,2],bins=np.geomspace(0.00001,.01,10)    ,color='orange',histtype="step",label="Data, reco",lw="2")

h_rl_unfolded    ,bins_rl_unfolded    ,_ = ax11test.hist(theta0_G_train[:,0], weights=myweights[niterations-1, 1, :], bins=np.linspace(0.00999,0.5,15),color='black',histtype="step",label="Unfolded",lw="2")
h_jetpt_unfolded ,bins_jetpt_unfolded ,_ = ax22test.hist(theta0_G_train[:,1], weights=myweights[niterations-1, 1, :], bins=[20,30,50,100]     ,color='black',histtype="step",label="Unfolded",lw="2")
h_weight_unfolded,bins_weight_unfolded,_ = ax33test.hist(theta0_G_train[:,2], weights=myweights[niterations-1, 1, :], bins=np.geomspace(0.00001,.01,10)    ,color='black',histtype="step",label="Unfolded",lw="2")

ct_rl     = np.array(h_rl_unfolded/h_rl_validate)
ct_jetpt  = np.array(h_jetpt_unfolded/h_jetpt_validate)
ct_weight = np.array(h_weight_unfolded/h_weight_validate)

ratio_rl_xaxis, ratio_jetpt_xaxis, ratio_weight_xaxis = [], [], []
for i in range(0,np.size(h_rl_validate)):
    delta_rl = (bins_rl_validate[i+1]-bins_rl_validate[i])
    ratio_rl_xaxis.append(bins_rl_validate[i] + delta_rl/2.)
ratio_rl_xaxis = np.array(ratio_rl_xaxis)

for i in range(0,np.size(h_jetpt_validate)):
    delta_jetpt = (bins_jetpt_validate[i+1]-bins_jetpt_validate[i])
    ratio_jetpt_xaxis.append(bins_jetpt_validate[i] + delta_jetpt/2.)
ratio_jetpt_xaxis = np.array(ratio_jetpt_xaxis)

for i in range(0,np.size(h_weight_validate)):
    delta_weight = (bins_weight_validate[i+1]-bins_weight_validate[i])
    ratio_weight_xaxis.append(bins_weight_validate[i] + delta_weight/2.)
ratio_weight_xaxis = np.array(ratio_weight_xaxis)

ax11.set(xlabel="$R_L$",ylabel="Unfolded/Truth",xscale='log')
ax11.scatter(ratio_rl_xaxis,ct_rl,marker='o')

ax22.set(xlabel="$p^{jet}_T$",ylabel="Unfolded/Truth")
ax22.scatter(ratio_jetpt_xaxis,ct_jetpt,marker='o')

ax33.set(xlabel="$weight$",ylabel="Unfolded/Truth", xscale='log')
ax33.scatter(ratio_weight_xaxis,ct_weight,marker='o')

ax11.grid(axis='y', color='gainsboro', linewidth=0.5)
ax11.set_ylim((0.9,1.1))
ax22.grid(axis='y', color='gainsboro', linewidth=0.5)
ax33.grid(axis='y', color='gainsboro', linewidth=0.5)

plt.suptitle('Iter {}, 3 layers , #Nodes = {} , #Epochs = {} , Batch size = {}'.format(niterations,nnodes,nepochs,nbatch_size))
plt.savefig("./plots/closuretest-3d-multifold-{}iterations-{}nodes-{}epochs-{}nbatchsize.pdf".format(niterations, nnodes, nepochs, nbatch_size), format="pdf", bbox_inches="tight")
# plt.show()
