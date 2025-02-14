import numpy as np
from matplotlib import pyplot as plt

from keras.layers import Dense, Input
from keras.models import Model

import omnifold as of
import os
import tensorflow as tf

import uproot

file      = uproot.open("../output-files/ntuple_e2c_unfolding.root")
file_data = uproot.open("../output-files/ntuple_corre2c.root")
root_tree      = file["ntuple_unfolding"]
root_tree_data = file_data["ntuple_data"]

# Regularize input 
theta0_G_rl_np = np.array(root_tree["R_L_truth"]) # Use as reference
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

# Calling layers of NN
nnodes      = 7
niterations = 5

inputs = Input((np.shape(theta0_G)[-1], ))
hidden_layer_1 = Dense(nnodes, activation='relu')(inputs) # relu = rectified linear unit activation function
hidden_layer_2 = Dense(nnodes, activation='relu')(hidden_layer_1)
hidden_layer_3 = Dense(nnodes, activation='relu')(hidden_layer_2)
outputs = Dense(1, activation='sigmoid')(hidden_layer_3) # Must have one output!
model = Model(inputs=inputs, outputs=outputs)

myweights = of.omnifold(theta0_G, theta0_S,theta_unknown_S,niterations,model)

# Visualize
fig, ((ax11),(ax22)) = plt.subplots(1, 2, layout='tight')
fig.set_figwidth(12)
fig.set_figheight(6)

_,_,_ = ax11.hist(theta0_G[:,0]       ,bins=np.linspace(0.00999,0.5,15),color='blue',histtype="step",label="MC, true",lw="2")
_,_,_ = ax11.hist(theta0_S[:,0]       ,bins=np.linspace(0.00999,0.5,15),color='green',histtype="step",label="MC, reco",lw="2")
_,_,_ = ax11.hist(theta_unknown_S[:,0],bins=np.linspace(0.00999,0.5,15),color='orange',histtype="step",label="Data, reco",lw="2")
_,_,_ = ax11.hist(theta0_G[:,0]       ,weights=myweights[niterations-1, 1, :], bins=np.linspace(0.00999,0.5,15),color='black',histtype="step",label="Unfolded",lw="2")

_,_,_ = ax22.hist(theta0_G[:,1]       ,bins=np.linspace(15,100,25),color='blue',histtype="step",label="MC, true",lw="2")
_,_,_ = ax22.hist(theta0_S[:,1]       ,bins=np.linspace(15,100,25),color='green',histtype="step",label="MC, reco",lw="2")
_,_,_ = ax22.hist(theta_unknown_S[:,1],bins=np.linspace(15,100,25),color='orange',histtype="step",label="Data, reco",lw="2")
_,_,_ = ax22.hist(theta0_G[:,1]       ,weights=myweights[niterations-1, 1, :], bins=np.linspace(15,100,25),color='black',histtype="step",label="Unfolded",lw="2")

# ax11.set(xlabel="$R_{L}$", ylabel="$N_{pair}$", xscale="log", yscale="log")
# ax22.set(xlabel="$p^{jet}_{T}$", ylabel="$N_{pair}$", yscale="log")

ax11.set(xlabel="$R_{L}$", ylabel="$N_{pair}$", xscale="log")
ax22.set(xlabel="$p^{jet}_{T}$", ylabel="$N_{pair}$")

ax11.legend(frameon=False,loc='upper left')
ax22.legend(frameon=False,loc='upper right')

plt.savefig("./plots/{}-iterations_{}-hlayers_{}-nodes_rl_jetpt.pdf".format(niterations, 3, nnodes), format="pdf", bbox_inches="tight")
plt.show()
# for ax in (ax21, ax22):
#     for i in range(2):

# # histo returns array , bins , patches
# var_id = 0
# h1,h1bins,_ = ax1.hist(theta0_G[:,var_id]       ,bins=np.linspace(0.00999,0.5,15),color='blue',histtype="step",label="MC, true",lw="2")
# h2,h2bins,_ = ax1.hist(theta0_S[:,var_id]       ,bins=np.linspace(0.00999,0.5,15),color='green',histtype="step",label="MC, reco",lw="2")
# h3,h3bins,_ = ax1.hist(theta_unknown_S[:,var_id],bins=np.linspace(0.00999,0.5,15),color='orange',histtype="step",label="Data, reco",lw="2")
# h4,h4bins,_ = ax1.hist(theta0_G[:,var_id]       ,weights=myweights[niterations-1, 1, :], bins=np.linspace(0.00999,0.5,15),color='black',histtype="step",label="Unfolded",lw="2")

# ax1.set(ylabel="$N_{pair}$", xscale="log", yscale="log")

# ratio_y=np.array(h3/h4)

# ratio_x_list = []
# for i in range(0,np.size(h1bins)-1):
#     delta_rl = (h1bins[i+1]-h1bins[i])
#     ratio_x_list.append(h1bins[0] + (i+1/2)*delta_rl)
# ratio_x = np.array(ratio_x_list)

# ax2.set(xlabel="$R_L$",ylabel="Data reco / OF")
# ax2.set_ylim((0.95,1.05))
# ax2.plot(ratio_x,ratio_y,marker='o')

# ax1.legend(frameon=False,loc='upper left')
