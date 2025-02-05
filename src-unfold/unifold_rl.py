import numpy as np
from matplotlib import pyplot as plt

from keras.layers import Dense, Input
from keras.models import Model

import omnifold as of
import os
import tensorflow as tf

import uproot

# Specs
niterations = 2
nnodes = 6

# Get root files
file      = uproot.open("../output-files/ntuple_e2c_unfolding.root")
file_data = uproot.open("../output-files/ntuple_corre2c.root")
root_tree      = file["ntuple_unfolding"]
root_tree_data = file_data["ntuple_data"]

# Synthetic
theta0_G = np.array(root_tree["R_L_truth"])  # Generator-level synthetic sample
theta0_S = np.array(root_tree["R_L"])  # Detector smearing for synthetic sample

# Data
# theta_unknown_G = np.random.normal(0,1, np.size(theta0_S))
theta_unknown_S = np.array(root_tree_data["R_L"])

# Regularize size of samples
n_to_remove = np.size(theta_unknown_S) - np.size(theta0_S)
for index in range(0,n_to_remove):
    theta_unknown_S = np.delete(theta_unknown_S, index)

# Calling layers of NN
inputs = Input((1, ))
hidden_layer_1 = Dense(nnodes, activation='relu')(inputs) # relu = rectified linear unit activation function
hidden_layer_2 = Dense(nnodes, activation='relu')(hidden_layer_1)
hidden_layer_3 = Dense(nnodes, activation='relu')(hidden_layer_2)
outputs = Dense(1, activation='sigmoid')(hidden_layer_3)
model = Model(inputs=inputs, outputs=outputs)

# First use of Omnifold
myweights = of.omnifold(theta0_G, theta0_S,theta_unknown_S,niterations,model)

# Visualize
fig, (ax1, ax2) = plt.subplots(2, sharex = True)

# histo returns array , bins , patches
h1,h1bins,_ = ax1.hist(theta0_G,bins=np.linspace(0.00999,0.5,15),color='blue',histtype="step",label="MC, true",lw="2")
h2,h2bins,_ = ax1.hist(theta0_S,bins=np.linspace(0.00999,0.5,15),color='green',histtype="step",label="MC, reco",lw="2")
h3,h3bins,_ = ax1.hist(theta_unknown_S,bins=np.linspace(0.00999,0.5,15),color='orange',histtype="step",label="Data, reco",lw="2")
h4,h4bins,_ = ax1.hist(theta0_G,weights=myweights[-1, 0, :], bins=np.linspace(0.00999,0.5,15),color='black',histtype="step",label="OmniFolded",lw="2")

ax1.set(ylabel="$N_{pair}$", xscale="log", yscale="log")

ratio_y=np.array(h3/h4)

ratio_x_list = []
for i in range(0,np.size(h1bins)-1):
    delta_rl = (h1bins[i+1]-h1bins[i])
    ratio_x_list.append(h1bins[0] + (i+1/2)*delta_rl)
ratio_x = np.array(ratio_x_list)

ax2.set(xlabel="$R_L$",ylabel="Data reco / OF")
ax2.set_ylim((0.95,1.05))
ax2.plot(ratio_x,ratio_y,marker='o',label="Iter {}".format(niterations))

ax1.legend(frameon=False,loc='upper left')
ax2.legend(frameon=False,loc='upper left')
# plt.savefig("./plots/{}-iterations_{}-hlayers_{}-nodes_of.pdf".format(niterations, 3, nnodes), format="pdf", bbox_inches="tight")
plt.show()