import numpy as np
from matplotlib import pyplot as plt
import os
import uproot

file      = uproot.open("../output-files/ntuple_eec_unfolding.root")
root_tree = file["ntuple_unfolding"]

# Regularize input 
theta0_G_rl_np = np.array(root_tree["R_L_truth"])
theta0_G_jp_np = np.array(root_tree["jet_pt_truth"])
theta0_G_w_np  = np.array(root_tree["weight_truth"])
theta0_S_rl_np = np.array(root_tree["R_L"])
theta0_S_jp_np = np.array(root_tree["jet_pt"])
theta0_S_w_np  = np.array(root_tree["weight"])

# Visualize
fig, ((ax11),(ax22),(ax33)) = plt.subplots(1, 3, layout='tight')
fig.set_figwidth(15)
fig.set_figheight(5)

hist_rl , xbins_rl , ybins_rl , _ = ax11.hist2d(theta0_G_rl_np ,theta0_S_rl_np, bins=[np.linspace(0.01,0.5,15),np.linspace(0.01,0.5,15)])#bins=[np.linspace(0.00999,0.5,15),np.linspace(0.00999,0.5,15)])
hist_jp , xbins_jp , ybins_jp , _ = ax22.hist2d(theta0_G_jp_np ,theta0_S_jp_np, bins=[[20,30,50,100],[20,30,50,100]])
hist_w , xbins_w  , ybins_w  , _ = ax33.hist2d(theta0_G_w_np  ,theta0_S_w_np , bins=[np.geomspace(0.000001,1,10),np.geomspace(0.000001,1,10)]          )

ax11.set(xlabel="Truth",ylabel="Reco",title="$R_{L}$")
ax22.set(xlabel="Truth",ylabel="Reco",title="$p_{T,jet}$")
ax33.set(xlabel="Truth",ylabel="Reco",title="$weight$",xscale='log',yscale='log')

# plt.suptitle('Iter {}, 3 layers , #Nodes = {} , #Epochs = {} , Batch size = {}'.format(niterations,nnodes,nepochs,nbatch_size))
plt.show()
