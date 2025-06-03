import numpy as np
from matplotlib import pyplot as plt

from keras.layers import Dense, Input
from keras.models import Model

import omnifold as of
import os
import tensorflow as tf

# gpus = tf.config.experimental.list_physical_devices('GPU')
# print(gpus)

N = 145778

#Synthetic
theta0_G = np.random.normal(0.2,0.8,N)  # Generator-level synthetic sample
theta0_S = np.array([(x + np.random.normal(0, 0.5)) for x in theta0_G])  # Detector smearing for synthetic sample

theta0 = np.stack([theta0_G, theta0_S], axis=1)

#Natural
theta_unknown_G = np.random.normal(0,1, N)
theta_unknown_S = np.array([(x + np.random.normal(0, 0.5)) for x in theta_unknown_G]) 

inputs = Input((1, ))

# Calling layers of NN
hidden_layer_1 = Dense(50, activation='relu')(inputs) # relu = rectified linear unit activation function
hidden_layer_2 = Dense(50, activation='relu')(hidden_layer_1)
hidden_layer_3 = Dense(50, activation='relu')(hidden_layer_2)
hidden_layer_4 = Dense(50, activation='relu')(hidden_layer_3)
hidden_layer_5 = Dense(50, activation='relu')(hidden_layer_4)
outputs = Dense(1, activation='sigmoid')(hidden_layer_5)
model = Model(inputs=inputs, outputs=outputs)

# First use of Omnifold
myweights = of.omnifold(theta0,theta_unknown_S,2,model)

_,_,_=plt.hist(theta0_G,bins=np.linspace(-3,3,20),color='blue',histtype="step",label="MC, true",lw="2")
_,_,_=plt.hist(theta_unknown_G,bins=np.linspace(-3,3,20),color='orange',histtype="step",label="Data, true",lw="2")
_,_,_=plt.hist(theta0_G,weights=myweights[-1, 0, :], bins=np.linspace(-3,3,20),color='black',histtype="step",label="OmniFolded",lw="2")
plt.xlabel("x")
plt.ylabel("events")
plt.legend(frameon=False)

plt.show()