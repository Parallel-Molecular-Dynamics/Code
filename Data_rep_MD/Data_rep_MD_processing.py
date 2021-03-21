import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## Importing data
molecular_data = pd.read_csv("molecular_data.csv",sep = ',')
parameters = pd.read_csv("parameters.csv", sep = ',')
energy_data = pd.read_csv("energy_data.csv",sep = ',')

N = parameters.N.to_numpy()[0]
nsteps = parameters.iters.to_numpy()[0]

## Reshaping to get a matrix where each column is a particle and each row a time step.
x = np.reshape(molecular_data.x.to_numpy(),(nsteps,N))
y = np.reshape(molecular_data.y.to_numpy(),(nsteps,N))
vx = np.reshape(molecular_data.vx.to_numpy(),(nsteps,N))
vy = np.reshape(molecular_data.vy.to_numpy(),(nsteps,N))
K = energy_data.K.to_numpy()
W = energy_data.W.to_numpy()


# Snapshot of particles at a given timestep            
fig,ax = plt.subplots(1,1)

ax.scatter(x[nsteps-1,:],y[nsteps-1,:], marker = 'x')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Particle positions at timestep=' + str(nsteps-1))

# Energy Plot
fig,ax = plt.subplots(1,1)

ax.plot(K+W)
ax.set_xlabel('Time step')
ax.set_ylabel('E')
ax.set_title('Energy')

if(N==2):
## Plots for a single pair of particles.
    fig,ax = plt.subplots(1,1)
    
    ax.plot(y[:,0], label = 'Particle 1')
    ax.plot(y[:,1], label = 'Particle 2')
    ax.set_xlabel('Time step')
    ax.set_ylabel('y')
    ax.set_title('y positions')
    ax.legend()
    