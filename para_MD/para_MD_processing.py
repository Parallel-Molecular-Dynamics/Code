import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from celluloid import Camera

## Importing data
imported_data = pd.read_csv("molecular_data.csv",sep = ',')
parameters = pd.read_csv("parameters.csv", sep = ',')
energy_data = pd.read_csv("energy_data.csv",sep = ',')

N = parameters.N.to_numpy()[0]
nsteps = parameters.iters.to_numpy()[0]
L  = parameters.L.to_numpy()[0]

## Reshaping to get a matrix where each column is a particle and each row a time step.
#N is the number of particles.
#nsteps is the number of iterations.
x = np.reshape(imported_data.x.to_numpy(),(nsteps,N))
y = np.reshape(imported_data.y.to_numpy(),(nsteps,N))
vx = np.reshape(imported_data.vx.to_numpy(),(nsteps,N))
vy = np.reshape(imported_data.vy.to_numpy(),(nsteps,N))

K = energy_data.K.to_numpy()
W = energy_data.W.to_numpy()


# Snapshot of particles at a given timestep
fig,ax = plt.subplots(1,1)

ax.scatter(x[nsteps-1,:],y[nsteps-1,:], marker = 'x')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Particle positions at timestep=' + str(nsteps-1))

## Energy Plot
fig,ax = plt.subplots(1,1)

ax.plot(K+W)
ax.set_xlabel('Time step')
ax.set_ylabel('E')
ax.set_title('Energy')


## Plots for a single pair of particles.
fig,ax = plt.subplots(1,1)

ax.plot(y[:,0], label = 'Particle 1')
ax.plot(y[:,1], label = 'Particle 2')
ax.set_xlabel('Time step')
ax.set_ylabel('y')
ax.set_title('y positions')
ax.legend()

fig,ax = plt.subplots(1,1)

ax.plot(x[:,0], label = 'Particle 1')
ax.plot(x[:,1], label = 'Particle 2')
ax.set_xlabel('Time step')
ax.set_ylabel('x')
ax.set_title('x positions')
ax.legend()

## Plots for a single pair of particles.
fig,ax = plt.subplots(1,1)

ax.plot(vy[:,0], label = 'Particle 1')
ax.plot(vy[:,1], label = 'Particle 2')
ax.set_xlabel('Time step')
ax.set_ylabel('y')
ax.set_title('y velocity')
ax.legend()

fig,ax = plt.subplots(1,1)

ax.plot(vx[:,0], label = 'Particle 1')
ax.plot(vx[:,1], label = 'Particle 2')
ax.set_xlabel('Time step')
ax.set_ylabel('x')
ax.set_title('x velocity')
ax.legend()

plt.show()

# create figure object
fig = plt.figure()
# load axis box
ax = plt.axes()
# set axis limit
#ax.set_ylim(0, 10)
#ax.set_xlim(0, 10)
camera = Camera(fig)
leap = 200
for i in range(int(nsteps/leap)):
    ax.set_xlim(0,L)
    ax.set_ylim(0,L)
    ax.set_title("N = "+str(N)+" particles")
    ax.scatter(x[leap*i,:],y[leap*i,:])
    ax.grid()
    plt.pause(0.01)
    camera.snap()

    plt.gca().clear()
#    ax.clear()

animation = camera.animate()
#animation.save('animation.gif', writer='PillowWriter', fps=2)
animation.save('animation.gif', writer='imagemagick', fps=2)
#animation.save('animation.gif' , writer='matplotlib.animation.PillowWriter', fps=2)
