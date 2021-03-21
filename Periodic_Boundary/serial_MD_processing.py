import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

## Importing data
imported_data = pd.read_csv("molecular_data.csv",sep = ',')
parameters = pd.read_csv("parameters.csv", sep = ',')

N = parameters.N.to_numpy()[0]
nsteps = parameters.iters.to_numpy()[0]

## Reshaping to get a matrix where each column is a particle and each row a time step.
#N is the number of particles.
#nsteps is the number of iterations.
x = np.reshape(imported_data.x.to_numpy(),(nsteps,N))
y = np.reshape(imported_data.y.to_numpy(),(nsteps,N))
vx = np.reshape(imported_data.vx.to_numpy(),(nsteps,N))
vy = np.reshape(imported_data.vy.to_numpy(),(nsteps,N))
K = np.sum(np.reshape(imported_data.K.to_numpy(),(nsteps,N)),axis =1)
W = np.sum(np.reshape(imported_data.W.to_numpy(),(nsteps,N)), axis = 1)


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



from matplotlib import pyplot as plt
from celluloid import Camera
import numpy as np


# create figure object
fig = plt.figure()
# load axis box
ax = plt.axes()
# set axis limit
#ax.set_ylim(0, 10)
#ax.set_xlim(0, 10)
L = np.max(np.max(x))
camera = Camera(fig)
for i in range(100):
    ax.set(xlim=(0, L), ylim=(0, L))
    ax.set_title("N = "+str(N)+" particles")
    ax.scatter(x[1000*i,:],y[1000*i,:])
    plt.pause(0.01)
    camera.snap()

    plt.gca().clear()
#    ax.clear()

animation = camera.animate()
#animation.save('animation.gif', writer='PillowWriter', fps=2)
animation.save('animation.gif') #, writer='imagemagick', fps=2)
animation.save('animation.gif' , writer='matplotlib.animation.PillowWriter', fps=2)







#sigma = 1
#epsilon = 1
#
#rc = 3*sigma
#delta = 0.1
#
#r = np.linspace(0.1,10*sigma,1000)
#
#A = np.array([[1,rc,rc**2,rc**3],[0,1,2*rc,3*rc**2],[1,rc+delta,np.power(rc+delta,2),np.power(rc+delta,3)],[0,1,2*(rc+delta),3*np.power(rc+delta,2)]])
#b = np.array([4*epsilon*(np.power(sigma/rc,12)-np.power(sigma/rc,6)),24*epsilon*((np.power(sigma,6)/np.power(rc,7))-(2*np.power(sigma,12)/np.power(rc,13))) , 0, 0])
#
#a,b,c,d = np.linalg.solve(A,b)
#
#LJ = 4*epsilon*(np.power(sigma/r,12)-np.power(sigma/r,6))
#polynomial =  a + b*r + c*r**2+ d*r**3
#
#fig,ax = plt.subplots(1,1)
#ax.plot(r,LJ, label = 'Lennard-Jones Potential')
#ax.plot(r,polynomial, label = 'Polynomial cut-off')
#ax.set_xlabel('r')
#ax.set_ylabel('potential')
#ax.axvline(x=rc,color = 'b')
#ax.axvline(x=rc+delta, color = 'r')
#ax.axhline(y=0, color ='k')
#ax.set_ylim(-1,1)
#ax.legend()
#ax.set_title('Illustration of smooth cut-off for the L-J potential')
#
#
#rc = 3.05#np.sqrt(8);
#
#pot_deriv1 =  24*epsilon*((np.power(sigma,6)/np.power(rc,7))-(2*np.power(sigma,12)/np.power(rc,13)))
#pot_deriv2 = b+2*c*rc+3*d*np.power(rc,2)
#
#fx = (pot_deriv2/rc)*(-3.05)
#fy = (pot_deriv2/rc)*(0)
