import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

imported_data = pd.read_csv("molecular_data.csv",sep = ' ')
parameters = pd.read_csv("parameters.csv", sep = ' ')

N = parameters.N.to_numpy()[0]
nsteps = parameters.iters.to_numpy()[0]

x = imported_data.x.to_numpy()
y = imported_data.y.to_numpy()

x_positions = np.reshape(x,(nsteps,N))
y_positions = np.reshape(y,(nsteps,N))

plt.plot(y_positions[:,0])
plt.plot(y_positions[:,1])

y_separation = abs(y_positions[:,1]-y_positions[:,0])


















































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
