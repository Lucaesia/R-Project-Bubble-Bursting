import numpy as np
from matplotlib import pyplot as plt



file_name = "stability/time=2.900.dat"
file_name_2 = "stability/time=3.000.dat"
try:
    x,u_1= np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=True)
    x,u_2= np.loadtxt(file_name_2, delimiter=' ', usecols=(0,1), unpack=True)
    #nx,ny,nz,nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3,nx4,ny4,nz4 = np.loadtxt("x-y-plane-omp.dat", delimiter=' ', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), unpack=True)
#generated_plots/x-y-plane_large_interaction
except FileNotFoundError:
    print("file not found"+file_name)





# plots and saves graphs
fig = plt.figure()
ax = fig.add_subplot()
ax.set_title("Minimum radius over time - dt=0.001")
ax.set_xlabel("time")
ax.set_ylabel("radius")

ax.plot(x,u_1-u_2,linestyle='-',linewidth=2)

#ax.scatter(0,0,100,marker='.',label="Black Hole")

#plt.legend(loc ="lower right")
#ax.set_aspect('equal', adjustable='box')
print(np.max(np.absolute(u_1-u_2)))
plt.savefig('difference.png')
plt.clf() 

