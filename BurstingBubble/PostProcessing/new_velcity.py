import numpy as np # type: ignore
from matplotlib import pyplot as plt # type: ignore
from matplotlib import collections  as mc # type: ignore


# sets time inverval between files
dt = 0.1

scale=1
fig, ax = plt.subplots(6,figsize=(6.4*scale, 6.4*scale))
#variable for tracking number of files
total_points = 0
for i in [6,7,8,9,10,11]:
    # finds file name and exits loop once finished
    R = 0.0071
    Level = 8
    t = 10.0
    
    file_name = "../DriverCode/Water-R0.0031-Level"+str(i)+"/jet_vel.dat"

    
    try:
        jet_vel = np.loadtxt(file_name, delimiter=' ', usecols=(0,1,2,3), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        
   
    #ax[i-6].plot(jet_vel[:,3],jet_vel[:,2])
    ax[i-6].plot(jet_vel[:,3],jet_vel[:,1], label="Level = "+str(i))
    ax[i-6].set_ylim([-1,10])
    ax[i-6].set_xlim([0,3])


### JET VEL ###


# swaps collums and then reagranges into line segments


# plots and saves graphs

# set axes limits manually because Collections do not take part in autoscaling

#ax.set_aspect("equal")  # to make the arcs look circular




#plt.legend()

plt.savefig('bubble_velocity.png')
plt.clf()


