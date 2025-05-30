import numpy as np # type: ignore
from matplotlib import pyplot as plt # type: ignore
from matplotlib import collections  as mc # type: ignore


# sets time inverval between files
dt = 0.1


#variable for tracking number of files
total_points = 0
while total_points==0:
    # finds file name and exits loop once finished
    R = 0.0071
    Level = 8
    t = 10.0
    
    file_name = "../DriverCode/Water-R0.0031-Level8/jet_vel.dat"
    try:
        jet_vel = np.loadtxt(file_name, delimiter=' ', usecols=(0,1,2,3), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
   
    total_points+=1

### JET VEL ###


# swaps collums and then reagranges into line segments


# plots and saves graphs
scale=1
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))
# set axes limits manually because Collections do not take part in autoscaling

#ax.set_aspect("equal")  # to make the arcs look circular
ax.set_ylim([-1,5])



ax.plot(jet_vel[:,3],jet_vel[:,2])
ax.plot(jet_vel[:,3],jet_vel[:,1], label="interpolated vel")

plt.legend()

plt.savefig('bubble_vel.png')
plt.clf()

scale=1
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))
# set axes limits manually because Collections do not take part in autoscaling

#ax.set_aspect("equal")  # to make the arcs look circular
ax.set_ylim([-4,5])
ax.plot(jet_vel[:,3],jet_vel[:,0])
plt.savefig('bubble_vel1.png')
plt.clf() 
