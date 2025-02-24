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
    
    file_name = "../DriverCode/Water-R0.0031-Level8/logdroplets.dat"
    try:
        droplet_dat = np.loadtxt(file_name, delimiter=' ', usecols=(0,1,2,3,4,5), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
   
    total_points+=1

print(np.unique(droplet_dat[:,2]))
drop_0 = droplet_dat [droplet_dat[:,2]==0]
drop_1 = droplet_dat [droplet_dat[:,2]==1]
drop_2 = droplet_dat [droplet_dat[:,2]==2]
drop_3 = droplet_dat [droplet_dat[:,2]==3]


# plots and saves graphs
scale=1
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))
# set axes limits manually because Collections do not take part in autoscaling

#ax.set_aspect("equal")  # to make the arcs look circular
ax.set_ylim([-1,5])



ax.plot(drop_1[:,1],drop_1[:,4],label="Drop 1")
ax.plot(drop_2[:,1],drop_2[:,4],label="Drop 2")
ax.plot(drop_3[:,1],drop_3[:,4],label="Drop 3")

plt.legend()

plt.savefig('outputs/bubble properties.png')
plt.clf()


