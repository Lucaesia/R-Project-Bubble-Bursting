import numpy as np
from matplotlib import pyplot as plt
from matplotlib import collections  as mc


# sets time inverval between files
dt = 0.1


#variable for tracking number of files
total_points = 0
while total_points==0:
    # finds file name and exits loop once finished

    file_name = "../DriverCode/Water-R0.0071-Level8/Interfaces/interfaceDrop-3.5.dat"
    try:
        data = np.loadtxt(file_name, delimiter=' ', usecols=(0,1), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
   
    total_points+=1

# swaps collums and then reagranges into line segments
data[:, [0, 1]] = data[:, [1, 0]]
data_mirrored = np.concatenate((data,data*[-1,1]))

lines = data_mirrored.reshape((len(data_mirrored)//2,2,2))


# create a LineCollection with the half-circles
line_collection = mc.LineCollection(lines)

# plots and saves graphs
fig, ax = plt.subplots(figsize=(6.4, 3.2))
# set axes limits manually because Collections do not take part in autoscaling
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_aspect("equal")  # to make the arcs look circular



ax.add_collection(line_collection)

plt.savefig('interface.png')
plt.clf() 
