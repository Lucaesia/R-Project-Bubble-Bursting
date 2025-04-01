import numpy as np # type: ignore
from matplotlib import pyplot as plt # type: ignore
from matplotlib import collections  as mc # type: ignore

# finds file name and exits loop once finished
R = "0.0011"
IC = "IC"
Level = 10
t = 10.0
#string = "%.2f" % float(total_points*0.01)
file_name = "../DriverCode/Water-"+IC+"-R"+R+"-Level"+str(Level)+"/interfacestats.dat"
try:
    data = np.loadtxt(file_name, delimiter=' ', usecols=(0,1,2), unpack=False)
except FileNotFoundError:
    print("file not found"+file_name)
    

scale=4
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))
# set axes limits manually because Collections do not take part in autoscaling
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_aspect("equal")  # to make the arcs look circular


print(data)
ax.scatter(data[:,0],data[:,1],c=data[:,2])


plt.savefig('initial interface.png')
plt.clf() 
plt.close()
