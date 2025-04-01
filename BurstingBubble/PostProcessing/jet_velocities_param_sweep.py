import numpy as np # type: ignore
from matplotlib import pyplot as plt # type: ignore
from matplotlib import collections  as mc # type: ignore
import matplotlib # type: ignore
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


# sets time inverval between files
dt = 0.1

scale=1.3
fig, ax = plt.subplots(figsize=(6.4*scale, 6.4*scale))

#variable for tracking number of files
colors=["r","g","b","purple"]
j=-1
for R in ["0.001","0.0008","0.0006","0.0004"]:
    j+=1
    total_points = 1
    previous=0
    velocities = []
    
    #R = "0.001"
    IC = "NN"
    Level = 10
    t = 10.0
    string = "%.1f" % float(total_points*0.1)
    file_name = "../DriverCode/Water-"+IC+"-R"+R+"-Level"+str(Level)+"/interfacestats.dat"
    try:
        data = np.loadtxt(file_name, delimiter=' ', usecols=(0,1,2,3,4,5,6,7), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
    IC = "IC"
    Level = 10
    t = 10.0
    string = "%.1f" % float(total_points*0.1)
    file_name = "../DriverCode/Water-"+IC+"-R"+R+"-Level"+str(Level)+"/interfacestats.dat"
    try:
        data2 = np.loadtxt(file_name, delimiter=' ', usecols=(0,1,2,3,4,5,6,7), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
    
      
    ### JET VEL ###



    """ current = 
    velocities.append((current-previous)/dt)
    previous = current """
    time = data[:,1]
    bottom_vals = data[:,3]
    top_vals = data[:,4]
    top_derv = np.zeros(np.size(time)-1)
    bottom_derv = np.zeros(np.size(time)-1)
    time_derv = np.zeros(np.size(time)-1)
    for i in range(np.size(time)-1):
        
      top_derv[i] = (top_vals[i+1]-top_vals[i])/(time[i+1]-time[i])
      bottom_derv[i] = (bottom_vals[i+1]-bottom_vals[i])/(time[i+1]-time[i])

      time_derv[i] = time[i]
    
    #ax.plot(data[:,1],data[:,3],label="min"+str(Level))
    ax.plot(data2[:,1],data2[:,4],label="Calculated, R="+R,color=colors[j],linewidth=3)
    ax.plot(data[:,1],data[:,4],label="Simplified, R="+R,linestyle=":",color=colors[j],linewidth=3)
    #ax.plot(data[:,1],data[:,7],label="interp"+str(Level))
    #ax.plot(time_derv,top_derv,label="min_derv"+str(Level))
    #ax.plot(time_derv,bottom_derv,label="max_derv"+str(Level))
    

plt.legend(fontsize=17)
ax.tick_params(axis='both', which='major', labelsize=17)
ax.tick_params(axis='both', which='minor', labelsize=8)
ax.set_xlabel('Time', fontsize=17)
ax.set_ylabel('Jet Height', fontsize=17)

plt.savefig('bubble_velocity_new.png')
plt.clf()

plt.close()
total_points+=1



