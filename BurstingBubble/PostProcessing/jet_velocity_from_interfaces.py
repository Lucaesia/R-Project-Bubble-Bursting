import numpy as np # type: ignore
from matplotlib import pyplot as plt # type: ignore
from matplotlib import collections  as mc # type: ignore
from scipy import interpolate
from scipy import signal

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


# sets time inverval between files
dt = 0.1

scale=1
fig, ax = plt.subplots(figsize=(6.4*scale, 6.4*scale))

#variable for tracking number of files
for i in [6,7,8,9,10,11]:
    total_points = 1
    previous=0
    velocities = []
    
    R = "0.001"
    IC = "NN"
    Level = i
    t = 10.0
    string = "%.2f" % float(total_points*0.01)
    file_name = "../DriverCode/Water-"+IC+"-R"+R+"-Level"+str(Level)+"/interfacestats.dat"
    try:
        data = np.loadtxt(file_name, delimiter=' ', usecols=(0,1,2,3,4,5,6,7), unpack=False)
    except FileNotFoundError:
        print("file not found"+file_name)
        break
    
      
    ### JET VEL ###



    """ current = 
    velocities.append((current-previous)/dt)
    previous = current """
    time = data[:,1]
    top_vals = data[:,3]
    bottom_vals = data[:,4]
    top_derv = np.zeros(np.size(time)-1)
    bottom_derv = np.zeros(np.size(time)-1)
    time_derv = np.zeros(np.size(time)-1)


    top_vals = signal.savgol_filter(top_vals,30,3)
    interp = interpolate.CubicSpline(time,top_vals)
    pp=interp.derivative()
    

    #top_derv=np.convolve(top_derv,np.random.normal(size=np.size(top_derv)),mode='same')
    
    #ax.plot(data[:,1],data[:,3],label="min"+str(Level))
    #ax.plot(data[:,1],data[:,3],label="max"+str(Level))
    #ax.plot(data[:,1],top_vals,label="max smooth"+str(Level))
    #ax.plot(data[:,1],pp(data[:,1]),label="max_smooth"+str(Level))
    ax.plot(data[:,1],data[:,7],label="interp"+str(Level))
    #ax.plot(time_derv,top_derv,label="min_derv"+str(Level))
    #ax.plot(time_derv,bottom_derv,label="max_derv"+str(Level))

ax.set_xlim([0.14,0.25])
ax.set_ylim([-1,50])
plt.legend()
plt.savefig('bubble_velocity_new.png')
plt.clf()

plt.close()
total_points+=1



