import numpy as np # type: ignore
import scipy.integrate as sci # type: ignore
from matplotlib import pyplot as plt # type: ignore
import matplotlib # type: ignore

matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def func_B(x, R_0):
  return np.sqrt((R_0**2 - (x)**2))

def func_B_inv(y,R_0):
  return np.sqrt((R_0**2 - (y+l)**2))

def cube(x):
  return x*x*x

def top_line(x):
  if x<x_bar:
    return func_B(x, R_0) -l
  return np.interp(x,over[:,0],over[:,1])

def bottom_line(x):
 
  return np.interp(x,under[:,0],under[:,1])

def func_for_y(y):
  if y > R_0 - l:
    return 0
  if y>f_x_bar:
    return func_B_inv(y,R_0)
  if y > over[np.size(over)//2 -1,1]:
    return np.interp(y,np.flip(over[:,1]),np.flip(over[:,0]))
  if y > under[0,1]:
    return np.interp(y,under[:,1],under[:,0])
  return 0
def f(y):
  return func_for_y(y)*func_for_y(y)
def in_or_out (X):
  x = X[0]
  y = X[1]

  if x> x_R:
    return 0
  
  if x>x_bar and x<x_R:
    val = np.interp(x,over[:,0],over[:,1])
    val2 = np.interp(x,under[:,0],under[:,1])
    
    if y<val and y>val2:
      return 1
  if x< x_bar:
    val_higher = func_B(x, R_0) -l
    val_lower = np.interp(x,under[:,0],under[:,1])
    if y<val_higher and y>val_lower:
      return 1
  return 0

def volume_calc():
  count = 0

  for x in np.linspace(0,1,1000):
    for y in np.linspace(-3,2,1000):
      count+= in_or_out([x,y])
  volume = 5*count/(1000*1000)
  print(volume)
  d2_area = 2.71*2.71*volume*2/(scale_x*scale_x)
  return 2.71*2.71*volume*2/(scale_x*scale_x)

Thetas = ["2.3","2.6","3.0","3.5","4.2","5.1","6.5","8.9", "13.5","27.0"]
volumes = []
for i in Thetas:
  string = "theta"+i
  menisc = np.loadtxt(string+"/menisc.dat")
  over  = np.loadtxt(string+"/over.dat")
  RL =  np.loadtxt(string+"/RL.dat")
  under =  np.loadtxt(string+"/under.dat")


  x_R = RL[0]
  x_far = RL[1]
  x_bar = RL[2]
  scale_x = RL[8]
  R_0 = RL[9]
  l = RL[10]
  h = RL[11]
  f_x_bar = np.interp(x_bar,over[:,0],over[:,1])
  volumes.append(np.pi*sci.quad(f, -3, 2)[0]*2.71*2.71*2.71/(scale_x*scale_x*scale_x))
  
count = 0
### FROM PAST RUN
print(volumes)


scale=1
fig, ax = plt.subplots(figsize=(6.4*scale, 5*scale))

#ax.plot(np.arange(1.6,31.6,0.2),array*0.3848)


A = volumes[-1]/0.04
x = np.linspace(0.2,2,100)

ax.plot(x,(4/3)*np.pi*cube(x),label="Spherical Volume")
ax.plot(np.arange(2,0,-0.2),volumes, label="Calculated Volume",linestyle=":")
""" y_list=[]

for i in np.linspace(-2,2,100):
  y_list.append(func_for_y(i))
ax.plot(np.linspace(-2,2,100),y_list, label="Calculated Volume",linestyle=":")
ax.plot(over[:,0],over[:,1])
ax.plot(np.linspace(-2,2,100),func_B(np.linspace(-2,2,100),R_0)-l) """
#ax.set_xlim([-3,3])
#ax.set_ylim([-0.8,0.7])
#ax.set_aspect('equal')
#np.savetxt("output.txt",[sol2.y[0,:50],sol2.t[:50]])
#ax.plot(men_sol.t, men_sol.y[0,:])
ax.set_ylabel('Volume', fontsize=15)
ax.set_xlabel('Major Radius of bubble', fontsize=15)
plt.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=8)

ax.set_xlabel('$R$', fontsize=17)
ax.set_ylabel('$V$', rotation=0, fontsize=17)

#plt.legend()
plt.savefig('volume.png')
plt.clf()
plt.close()  