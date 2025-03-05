import numpy as np # type: ignore
import scipy.integrate as sci # type: ignore
from scipy.optimize import fsolve # type: ignore
from matplotlib import pyplot as plt # type: ignore


Theta = 4
delta = 0.001

def func_wrt_x(x,y):
  z=y[0]
  sin_theta = y[1]
  if x==0:
    return np.array([sin_theta/np.sqrt(1-sin_theta**2),Theta/2])
  return np.array([sin_theta/np.sqrt(1-sin_theta**2),2*z-sin_theta/x + Theta])

def func_wrt_z(z,y):
  x=y[0]
  cos_theta = y[1]
  if x==0:
    return np.array([cos_theta/np.sqrt(1-cos_theta**2),2*z + Theta])
  return np.array([cos_theta/np.sqrt(1-cos_theta**2),-2*z + np.sqrt(1-cos_theta**2)/x - Theta])

def func_C_eq(x,y):
  z=y[0]
  sin_theta = y[1]
  if x==0:
    return np.array([sin_theta/np.sqrt(1-sin_theta**2),0])
  return np.array([sin_theta/np.sqrt(1-sin_theta**2),2*z-sin_theta/x])

def func_B(x, R_0):
  return np.sqrt((R_0**2 - x**2))

def func_B_derv(x, R_0):
  return (-x)/(np.sqrt(R_0**2 - x**2))

def Full_system(x):
  z_bar = x[0]
  x_bar = x[1]
  h     = x[2]
  l     = x[3]
  R_0   = x[4]
  x_0   = x[5]
  #sol = sci.solve_ivp(func_wrt_x, [0,4], np.array([0,0]), method='RK45', rtol=1e-10, atol=1e-10)
  #next_x = sol.t[np.size(sol.t)-1]
  #next_z = sol.y[0,np.size(sol.y[0,:])-1]
  #sol2 = sci.solve_ivp(func_wrt_z, [next_z,next_z+4], np.array([next_x,0]), method='RK45', rtol=1e-10, atol=1e-10)
  f_1 = np.interp(x_bar,f_1_x,f_1_y)
  f_1_derv = (-np.interp(x_bar,f_1_x,f_1_y) + np.interp(x_bar+delta,f_1_x,f_1_y))/delta
  if x_0 == np.nan:
    x_0 = 13

  men_x_initial = x_0
  men_z_initial = 0.00001
  sintheta_initial = 0.00001
  men_sol = sci.solve_ivp(func_C_eq, [men_x_initial,0], np.array([men_z_initial,sintheta_initial]), method='RK45', rtol=1e-9, atol=1e-9)
  f_2_x = np.flip(men_sol.t)
  f_2_y = np.flip(men_sol.y[0,:])
  f_2 = np.interp(x_bar,f_2_x,f_2_y)
  f_2_derv = (-np.interp(x_bar,f_2_x,f_2_y)+np.interp(x_bar+delta,f_2_x,f_2_y))/delta
  f_3 = func_B(x_bar,R_0)
  f_3_derv = func_B_derv(x_bar,R_0)
  return np.array([f_1-h-z_bar,f_2-z_bar, f_3-l-z_bar,f_1_derv-f_3_derv, f_2_derv-f_3_derv, Theta-4/R_0 +2*h])

def Half_system(x,xmin, xmax):
  if (x[0]< xmin) or (x[0]> xmax):
    return np.array([1000,1000])
  x_bar = x[0]
  R_0   = x[1]
  h     = 2/R_0-Theta/2
  f_1 = np.interp(x_bar,f_1_x,f_1_y)
  f_1_derv = (-np.interp(x_bar,f_1_x,f_1_y) + np.interp(x_bar+delta,f_1_x,f_1_y))/delta
  f_3 = func_B(x_bar,R_0)
  f_3_derv = func_B_derv(x_bar,R_0)
  print((f_1-h),f_3,f_1_derv,f_3_derv)
  # need an l
  return np.array([(f_1-h)-(f_3),f_1_derv-f_3_derv]) 

scale=1
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))

#for Theta in [20]:#np.linspace(1,20,0.1):
sol = sci.solve_ivp(func_wrt_x, [0,4], np.array([0,0]), method='RK45', rtol=1e-5, atol=1e-10)
next_x = sol.t[np.size(sol.t)-1]
next_z = sol.y[0,np.size(sol.y[0,:])-1]

sol2 = sci.solve_ivp(func_wrt_z, [next_z,next_z+4], np.array([next_x,0]), method='RK45', rtol=1e-10, atol=1e-10)


    
    
    #ax.plot(sol.t,sol.y[0,:])
    #ax.plot(sol2.y[0,:], sol2.t)


""" SOLVING FOR A """

""" sol = sci.solve_ivp(func_wrt_x, [0,4], np.array([0,0]), method='RK45', rtol=1e-10, atol=1e-10)
next_x = sol.t[np.size(sol.t)-1]
next_z = sol.y[0,np.size(sol.y[0,:])-1]

sol2 = sci.solve_ivp(func_wrt_z, [next_z,next_z+4], np.array([next_x,0]), method='RK45', rtol=1e-10, atol=1e-10)
x_A = np.concatenate((sol.t,sol2.y[0,:]))
z_A = np.concatenate((sol.y[0,:],sol2.t)) """
""" SOLVING FOR C """
# Ideally we would have x=inf, z=0, phi=0
men_x_initial = 4
men_z_initial = 0.00001
sintheta_initial = 0.00001
#men_sol = sci.solve_ivp(func_C_eq, [men_x_initial,0], np.array([men_z_initial,sintheta_initial]), method='RK45', rtol=1e-9, atol=1e-9)
#ax.plot(men_sol.t, men_sol.y[0,:])

f_1_x = np.flip(sol2.y[0,:])
f_1_y = np.flip(sol2.t[:]) 


### z_bar, x_bar, h, l, R, x_0
initial_guess = np.array([
  0.13,
  0.284,
  0.592,
  0.592,
  0.772,
  0.596
])
root = fsolve(Full_system,initial_guess)
x_bar = root[1]
h = root[2]
l = root[3]
r_0 = root[4]
print(root)
""" 
 #[0.016,0.196]
root = root(Half_system,np.array([0.655,1.169]), args=(f_1_x[0],f_1_x[np.size(f_1_x)-1]),method='lm')
print(root.x)
print(Theta)
x_bar, r_0 = root.x
h     = 2/r_0-Theta/2"""
x = np.linspace(0,r_0,100)
test = np.linspace(0,r_0,100)
gap = 0#0.5035416392877905 - 1.174933217213838 
#ax.plot(x,func_B(x,r_0)+h-l)
#ax.plot(x_bar,func_B(x_bar,r_0)+h-l,'ro')

#############################################################################################
men_x_initial = root[5]
men_z_initial = 0.00001
sintheta_initial = 0.00001
men_sol = sci.solve_ivp(func_C_eq, [men_x_initial,0], np.array([men_z_initial,sintheta_initial]), method='RK45', rtol=1e-9, atol=1e-9)
f_2_x = np.flip(men_sol.t)
f_2_y = np.flip(men_sol.y[0,:])

##################################################################################################

ax.plot(np.concatenate((sol.t,sol2.y[0,:])),np.concatenate((sol.y[0,:],sol2.t))-h, label = "Rounded bottom")
ax.plot(x,func_B(x,r_0)-l, label="cap")
ax.plot(x_bar,func_B(x_bar,r_0)-l,'ro')
ax.plot(f_2_x,f_2_y)
#ax.set_xlim([0,1])
#np.savetxt("output.txt",[sol2.y[0,:50],sol2.t[:50]])
#ax.plot(men_sol.t, men_sol.y[0,:])
plt.legend()
plt.savefig('interface.png')
plt.clf() 
