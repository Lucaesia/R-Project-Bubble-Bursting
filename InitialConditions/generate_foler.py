import numpy as np # type: ignore
import scipy.integrate as sci # type: ignore
from scipy.optimize import fsolve # type: ignore
from matplotlib import pyplot as plt # type: ignore
from scipy.optimize import root_scalar # type: ignore
from scipy.optimize import minimize # type: ignore
from scipy.interpolate import PchipInterpolator # type: ignore
import matplotlib # type: ignore
matplotlib.rcParams['mathtext.fontset'] = 'cm'

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

def Full_system_2_var(x):
  x_bar = x[0]
  x_0 = x[1]
  
  men_x_initial = x_0
  men_z_initial = 0.00001
  sintheta_initial = 0.00001
  men_sol = sci.solve_ivp(func_C_eq, [men_x_initial,0], np.array([men_z_initial,sintheta_initial]), method='RK45', rtol=1e-9, atol=1e-9)
  f_2_x = np.flip(men_sol.t)
  f_2_y = np.flip(men_sol.y[0,:])
  f_2 = np.interp(x_bar,f_2_x,f_2_y)
  f_2_derv = (-np.interp(x_bar,f_2_x,f_2_y)+np.interp(x_bar+delta,f_2_x,f_2_y))/delta

  f_1 = np.interp(x_bar,f_1_x,f_1_y)
  f_1_derv = (-np.interp(x_bar,f_1_x,f_1_y) + np.interp(x_bar+delta,f_1_x,f_1_y))/delta

  h = f_1 - f_2
  R_0 = 2/(h + Theta/2)
  #l = func_B(x_bar,R_0)-f_1

  
  #f_3 = func_B(x_bar,R_0)
  f_3_derv = func_B_derv(x_bar,R_0)
  return abs(f_1_derv-f_3_derv)**2 + abs(f_2_derv-f_3_derv)**2

def variable_from_two(x):
  x_bar = x[0]
  x_0 = x[1]
  
  
  men_x_initial = x_0
  men_z_initial = 0.00001
  sintheta_initial = 0.00001
  men_sol = sci.solve_ivp(func_C_eq, [men_x_initial,0], np.array([men_z_initial,sintheta_initial]), method='RK45', rtol=1e-9, atol=1e-9)
  f_2_x = np.flip(men_sol.t)
  f_2_y = np.flip(men_sol.y[0,:])
  f_2 = np.interp(x_bar,f_2_x,f_2_y)
  f_2_derv = (-np.interp(x_bar,f_2_x,f_2_y)+np.interp(x_bar+delta,f_2_x,f_2_y))/delta

  f_1 = np.interp(x_bar,f_1_x,f_1_y)
  f_1_derv = (-np.interp(x_bar,f_1_x,f_1_y) + np.interp(x_bar+delta,f_1_x,f_1_y))/delta

  h = f_1 - f_2
  R_0 = 2/(h + Theta/2)
  l = func_B(x_bar,R_0)-f_2
  #print("z=",f_1,f_2,func_B(x_bar,R_0))

  return np.array([h,l,R_0])

def Full_system_1_var(x_0):
  men_x_initial = x_0
  men_z_initial = 0.00001
  sintheta_initial = 0.00001
  men_sol = sci.solve_ivp(func_C_eq, [men_x_initial,0], np.array([men_z_initial,sintheta_initial]), method='RK45', rtol=1e-9, atol=1e-9)
  f_2_x = np.flip(men_sol.t)
  f_2_y = np.flip(men_sol.y[0,:])
  f_2 = np.interp(x_bar,f_2_x,f_2_y)

  f_1_derv_func = lambda t : (-np.interp(t,f_1_x,f_1_y) + np.interp(t+delta,f_1_x,f_1_y))/delta
  f_2_derv_func = lambda t : (-np.interp(t,f_2_x,f_2_y)+np.interp(t+delta,f_2_x,f_2_y))/delta
  root_func = lambda t : f_1_derv_func(t)-f_2_derv_func(t)

  x_bar = root_scalar(root_func,bracket=[f_2_x[0]+0.0001, f_2_x[np.size(f_2_x)-1]-0.0001])
  f_2_derv = (-np.interp(x_bar,f_2_x,f_2_y)+np.interp(x_bar+delta,f_2_x,f_2_y))/delta

  f_1 = np.interp(x_bar,f_1_x,f_1_y)
  f_1_derv = (-np.interp(x_bar,f_1_x,f_1_y) + np.interp(x_bar+delta,f_1_x,f_1_y))/delta

  h = f_1 - f_2
  R_0 = 2/(h + Theta/2)
  #l = func_B(x_bar,R_0)-f_1

  
  #f_3 = func_B(x_bar,R_0)
  f_3_derv = func_B_derv(x_bar,R_0)
  return np.array([f_1_derv-f_3_derv])

#=================================#
# Printing out graphs of Under    #
#=================================#

""" scale=1.2
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))
for Theta in [4]:
    sol = sci.solve_ivp(func_wrt_x, [0,4], np.array([0,0]), method='RK45', rtol=1e-10, atol=1e-10)
    next_x = sol.t[np.size(sol.t)-1]
    next_z = sol.y[0,np.size(sol.y[0,:])-1]

    sol2 = sci.solve_ivp(func_wrt_z, [next_z,next_z+4], np.array([next_x,0]), method='RK45', rtol=1e-10, atol=1e-10)
    f_1_x = np.flip(sol2.y[0,:])
    f_1_y = np.flip(sol2.t[:]) 

    
    
    #ax.plot(np.concatenate ((sol.t,sol2.y[0,:])),np.concatenate ((sol.y[0,:],sol2.t)),label="Theta = "+str(Theta))
    #ax.plot(sol2.y[0,:], sol2.t)
 """
#=================================#
# Printing out graphs of meniscus #
#=================================#
""" for x_bar_loop in [4,4.2,4.4,4.7,5,5.5,6]:
  men_x_initial = x_bar_loop
  men_z_initial = 0.00001
  sintheta_initial = 0.00001
  men_sol = sci.solve_ivp(func_C_eq, [men_x_initial,0], np.array([men_z_initial,sintheta_initial]), method='RK45', rtol=1e-9, atol=1e-9)
  f_2_x = np.flip(men_sol.t)
  f_2_y = np.flip(men_sol.y[0,:])
  ax.plot(f_2_x, f_2_y, label="$\overline{x}_0 = " + str(x_bar_loop)+"$")

ax.set_aspect('equal')
ax.set_xlabel('$\overline{x}$', fontsize=12)
ax.set_ylabel('$\overline{z}$', rotation=0, fontsize=12)
ax.set_xlim([0,0.1])
plt.legend()
plt.subplots_adjust(bottom=0.15)
plt.savefig('menisdus_interface.png')
plt.clf()
plt.close() 
exit()  """

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


initial_guess = np.array([##for 1.6
  0.8034748,  8.79743954
])

initial_guesses = np.array([
  [0.8034748,  8.79743954],
  [0.58194697, 8.34403281],
  [0.41497931, 7.89626882],
  [0.3088965,  7.51300205],
  [0.1215658,  6.55189027],
  [0.03871182, 5.78062276]
])
Theta = 2.3
scale=2.1
fig, ax = plt.subplots(figsize=(4*scale, 4*scale))

initial_guess = np.array([##for 2.3
  0.58194697, 8.34403281
])
""" initial_guess = np.array([##for 2.63
  0.51888806, 8.16269956
])  """
""" initial_guess = np.array([##for 3
  0.34733607, 7.76091227
])  """
""" initial_guess = np.array([##for 3.55
  0.3389654,  7.63950984
])  """
""" initial_guess = np.array([##for 4.25
  0.3389654,  7.43950984
])  """ 
""" initial_guess = np.array([##for 5.1
  0.23764154, 7.37542487
])  """
""" initial_guess = np.array([##for 6.5
  0.1772592,  7.13220278
])   """
""" initial_guess = np.array([##for 8.9
  0.11988799, 6.73349939
])   """ 
""" initial_guess = np.array([##for 13.5
  0.07808964, 6.50660065
])   """
""" initial_guess = np.array([##for 27
  0.01110725, 4.50652028
])  """

sol = sci.solve_ivp(func_wrt_x, [0,4], np.array([0,0]), method='RK45', rtol=1e-10, atol=1e-10)
next_x = sol.t[np.size(sol.t)-1]
#print(next_x*2.71)
#exit()
next_z = sol.y[0,np.size(sol.y[0,:])-1]
sol2 = sci.solve_ivp(func_wrt_z, [next_z,next_z+4], np.array([next_x,0]), method='RK45', rtol=1e-10, atol=1e-10)
f_1_x = np.flip(sol2.y[0,:])
f_1_y = np.flip(sol2.t[:]) 
string = "%.2f" % float(Theta) ## method="trust-constr"
sol_root = minimize(Full_system_2_var,initial_guess, bounds = ((f_1_x[0]+0.001, f_1_x[np.size(f_1_x)-1]-0.001), (1, None)))
#############################################################################################
men_x_initial = sol_root.x[1]
men_z_initial = 0.00001
sintheta_initial = 0.00001
men_sol = sci.solve_ivp(func_C_eq, [men_x_initial,0], np.array([men_z_initial,sintheta_initial]), method='RK45', rtol=1e-9, atol=1e-9)
f_2_x = np.flip(men_sol.t)
f_2_y = np.flip(men_sol.y[0,:])


x_bar = sol_root.x[0]
h = variable_from_two(sol_root.x)[0]
l = variable_from_two(sol_root.x)[1]
r_0 = variable_from_two(sol_root.x)[2]
x = np.linspace(0,r_0,100)
##################################################################################################

ax.plot(np.concatenate((sol.t,sol2.y[0,:])),np.concatenate((sol.y[0,:],sol2.t))-h, label = "Rounded bottom")
ax.plot(x,func_B(x,r_0)-l, label="cap")
ax.plot(x_bar,func_B(x_bar,r_0)-l,'ro')
ax.plot(f_2_x,f_2_y)

""" for i in range(10):
  men_x_initial += 0.5
  men_sol = sci.solve_ivp(func_C_eq, [men_x_initial,0], np.array([men_z_initial,sintheta_initial]), method='RK45', rtol=1e-9, atol=1e-9)
  f_2_x = np.flip(men_sol.t)
  f_2_y = np.flip(men_sol.y[0,:])
  ax.plot(f_2_x,f_2_y,label=str(men_x_initial)) """
ax.set_xlim([0,1.51])
ax.set_aspect('equal')
#np.savetxt("output.txt",[sol2.y[0,:50],sol2.t[:50]])
#ax.plot(men_sol.t, men_sol.y[0,:])
plt.legend()
plt.savefig('interface.png')
plt.clf()
plt.close() 
print(sol_root.x)  
print(Full_system_2_var(sol_root.x))

####################################### POLYNOMIAL CREATION ######################################
print("polystart")
f_1_x = np.flip(sol2.y[0,:])
f_1_y = np.flip(sol2.t[:]) 

under_coof = PchipInterpolator(sol.t,sol.y[0,:]-h)
over_coof = PchipInterpolator(f_1_x,f_1_y-h)
miniscus_coof = PchipInterpolator(f_2_x,f_2_y)
string = "%.1f" % float(Theta)

scale = 1/next_x
np.savetxt("theta"+string+"/under.dat",scale*np.transpose(np.stack((sol.t,sol.y[0,:]-h))))

np.savetxt("theta"+string+"/over.dat",scale*np.transpose(np.stack((f_1_x,f_1_y-h))))

np.savetxt("theta"+string+"/menisc.dat",scale*np.transpose(np.stack((f_2_x,f_2_y))))
#                                                                x_R          menisc_x    x_bar  under_N, Over_N, menisc_N
np.savetxt("theta"+string+"/RL.dat",scale*np.transpose(np.array([sol2.y[0,0],men_sol.t[0],x_bar])))

with open("theta"+string+"/RL.dat","a") as my_file:
    my_file.write(str(np.size(sol.t))+"\n")
    my_file.write(str(np.size(f_1_x))+"\n")
    my_file.write(str(np.size(f_2_x))+"\n")
    my_file.write(str(scale*(sol.y[0,0]-h))+"\n")
    my_file.write(str(scale*(f_2_y[0]))+"\n")
    my_file.write(str(scale)+"\n")
    my_file.write(str(scale*r_0)+"\n")
    my_file.write(str(scale*l)+"\n")
    my_file.write(str(scale*h)+"\n")

right_most_x = sol.t[np.size(sol.t)-1]
under_coof_x_list = np.linspace(0,right_most_x,100)
over_coof_x_list = np.linspace(x_bar,right_most_x,100)
miniscus_coof_x_list = np.linspace(x_bar,4,100)


scale=1
fig, ax = plt.subplots(figsize=(6.4*scale, 3.2*scale))

ax.plot(under_coof_x_list,under_coof(under_coof_x_list), label = "Rounded bottom")
ax.plot(over_coof_x_list,over_coof(over_coof_x_list), label = "Rounded side")
ax.plot(miniscus_coof_x_list,miniscus_coof(miniscus_coof_x_list), label = "minscus")
ax.plot(x,func_B(x,r_0)-l, label="cap")
ax.plot(x_bar,func_B(x_bar,r_0)-l,'ro')

ax.set_xlim([0,1.11])
ax.set_aspect('equal')
#np.savetxt("output.txt",[sol2.y[0,:50],sol2.t[:50]])
#ax.plot(men_sol.t, men_sol.y[0,:])
plt.legend()
#plt.savefig('interface_poly.png')
plt.clf() 

#######################################---------------------######################################
