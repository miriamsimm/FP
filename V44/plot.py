import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp

from scipy.optimize import curve_fit

def gauss(x,mu,A,sigma_sqr):
    return A*np.exp(- (   (x -   mu)**2 )   /     (    2 *sigma_sqr  )      )

data = np.genfromtxt("data/detectorscan.UXD")
x = data[:,0]
y = data[:,1]
xplot = np.linspace(-0.6,0.6,10000)
param,cova = curve_fit(gauss, x, y)
error = np.sqrt(np.diag(cova))
plt.plot(xplot,gauss(xplot,param[0],param[1],param[2]),label = "Fitkurve",color = "orange")
plt.plot(x,y,"x",label = "Datenpunkte",color = "skyblue")
plt.xlim(-0.52,0.52)
plt.ylim(-2e4,1e6)
plt.xlabel(r"$\theta$[°]")
plt.ylabel(r"Intensität  [rel. Ein.]")
sigma = unp.uarray(param[2],error[2])
sigma = unp.sqrt(sigma)
print("mu = ",param[0],r"$\pm$",error[0],error[0]/param[1])
print("A = ",param[1],r"$\pm$",error[1],error[1]/param[1])
print("sigma = ",sigma,unp.std_devs(sigma)/unp.nominal_values(sigma))
print("FWHM= ",2*np.sqrt(2*np.log(2))*sigma  )
plt.legend()
plt.savefig("plots/detektor.pdf")
plt.close()

data = np.genfromtxt("data/z2.UXD")
x = data[:,0]
y = data[:,1]
plt.plot(x,y,ms = 1.5,label = "z-scan",color = "blue")
plt.vlines(-0.08, -20000, 1e7,color = "red",linestyles="dotted")
plt.vlines(0.16, -20000, 1e7,color = "red",linestyles="dotted")
plt.ylim(-20000,1e6)
plt.xlim(-0.5,0.5)
plt.xlabel(r"$\theta$[°]")
plt.ylabel(r"Intensität  [rel. Ein.]")
plt.legend()
plt.grid()
plt.savefig("plots/z_scan.pdf")
plt.close()

data = np.genfromtxt("data/rocking1.UXD")
x = data[:,0]
y = data[:,1]
plt.plot(x,y,"x",ms = 1.5,label = "rocking",color = "blue")
plt.vlines(-0.7, -20000, 1e7,color = "red",linestyles="dotted")
plt.vlines(0.7, -20000, 1e7,color = "red",linestyles="dotted")
plt.ylim(-10000,0.55e6)
plt.xlim(-1,1)
plt.grid()
plt.xlabel(r"$\theta$[°]")
plt.ylabel(r"Intensität  [rel. Ein.]")
plt.legend()
plt.savefig("plots/rocking.pdf")
plt.close()


theta_g = 0.7
d_0 = 0.24
D = d_0/np.sin(np.deg2rad(theta_g)) 

data = np.genfromtxt("data/scan1.UXD")
x = data[:,0]
y = data[:,1]
data2 = np.genfromtxt("data/scan2-diffus.UXD")
x2 = data2[:,0]
y2 = data2[:,1]

plt.plot(x,y,lw = 1,label = "Reflektivitäts Scan",color = "blue")
plt.plot(x2,y2,lw = 1,label = "Diffuser Scan",color = "orange")
y = y -y2
plt.plot(x,y,lw = 1,label = "Differenz",color = "red")
plt.yscale("log")
plt.ylim(0.7,1e7)
plt.xlim(0,2.5)
plt.xlabel(r"$\theta$[°]")
plt.ylabel(r"Intensität  [rel. Ein.]")
plt.legend()
plt.savefig("plots/refklekt.pdf")
plt.close()

lam = 1.54*1e-10
xplot = np.linspace(0,2.5,10000)
def silizium(alpha):
    alpha = np.deg2rad(alpha)
    n0 = 1
    n1 = 1- 7.6e-6 + 1.73j*1e-7
    k0 = 2*np.pi/lam *np.sqrt(n0**2 - np.cos(alpha)**2)
    k1 = 2*np.pi/lam *np.sqrt(n1**2 -np.cos(alpha)**2 )
    return np.abs((k0 - k1 )/( k0 + k1 ))**2
plt.plot(xplot,param[1]*silizium(xplot),label = "Ideale Siliziumschicht",color = "orange")
plt.plot(x,y,label = "Differenz",color = "blue")
plt.xlim(0,2.5)
plt.ylim(1e-1,1e7)
plt.xlabel(r"$\theta$[°]")
plt.ylabel(r"Intensität  [rel. Ein.]")
plt.yscale("log")
plt.legend(loc = "best")
plt.savefig("plots/ideal.pdf")
plt.close()


for i in range(np.size(x)):
    if(x[i] < theta_g):
        if(np.sin(np.deg2rad(x[i]))  != 0):
            y[i] = np.sin(np.deg2rad(theta_g))*y[i]/np.sin(np.deg2rad(x[i]))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Major ticks every 20, minor ticks every 5
major_ticks = np.arange(0, 1, 0.1)
minor_ticks = np.arange(0, 1, 0.01)

ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(major_ticks)
ax.set_yticks(minor_ticks, minor=True)

# And a corresponding grid

ax.plot(x,y,lw = 1,label = "Differenz",color = "red")
plt.yscale("log")
plt.ylim(1e2,2e4)
plt.xlim(0.4,1.0)
plt.xlabel(r"$\theta$[°]")
plt.ylabel(r"Intensität  [rel. Ein.]")
ax.grid(which='both')
plt.savefig("plots/refklekt_norm_zoom.pdf")
plt.close()



theta_min = np.array([0.44,0.49,0.54,0.59,0.64,0.7,0.75,0.8,0.85,0.91,0.96])
theta_min = np.deg2rad(theta_min)
print("Number of min: ",np.size(theta_min))
print("minima: ",theta_min)
print("Delta_Theta","Delta_q","d")
delta_theta = np.zeros(np.size(theta_min)-1)
delta_q = np.zeros(np.size(theta_min)-1)
d = np.zeros(np.size(theta_min)-1)
for i in range(np.size(theta_min)-1):
    delta_theta[i] = theta_min[i+1]-theta_min[i]
    delta_q[i] = 4*np.pi*( delta_theta[i]  ) / lam
    d[i] = lam/(2 *delta_theta[i] )
    print(np.rad2deg(delta_theta[i])," " ,  delta_q[i], " ",d[i])
    
z0 = np.sum(d)/np.size(d)

print(np.rad2deg(np.mean(delta_theta)),np.rad2deg(np.std(delta_theta)))
print(np.mean(delta_q),np.std(delta_q))
print(np.mean(d),np.std(d),np.std(d)/np.mean(d))
#,np.sum(delta_q)/np.size(delta_q),np.sum(d)/np.size(d))




#factor = y[x == 0.1950]   # Faktor um die Messwerte zu normalisieren
factor = param[1]
y = y/factor
def parrat(alpha,z,delta1,delta2,beta1,beta2,sigma0,sigma1):
    alpha = np.deg2rad(alpha)
    n0 = 1
    n1 = 1- delta1 + 1j*beta1
    n2 = 1- delta2 + 1j*beta2
    lam = 1.54*1e-10
    k0 = 2*np.pi/lam *np.sqrt(n0**2 - np.cos(alpha)**2)
    k1 = 2*np.pi/lam *np.sqrt(n1**2 -np.cos(alpha)**2 )
    k2 = 2*np.pi/lam *np.sqrt(n2**2 -np.cos(alpha)**2 )

    r0 = (k0 - k1 )/( k0 + k1 )*np.exp(    -2*k0*k1*sigma0**2     )
    r1 = (k1 - k2 )/( k1 + k2 )*np.exp(    -2*k1*k2*sigma1**2     )
    X1 = np.exp(-2j*k1*z )*r1
    X0 = (  r0  + X1   )/  ( 1+ r0*X1   )
    return np.abs(X0)**2
literatur = np.array([3.5e-6,7.6e-6,4/(1e-2),141/(1e-2),0.153,0.223])
literatur[2] = lam*literatur[2]/(4*np.pi)
literatur[3] = lam*literatur[3]/(4*np.pi)
print("Literatur:",literatur)
print(" FIT 1:")

bounds = ( [1e-10,1e-7,1e-8,1e-8,1e-8,1e-10,1e-10], [1e-6,1e-5,1e-5,1e-5,1e-5,1e-9,1e-9])
param_0 = [9e-8,1e-6,1e-7,2e-6,1e-7,3e-10,8e-10]
xplot = np.linspace(0,3,10000)
par_param,par_cova = curve_fit(parrat, x, y ,p0 = param_0,bounds= bounds)
error_par = np.sqrt(np.diag(par_cova))
for i in range (0,np.size(par_param) ):
    print(par_param[i],"+-",error_par[i], "rel. err:",error_par[i]/par_param[i])
delta1 = unp.uarray(par_param[1],error_par[1])
delta2 = unp.uarray(par_param[2],error_par[2])
theta_c_1 = 180/np.pi*(unp.sqrt(2*delta1))
theta_c_2 = 180/np.pi*(unp.sqrt(2*delta2))
print("theta_c_1:",theta_c_1,"LIT_Abw:",(theta_c_1 - literatur[4] )/literatur[4] )
print("theta_c_2:",theta_c_2,"LIT_Abw:",(theta_c_2 - literatur[5] )/literatur[5] )
print("Lit_Abw:","delta1",(par_param[0] -literatur[0])/literatur[0] )
print("Lit_Abw:","delta2",(par_param[1] -literatur[1])/literatur[1] )
print("Lit_Abw:","beta1",(par_param[2] -literatur[2])/literatur[2] )
print("Lit_Abw:","beta1",(par_param[3] -literatur[3])/literatur[3] )

plt.plot(xplot,parrat(xplot,*par_param),label = "Fitkurve",color = "orange")
plt.plot(x,y,label = "normierte Messkurve",color = "blue")
plt.xlim(0,2.5)
plt.ylim(1e-7,100)
plt.yscale("log")
plt.xlabel(r"$\theta$[°]")
plt.ylabel(r"Intensität  [rel. Ein.]")
plt.legend(loc = "best")
plt.savefig("plots/parrat.pdf")
plt.close()

bounds2 = ( [1e-10,1e-7,1e-8,1e-8,1e-8,1e-10], [1e-6,1e-5,1e-5,1e-5,1e-5,1e-9])
param_0 = [9e-8,1e-6,1e-7,2e-6,1e-7,3e-10]
def parrat_2(alpha,z,delta1,delta2,beta1,beta2,sigma1):
    return parrat(alpha,z,delta1,delta2,beta1,beta2,sigma1,3*sigma1)

print(" FIT 2:")
par_param,par_cova = curve_fit(parrat_2, x, y ,p0 = param_0,bounds= bounds2)
error_par = np.sqrt(np.diag(par_cova))
for i in range (0,np.size(par_param) ):
    print(par_param[i],"+-",error_par[i], "rel. err:",error_par[i]/par_param[i])
delta1 = unp.uarray(par_param[1],error_par[1])
delta2 = unp.uarray(par_param[2],error_par[2])
theta_c_1 = 180/np.pi*(unp.sqrt(2*delta1))
theta_c_2 = 180/np.pi*(unp.sqrt(2*delta2))
print("theta_c_1:",theta_c_1,"LIT_Abw:",(theta_c_1 - literatur[4] )/literatur[4] )
print("theta_c_2:",theta_c_2,"LIT_Abw:",(theta_c_2 - literatur[5] )/literatur[5] )
print("Lit_Abw:","delta1",(par_param[0] -literatur[0])/literatur[0] )
print("Lit_Abw:","delta2",(par_param[1] -literatur[1])/literatur[1] )
print("Lit_Abw:","beta1",(par_param[2] -literatur[2])/literatur[2] )
print("Lit_Abw:","beta1",(par_param[3] -literatur[3])/literatur[3] )
plt.plot(xplot,parrat_2(xplot,*par_param),label = "Fitkurve",color = "orange")
plt.plot(x,y,label = "normierte Messkurve",color = "blue")
plt.xlim(0,2.5)
plt.ylim(1e-7,100)
plt.yscale("log")
plt.xlabel(r"$\theta$[°]")
plt.ylabel(r"Intensität  [rel. Ein.]")
plt.legend(loc = "best")
plt.savefig("plots/parrat3.pdf")
plt.close()



bounds2 = ( [1e-10,1e-7,1e-8,1e-8,1e-8], [1e-6,1e-5,1e-5,1e-5,1e-5])
param_0 = [9e-8,1e-6,1e-7,2e-6,1e-7]
def parrat_3(alpha,z,delta1,delta2,beta1,beta2):
    return parrat(alpha,z,delta1,delta2,beta1,beta2,0,0)

print(" FIT 3:")
par_param,par_cova = curve_fit(parrat_3, x, y ,p0 = param_0,bounds= bounds2)
error_par = np.sqrt(np.diag(par_cova))
for i in range (0,np.size(par_param) ):
    print(par_param[i],"+-",error_par[i], "rel. err:",error_par[i]/par_param[i])
delta1 = unp.uarray(par_param[1],error_par[1])
delta2 = unp.uarray(par_param[2],error_par[2])
theta_c_1 = 180/np.pi*(unp.sqrt(2*delta1))
theta_c_2 = 180/np.pi*(unp.sqrt(2*delta2))
print("theta_c_1:",theta_c_1,"LIT_Abw:",(theta_c_1 - literatur[4] )/literatur[4] )
print("theta_c_2:",theta_c_2,"LIT_Abw:",(theta_c_2 - literatur[5] )/literatur[5] )
print("Lit_Abw:","delta1",(par_param[0] -literatur[0])/literatur[0] )
print("Lit_Abw:","delta2",(par_param[1] -literatur[1])/literatur[1] )
print("Lit_Abw:","beta1",(par_param[2] -literatur[2])/literatur[2] )
print("Lit_Abw:","beta1",(par_param[3] -literatur[3])/literatur[3] )
plt.plot(xplot,parrat_3(xplot,*par_param),label = "Fitkurve",color = "orange")
plt.plot(x,y,label = "normierte Messkurve",color = "blue")
plt.xlim(0,2.5)
plt.ylim(1e-7,100)
plt.yscale("log")
plt.xlabel(r"$\theta$[°]")
plt.ylabel(r"Intensität  [rel. Ein.]")
plt.legend(loc = "best")
plt.savefig("plots/parrat3.pdf")
plt.close()



