import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.umath as unp
N = np.genfromtxt("data/muon_data.txt")
print(np.size(N) )

kali =np.genfromtxt("data/kalibration.txt") 

c = kali[:,1]
t = kali[:,0]
print(c)
print(t)
k_par,k_cov = np.polyfit(c,t,1,cov = True)
k_error = np.sqrt(np.diag(k_cov))
kali_A = ufloat(k_par[0],k_error[0])
kali_B = ufloat(k_par[1],k_error[1])
p = np.poly1d(k_par)
print("A",kali_A)
print("B",kali_B)

xplot = np.linspace(0,250,10000)
plt.plot(c,t,"x",label = "Messwerte")
plt.plot(xplot,p(xplot),label = "Fitgerade")
plt.xlim(0,250)
plt.ylim(0,12)
plt.xlabel("Channel")
plt.ylabel(r"$t / \mu s$")
plt.legend(loc = "best")
plt.savefig("kalibration.pdf")
plt.close()

Events = np.zeros(224)
t = np.zeros(224)
for i in range(224):
    t[i] = p(i+5)
    if(N[i+4] == 1):
        Events[i] = 1
    else:
        Events[i] =N[i+4] -1
        print(t[i],Events[i])


par,cov = np.polyfit(t,np.log(Events),1,cov = True)
error = np.sqrt(np.diag(cov))
p2 = np.poly1d(par)
print(p2)
tau = ufloat(par[0],error[0])
tau = 1/tau
N_0 = ufloat(par[1],error[1])
N_0 = unp.exp(N_0)
plt.plot(t,Events,"x",label = "Messdaten")
plt.plot(xplot,np.exp(p2(xplot)),label = "Fitgerade")
plt.yscale("log")
plt.xlim(0.6,10.5)
plt.ylim(1,3*100)
plt.xlabel(r"$t / \mu s $")
plt.ylabel(r"Counts")
plt.legend(loc = "best")
plt.savefig("lebensdauer.pdf")
plt.close()
print("Lebensdauer",tau)
print("N_0=",N_0)
def Gauss(x,A,mu,sigma):
    return A*np.exp(-   (x-mu)**2 / (2*sigma**2)       )
koin =np.genfromtxt("data/koinzidenz.txt")
ydata = koin[:,1]
xdata = koin[:,0]
print(xdata)

par_koin,cov_koin  = curve_fit(Gauss, xdata, ydata) 
err_koin = np.sqrt(np.diag(cov_koin))
sigma = ufloat(par_koin[2],err_koin[2])

print(par_koin[0],"+-",err_koin[0])
print(par_koin[1],"+-",err_koin[1])
print(par_koin[2],"+-",err_koin[2])
print(2*np.sqrt(2*np.log(2))*sigma)
xplot = np.linspace(-20,20,10000)

plt.plot(xdata,ydata,"x",label = "Messwerte")
#plt.plot(xplot,Gauss(xplot, *par_koin),color = "blue",label = "Fitkurve ")
plt.xlabel(r"Verzögerungszeit $T_V$ / ns")
plt.ylabel("Counts")
plt.xlim(-13,13)
plt.ylim(0,300)
xplot1 = np.linspace(-13,-2,1000)
xplot2 = np.linspace(-2,4,1000)
xplot3 = np.linspace(4,13,1000)

par_min,cov_par_min = np.polyfit(xdata[0:6],ydata[0:6],1,cov=True)
err_L = np.sqrt(np.diag(cov_par_min))
links_m = ufloat(par_min[0],err_L[0])
links_c = ufloat(par_min[1],err_L[1])
par_const,cov_const = np.polyfit(xdata[6:18],ydata[6:18],0,cov=True)
err_C = np.sqrt(np.diag(cov_const))
plat = ufloat(par_const,err_C)
par_pos,cov_pos = np.polyfit(xdata[17:22],ydata[17:22],1,cov=True)
err_R = np.sqrt(np.diag(cov_pos))
rechts_m = ufloat(par_pos[0],err_R[0])
rechts_c = ufloat(par_pos[1],err_R[1])

plt.plot(xplot1,par_min[0]*xplot1 + par_min[1],label = "lin. Fit an die Linke Flanke ")
plt.plot(xplot2,par_const + xplot2*0,label = "konstanter Fit an das Plateau ")
plt.plot(xplot3,par_pos[0]*xplot3 + par_pos[1],label = "lin. Fit an die Rechte Flanke  ")
print("plateau:",par_const ,"+-", err_C  )
print("Links:",par_min ,"+-", err_L  )
print("Rechts:",par_pos ,"+-", err_R  )
print("t_Links:", (plat/2 - links_c)/links_m )
print("t_rechts:", (plat/2 - rechts_c)/rechts_m )
print("t_1/2 : ",(plat/2 - rechts_c)/rechts_m -(plat/2 - links_c)/links_m  )
plt.legend(loc = "best")
plt.savefig("plot.pdf")
plt.close()











start = 3604339
messzeit = 181285
suchzeit = 1e-5
lam = start*suchzeit/messzeit
print("lambda:",lam)
U = lam*np.exp(-lam)*start 
print("gesamt",U)
print("cahnnel",U/511)

for i in range(57):
    print(i,"&",np.round(p(i),1),"&",np.round(N[i],1),"&",i+57,"&",np.round(p(i+57),1),"&",np.round(N[i+57],1),"&",i+57+57,"&",np.round(p(i+57+57),1),"&",np.round(N[i+57+57],1),"&",i+171,"&",np.round(p(i+171),1),"&",np.round(N[i+171],1), " \\\\ " )
