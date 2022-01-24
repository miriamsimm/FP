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
plt.plot(xplot,Gauss(xplot, *par_koin),color = "blue",label = "Fitkurve ")
plt.xlabel(r"Verzögerungszeit $T_V$ / ns")
plt.ylabel("Counts")
plt.xlim(-20,20)
plt.ylim(0,250)
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
