import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from uncertainties import ufloat
import matplotlib.pyplot as plt

import csv

from matplotlib import rc
plt.rc('font', family='serif')
rc('text', usetex=True)

#----------Messwerte auslesen----------
#t_2 = []
#U_2 = []

#with open('Messwerte/T_2.csv', 'r') as file:
#    data = csv.reader(file, delimiter = ',', quoting = csv.QUOTE_NONNUMERIC)
#    for row in data:
#        t_2.append(row[0])
#        U_2.append(row[1])

n, y = np.genfromtxt('Messwerte/T_2.csv', delimiter=',', unpack=True)

t_2 = []
for k in n:
    t_2.append(2*k*0.022)

U_2 = []
for p in y:
    U_2.append((300-p)*0.005)

#----------Fitfunktion für T_2----------
def exp(t, T_2, M_0, M_1):
    return M_0 * np.exp(-t/T_2) + M_1

#----------Bestimmung von T_2----------
parameters, covariance_matrix = curve_fit(exp, t_2, U_2)
uncertainties = np.sqrt(np.diag(covariance_matrix))

T_2 = ufloat(parameters[0], uncertainties[0])
M_0 = ufloat(parameters[1], uncertainties[1])
M_1 = ufloat(parameters[2], uncertainties[2])

print('Parameter Messung T_2:', '\n', 'T_2 = ', T_2, ' s', '\n', 'M_0 = ', M_0, ' V', '\n', 'M_1 = ', M_1, ' V', '\n', )

t_values = np.linspace(t_2[0]-0.05, t_2[-1]+0.05, 1000)

plt.figure()
plt.plot(t_2, U_2, 'x', color='darkgreen', label = r'Messwerte')
for i in range(len(U_2)):
    plt.vlines(2*2*(i+1)*0.022, 0, U_2[i], linestyles='solid')
plt.plot(t_values, exp(t_values, *parameters), 'g-', label = r'Fitfunktion')
plt.xlabel(r'$\tau / \, \mathrm{s}$')
plt.ylabel(r'$U / \, \mathrm{V}$')
plt.legend(loc = 'best')
plt.grid()
plt.savefig('figures/T_2.pdf')
plt.clf()
