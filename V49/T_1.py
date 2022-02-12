import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import matplotlib.pyplot as plt

import csv

from matplotlib import rc
plt.rc('font', family='serif')
rc('text', usetex=True)

#----------Messwerte auslesen----------
t_1 = []
U_1 = []

with open('Messwerte/T_1.csv', 'r') as file:
    data = csv.reader(file, delimiter = ',', quoting = csv.QUOTE_NONNUMERIC)
    for row in data:
        t_1.append(row[0])
        U_1.append(row[1])

#----------Fitfunktion für T_1----------
def exp(t, T_1, M_0):
    return M_0 * (1 - 2*np.exp(-t/T_1))

#----------Bestimmung von T_1----------
parameters, covariance_matrix = curve_fit(exp, t_1, U_1)
uncertainties = np.sqrt(np.diag(covariance_matrix))

T_1 = ufloat(parameters[0], uncertainties[0])
M_0 = ufloat(parameters[1], uncertainties[1])

print('Parameter Messung T_1:', '\n', 'M_0 = ', M_0, ' V', '\n', 'T_1 = ', T_1, ' s', '\n')

t_values = np.linspace(t_1[0]-0.01, t_1[-1]+1, 1000)

plt.figure()
plt.plot(t_1, U_1, 'x', label = r'Messwerte')
plt.plot(t_values, exp(t_values, *parameters), 'g-', label = r'Fitfunktion')
plt.xscale('log')
plt.xlabel(r'$\tau / \, \mathrm{s}$')
plt.ylabel(r'$U / \, \mathrm{V}$')
plt.legend(loc = 'best')
plt.grid()
plt.savefig('figures/T_1.pdf')
plt.clf()

#----------Fitfunktion für T_2----------

#----------Bestimmung von T_2----------



#----------Diffusionsmessung----------

#----------Diskussion----------