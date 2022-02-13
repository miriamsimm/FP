import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import matplotlib.pyplot as plt
import scipy.constants as const

import csv

#----------Messwerte auslesen----------
t_1, U_1 = np.genfromtxt('Messwerte/T_1.csv', delimiter=',', unpack=True)

#----------Fitfunktion für T_1----------
def exp(t, T_1, M_0):
    return M_0 * (1 - 2*np.exp(-t/T_1))

#----------Bestimmung von T_1----------
parameters, covariance_matrix = curve_fit(exp, t_1, U_1)
uncertainties = np.sqrt(np.diag(covariance_matrix))

T_1 = ufloat(parameters[0], uncertainties[0])
M_0 = ufloat(parameters[1], uncertainties[1])

#print('Parameter Messung T_1:', '\n', 'M_0 = ', M_0, ' V', '\n', 'T_1 = ', T_1, ' s', '\n')

file = open('Messwerte/Ergebnisse.txt', 'a')
file.write('Parameter Messung T_1: \n')
input_M_O = '{} {} {} \n'.format('M_0 = ', M_0, ' V')
input_T_1 = '{} {} {} \n'.format('T_1 = ', T_1, ' s')
blankspace = '{} \n'.format(' ')
file.write(input_M_O)
file.write(input_T_1)
file.write(blankspace)
file.close()

t_values = np.linspace(t_1[0]-0.01, t_1[-1]+1, 1000)

plt.figure(figsize=(6.4,4.0))
plt.plot(t_1, U_1, 'x', label = r'Messwerte')
plt.plot(t_values, exp(t_values, *parameters), 'g-', label = r'Ausgleichsrechnung')
plt.xscale('log')
plt.xlabel(r'$\tau / \, \mathrm{s}$')
plt.ylabel(r'$U / \, \mathrm{V}$')
plt.legend(loc = 'best')
plt.grid()
plt.savefig('figures/T_1.pdf')
plt.clf()
