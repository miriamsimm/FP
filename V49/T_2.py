import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import matplotlib.pyplot as plt
import scipy.constants as const

import csv

#----------Messwerte auslesen----------
n, y = np.genfromtxt('Messwerte/T_2.csv', delimiter=',', unpack=True)

#Umrechnen der ausgemessenen Werte
t_2 = []
for k in n:
    t_2.append(2*k*0.022)

U_2 = []
for p in y:
    U_2.append((150-p)*10*1e-3)

file = open('Messwerte/T_2_Werte.csv', 'w')
write = csv.writer(file)
for n in range (len(t_2)):
    input = [t_2[n], U_2[n]]
    write.writerow(input)


#----------Fitfunktion für T_2----------
def exp(t, T_2, M_0, M_1):
    return M_0 * np.exp(-t/T_2) + M_1

#----------Bestimmung von T_2----------
parameters, covariance_matrix = curve_fit(exp, t_2, U_2)
uncertainties = np.sqrt(np.diag(covariance_matrix))

T_2 = ufloat(parameters[0], uncertainties[0])
M_0 = ufloat(parameters[1], uncertainties[1])
M_1 = ufloat(parameters[2], uncertainties[2])

#print('Parameter Messung T_2:', '\n', 'T_2 = ', T_2, ' s', '\n', 'M_0 = ', M_0, ' V', '\n', 'M_1 = ', M_1, ' V', '\n', )

file = open('Messwerte/Ergebnisse.txt', 'a')
file.write('Parameter Messung T_2: \n')
input_M_O = '{} {} {} \n'.format('M_0 = ', M_0, ' V')
input_M_1 = '{} {} {} \n'.format('M_1 = ', M_1, ' V')
input_T_2 = '{} {} {} \n'.format('T_2 = ', T_2, ' s')
blankspace = '{} \n'.format(' ')
file.write(input_M_O)
file.write(input_M_1)
file.write(input_T_2)
file.write(blankspace)
file.close()


t_values = np.linspace(t_2[0]-0.05, t_2[-1]+0.05, 1000)

plt.figure(figsize=(6.4,4.0))
plt.plot(t_values, exp(t_values, *parameters), '-', color='darkorange', label = r'Fitfunktion')
plt.plot(t_2, U_2, 'x', color='green', label = r'Messwerte')
for i in range(len(U_2)):
    plt.vlines(2*2*(i+1)*0.022, 0, U_2[i], linestyles='solid')
plt.xlabel(r'$\tau / \, \mathrm{s}$')
plt.ylabel(r'$U / \, \mathrm{V}$')
plt.legend(loc = 'best')
plt.grid()
plt.savefig('figures/T_2.pdf')
plt.tight_layout()
plt.clf()
