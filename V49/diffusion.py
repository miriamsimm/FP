import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import matplotlib.pyplot as plt
import scipy.constants as const

import csv

#----------Messwerte auslesen----------
times, real, imag = np.genfromtxt('Messwerte/fourier.csv', delimiter=',', unpack=True)
t_echo, U_echo = np.genfromtxt('Messwerte/echo.csv', delimiter=',', unpack=True)

T_2 = 1.37

#----------Fitfunktion Echo----------
def echo(t, D, M_0, M_1):
    return M_0 * np.exp(-2*t/T_2) * np.exp(-t**3/D) + M_1

#----------Bestimmung der Diffusionskonstante----------
parameters, covariance_matrix = curve_fit(echo, t_echo, U_echo, p0=[1e-6,1,1e-3])
uncertainties = np.sqrt(np.diag(covariance_matrix))

T_D = ufloat(parameters[0], uncertainties[0])
M_0 = ufloat(parameters[1], uncertainties[1])
M_1 = ufloat(parameters[2], uncertainties[2])

#print('Parameter Messung Diffusion:', '\n', 'T_D = ', T_D, ' s', '\n', 'M_0 = ', M_0, ' V', '\n', 'M_1 = ', M_1, ' V', '\n', )

file = open('Messwerte/Ergebnisse.txt', 'a')
file.write('Parameter Messung Diffusion: \n')
input_M_O = '{} {} {} \n'.format('M_0 = ', M_0, ' V')
input_M_1 = '{} {} {} \n'.format('M_1 = ', M_1, ' V')
input_T_D = '{} {} {} \n'.format('T_D = ', T_D, ' s')
blankspace = '{} \n'.format(' ')
file.write(input_M_O)
file.write(input_M_1)
file.write(input_T_D)
file.write(blankspace)
file.close()

t_values = np.linspace(t_echo[0]-1e-5, t_echo[-1]+1e-4, 1000)
t_values_log = np.linspace(t_echo[0]-1e-5, t_echo[-1]+1e-4, 1000)

plt.figure()
plt.subplot(2, 1, 1)
plt.title(r'Logarithmische Darstellung')
plt.xlabel(r'$\tau^3 / \, \si{\micro\second}$')
plt.ylabel(r'$\ln (U(\tau)) -2\tau / \, T_2$')
plt.plot(t_echo**3*1e6, np.log(U_echo)-2*t_echo/T_2, 'gx', label = r'Messwerte')
plt.plot(t_values_log**3*1e6, np.log(echo(t_values_log, *parameters))-2*t_values_log/T_2, '-', label = r'Fitfunktion')
plt.xscale('log')
plt.legend(loc = 'best')
plt.grid()

plt.subplot(2, 1, 2)
plt.title(r'Standard-Darstellung')
plt.xlabel(r'$\tau^3 / \, \si{\micro\second}$')
plt.ylabel(r'$\ln (U(\tau)) -2\tau / \, T_2$')
plt.plot(t_echo**3*1e6, np.log(U_echo)-2*t_echo/T_2, 'gx', label = r'Messwerte')
plt.plot(t_values**3*1e6, np.log(echo(t_values, *parameters))-2*t_values/T_2, '-', label = r'Fitfunktion')
plt.legend(loc = 'best')
plt.grid()

plt.tight_layout()
plt.savefig('figures/echo.pdf')
plt.clf()

#----------Plot der Messwerte für die Fouriertransformation----------
#Einzeichnen der Linie bei der das Echo-Minimum abgeschnitten wird
start = np.argmin(real) #argmin statt argmax da Minimum statt Maximum des Realteils

plt.figure()
plt.subplot(2, 1, 1)
plt.title(r'Messwerte ohne Phasenkorrektur')
plt.vlines(times[start]*1e3, 0, real[start], linestyles='dashed', color='black', linewidth=1,  label=r'Beginn des Spektrums')
plt.plot(times*1e3, real, 'b-', label = r'Realteil')
plt.plot(times*1e3, imag, 'r-', label = r'Imaginärteil')
plt.xlabel(r'$\tau / \, \mathrm{ms}$')
plt.ylabel(r'$U / \, \mathrm{V}$')
plt.legend(loc = 'best', fontsize = 'small')
plt.grid()

plt.subplot(2, 1, 2)
plt.title(r'Messwerte mit Phasenkorrektur')
phase = np.arctan2(imag[start], real[start])
real = real*np.cos(phase)+imag*np.sin(phase)*(-1)
imag = -real*np.sin(phase)+imag*np.cos(phase)*(-1)
plt.vlines(times[start]*1e3, 0, real[start], linestyles='dashed', color='black', linewidth=1,  label=r'Beginn des Spektrums')
plt.plot(times*1e3, real, 'b-', label = r'Realteil')
plt.plot(times*1e3, imag, 'r-', label = r'Imaginärteil')
plt.xlabel(r'$\tau / \, \mathrm{ms}$')
plt.ylabel(r'$U / \, \mathrm{V}$')
plt.legend(loc = 'best', fontsize = 'small')
plt.grid()

plt.tight_layout()
plt.savefig('figures/fourier_plot.pdf')
plt.clf()

#----------Fouriertransformation----------
#Suchen des Echo-Minimums und alle Daten davor abschneiden
times = times[start:]
real = real[start:]
imag = imag[start:]
#Phasenkorrektur - der Imaginärteil bei t=0 muss 0 sein
#phase = np.arctan2(imag[0], real[0])
#Daten in komplexes Array mit Phasenkorrektur speichern
compsignal = real + imag*1j
#compsignal = (real*np.cos(phase)+imag*np.sin(phase))+ (-real*np.sin(phase)+imag*np.cos(phase))*1j
#Offsetkorrektur, ziehe den Mittelwert der letzten 512 Punkte von allen Punkten ab
compsignal = compsignal - compsignal[-512:-1].mean()
#Der erste Punkt einer FFT muss halbiert werden
compsignal[0] = compsignal[0]/2.0
#Anwenden einer Fensterfunktion (siehe z. Bsp.
#https://de.wikipedia.org/wiki/Fensterfunktion )
#Hier wird eine Gaußfunktion mit sigma = 100 Hz verwendet
apodisation = 100.0*2*np.pi
compsignal = compsignal*np.exp(-1.0/2.0*((times-times[0])*apodisation)**2)
#Durchführen der Fourier-Transformation
fftdata = np.fft.fftshift(np.fft.fft(compsignal))
#Generieren der Frequenzachse
freqs = np.fft.fftshift(np.fft.fftfreq(len(compsignal), times[1]-times[0]))
#Speichern des Ergebnisses als txt-Datei
np.savetxt('Messwerte/fourier_echo_gradient.txt', np.array([freqs, np.real(fftdata), np.imag(fftdata)]).transpose())
#Erstellen eines Plots
plt.figure(figsize=(6.4,4.0))
plt.hlines(0, -8.825151682294324928, 6.618863761720744151, linestyles='dashed', color='black', linewidth=1,  label=r'Durchmesser $d_f$')
plt.plot(freqs[(freqs>-15000) & (freqs<15000)]*1e-3, np.real(fftdata)[(freqs>-15000) & (freqs<15000)], 'x', label="Fourier-Transformation")
plt.xlabel(r'$f / \, \mathrm{kHz}$')
plt.ylabel(r'Anzahl Protonen')
plt.legend(loc='upper right', fontsize = 'small')
plt.grid()
plt.savefig("figures/fourier_trafo.pdf")
plt.clf()

d_f = (8.825151682294324928 + 6.618863761720744151)*1e3 #[Hz]
d_Probe = 4.2e-3 #[m]
gamma = const.value('proton gyromag. ratio') #[1/sT], gyromagnetisches Verhältnis
g = 2 * np.pi * d_f /d_Probe /gamma #[T/m]
tau = 0.2e-3 #[s]
D = 3/2 /T_D /gamma**2 /g**2 #[m^2/s]
#print('Durchmesser: ', d_f, ' Hz')
#print('Gradientenstärke: ', g, ' T/m')
#print('Diffusionskonstante: ', D, 'm^2/s')

file = open('Messwerte/Ergebnisse.txt', 'a')
file.write('Parameter Messung Diffusion: \n')
input_Durchmesser = '{} {} {} \n'.format('Durchmesser: ', d_f, ' Hz')
input_Gradient = '{} {} {} \n'.format('Gradientenstärke: ', g, ' T/m')
input_Diffusionskonstante = '{} {} {} \n'.format('Diffusionskonstante: ', D, ' s')
blankspace = '{} \n'.format(' ')
file.write(input_Durchmesser)
file.write(input_Gradient)
file.write(input_Diffusionskonstante)
file.write(blankspace)
file.close()
