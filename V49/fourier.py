import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import matplotlib.pyplot as plt

import csv

#----------Messwerte auslesen----------
data = np.genfromtxt('Messwerte/fourier.csv', delimiter=',', unpack=True)
times = data[0]
real = data[1]
imag = data[2]

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
plt.figure()
plt.hlines(0, -8.825151682294324928, 6.618863761720744151, linestyles='dashed', color='black', linewidth=1,  label=r'Durchmesser $d_f$')
plt.plot(freqs[(freqs>-15000) & (freqs<15000)]*1e-3, np.real(fftdata)[(freqs>-15000) & (freqs<15000)], 'x', label="Fourier-Transformation")
plt.xlabel(r'$f / \, \mathrm{kHz}$')
plt.ylabel(r'Anzahl Protonen')
plt.legend(loc='upper right', fontsize = 'small')
plt.grid()
plt.savefig("figures/fourier_trafo.pdf")
plt.clf()

d_f = (8.825151682294324928 + 6.618863761720744151)*1e3 #[Hz]
d_Probe = 4.2*1e-3 #[m]
gamma = 267.5 * 1e6 #[1/sT], gyromagnetisches Verhältnis
g = 2 * np.pi * d_f /(d_Probe * gamma)
print('Durchmesser: ', d_f, ' Hz')
print('Gradientenstärke: ', g, ' T/m')