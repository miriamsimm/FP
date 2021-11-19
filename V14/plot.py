import numpy as np
from scipy.optimize import curve_fit
import scipy.constants as const
import matplotlib.pyplot as plt
from matplotlib import rc
from uncertainties import ufloat
import uncertainties.unumpy as unp
import math
rc('text', usetex=True)

#--------------------Daten auslesen--------------------

Leermessung = np.genfromtxt('data/leermessung.Spe', unpack=True)/300

#Ausgangsintensitäten 
I_W3 = np.zeros(12)
I_W3[[0, 1, 2, 3, 4, 5]] = 1301/241.24
I_W3[[7, 10]] = 594/241.24
I_W3[[6, 8, 9, 11]] = 1164/241.24

#I_W3 = np.zeros(3)
#I_W3[[0]] = 1301/241.24
#I_W3[[2]] = 594/241.24
#I_W3[[1]] = 1164/241.24

I_W4 = np.zeros(12)
I_W4[[0]] = 30557/241.80
I_W4[[1]] = 1612/241.46
I_W4[[2]] = 30303/241.78
I_W4[[3]] = 10235/241.46
I_W4[[4]] = 9544/241.46
I_W4[[5]] = 10361/241.46
I_W4[[6]] = 6993/241.46
I_W4[[7]] = 6196/241.46
I_W4[[8]] = 9419/241.46
I_W4[[9]] = 7478/241.46
I_W4[[10]] = 6377/241.46
I_W4[[11]] = 9861/241.46

#Fehler 
I_W3 = unp.uarray(I_W3, np.sqrt(I_W3))
I_W4 = unp.uarray(I_W4, np.sqrt(I_W4))

#--------------------Matrix A--------------------#

s = np.sqrt(2)
A = np.matrix([[1,1,1,0,0,0,0,0,0],  #1
               [0,0,0,1,1,1,0,0,0],  #2
               [0,0,0,0,0,0,1,1,1],  #3
               [1,0,0,1,0,0,1,0,0],  #4
               [0,1,0,0,1,0,0,1,0],  #5
               [0,0,1,0,0,1,0,0,1],  #6
               [0,s,0,0,0,s,0,0,0],  #7
               [s,0,0,0,s,0,0,0,s],  #8
               [0,0,0,s,0,0,0,s,0],  #9
               [0,s,0,s,0,0,0,0,0],  #10
               [0,0,s,0,s,0,s,0,0],  #11
               [0,0,0,0,0,s,0,s,0]]) #12

#A_2 = np.matrix([[3, 0, 0], [0, 2*s, 0], [0, 0, 3*s]])


#--------------------Funktionen definieren--------------------#

def Absorptionskoeffizient(A, I):
    A_T = np.transpose(A)
    A_I = np.linalg.inv(A_T @ A)
    A_P = A_I @ A_T
    mu = np.dot(A_P, I) #d = 1cm implizit
    mu = np.transpose(mu) #Ausgabe
    return np.array(mu)
    #mu = np.linalg.lstsq(A,I)

def Intensitäten(I_W, I_Aluminium):
    Intensitäten = unp.log(I_Aluminium/I_W)
    return np.array(Intensitäten)

#------------------------------------------------------------------------------------

#--------------------Leermessung und Plot des Absorptionsspektrums--------------------

Kanal = np.linspace(0, len(Leermessung)-1, len(Leermessung))

plt.plot(Kanal[0:200], Leermessung[0:200], linewidth = 1, label = 'Messwerte')
plt.xlabel(r'$\mathrm{Kanalnummer}$')
plt.ylabel(r'$\textrm{Zählrate}\, / \, s^{-1}$')
plt.xlim(20, 200)
plt.grid()
plt.legend(loc="best")
plt.savefig('Leermessung.pdf')
plt.clf()

#--------------------Würfel 1: Aluminium--------------------#

C_Aluminium = np.ones(12)*37660
t_Aluminium = 241.69 #s
I_Aluminium = C_Aluminium/t_Aluminium #Eingangsintensität
#I_Aluminium = np.array(I_Aluminium)
I_Aluminium = unp.uarray(I_Aluminium, np.sqrt(I_Aluminium))

#C_Aluminium_2 = np.ones(3)*37660
#t_Aluminium_2 = 241.69 #s
#I_Aluminium_2 = C_Aluminium_2/t_Aluminium_2 #Eingangsintensität
#I_Aluminium_2 = unp.uarray(I_Aluminium_2, np.sqrt(I_Aluminium_2))

#--------------------Würfel 3: Schweres Material--------------------#

Int_W3 = Intensitäten(I_W3, I_Aluminium)
AbK_W3 = Absorptionskoeffizient(A, Int_W3)
print(AbK_W3)
#print(AbK_W3[0])

#--------------------Würfel 4: Unbekanntes Material--------------------#

Int_W4 = Intensitäten(I_W4, I_Aluminium)
AbK_W4 = Absorptionskoeffizient(A, Int_W4)
print(AbK_W4)
#print(AbK_W4[0])