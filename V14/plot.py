import numpy as np
import uncertainties
from scipy.optimize import curve_fit
import scipy.constants as const
import matplotlib.pyplot as plt
from matplotlib import rc
from uncertainties import ufloat
import uncertainties.unumpy as unp
from uncertainties.unumpy import std_devs as stds
from uncertainties.unumpy import nominal_values as noms
import math
rc('text', usetex=True)

#--------------------Daten auslesen--------------------#


Leermessung = np.genfromtxt('data/leermessung.Spe', unpack=True)/300 #wird nicht gebraucht

#Ausgangsintensitäten 
N_W3 = np.array([1301, 594, 1164])
T_W3 = np.array([241.24, 241.24, 241.24])

I_W3 = np.zeros(12)
I_W3[[0, 1, 2, 3, 4, 5]] = N_W3[0]/T_W3[0]   
I_W3[[7, 10]] = N_W3[1]/T_W3[1]              
I_W3[[6, 8, 9, 11]] = N_W3[2]/T_W3[2]        
I_W3 = unp.uarray(I_W3, np.sqrt(I_W3))

print(f"\n Intensitäten Würfel 3: \n {I_W3}")

N_W4 = np.array([30557, 1612, 30303, 10235, 9544, 10361, 6993, 6196, 9419, 7478, 6377, 9861]) 
T_W4 = np.array([241.80, 241.42, 241.78, 241.46, 241.46, 241.48, 241.48, 241.46, 241.56, 241.50, 241.50, 241.54]) 

I_W4 = N_W4/T_W4
I_W4 = unp.uarray(I_W4, np.sqrt(I_W4))

print(f"\n Intensitäten Würfel 4: \n {I_W4}")


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

A_T = np.transpose(A)

#--------------------Funktionen definieren--------------------#

def Absorptionskoeffizient(V_I, V_Abk, I):
    V_I_Inv = np.linalg.inv(V_I)
    Mat_Prod = (V_Abk * A_T) * V_I_Inv
    mu = np.dot(Mat_Prod, I)
    mu = np.transpose(mu)
    return np.array(mu)

def Intensitäten(I_W, I_Aluminium):
    Intensitäten = unp.log(I_Aluminium/I_W)
    return np.array(Intensitäten)

def Kovarianz_Int(I):
    V = np.zeros((len(I), len(I)))
    np.fill_diagonal(V, I**2)
    return np.matrix(V)

def Kovarianz_Abk(V):
    V_Inv = np.linalg.inv(V)
    Mat_Prod = (A_T * (V_Inv * A))
    return np.linalg.inv(Mat_Prod)

#-------------------------------------------------------------------------------------#

#--------------------Leermessung und Plot des Absorptionsspektrums--------------------#

Kanal = np.linspace(0, len(Leermessung)-1, len(Leermessung))

plt.plot(Kanal[0:200], Leermessung[0:200], linewidth = 1, label = 'Messwerte')
plt.xlabel(r'$\mathrm{Kanalnummer}$')
plt.ylabel(r'$\textrm{Zählrate}\, / \, s^{-1}$')
plt.axvline(128, linewidth=1, ymin=0, ymax=0.953, color='g', label = '662 keV')
#plt.plot(128, 19, marker='o', markersize=3, color = 'g')
plt.xlim(20, 200)
plt.grid()
plt.legend(loc="best")
plt.savefig('Leermessung.pdf')
plt.clf()

#--------------------Würfel 1: Aluminium--------------------#

C_Aluminium = np.ones(12)*37660
t_Aluminium = 241.69 
I_Aluminium = C_Aluminium/t_Aluminium #Eingangsintensität
I_Aluminium = unp.uarray(I_Aluminium, np.sqrt(I_Aluminium))

#--------------------Würfel 3: Schweres Material--------------------#

Int_W3 = Intensitäten(I_W3, I_Aluminium)

#Messunsicherheiten und Varianzen
V_I_W3 = Kovarianz_Int(stds(Int_W3)) #Matrix mit dem Quadrat der Standardabweichung = Varianz auf der Diagonalen
V_Abk_W3 = Kovarianz_Abk(V_I_W3)
Abk_W3 = Absorptionskoeffizient(V_I_W3, V_Abk_W3, noms(Int_W3)) #mu
Abk_W3_Err = np.array(np.sqrt(np.diag(V_Abk_W3)))

print(f"\n Absorptionskoeffizienten Würfel 3: \n {Abk_W3} \n {Abk_W3_Err}")

#--------------------Würfel 4: Unbekanntes Material--------------------#

Int_W4 = Intensitäten(I_W4, I_Aluminium)

#Messunsicherheiten und Varianzen
V_I_W4 = Kovarianz_Int(stds(Int_W4)) #Matrix mit dem Quadrat der Standardabweichung = Varianz auf der Diagonalen
V_Abk_W4 = Kovarianz_Abk(V_I_W4)
Abk_W4 = Absorptionskoeffizient(V_I_W4, V_Abk_W4, noms(Int_W4)) #mu
Abk_W4_Err = np.array(np.sqrt(np.diag(V_Abk_W4)))

print(f"\n Absorptionskoeffizienten Würfel 4: \n {Abk_W4} \n {Abk_W4_Err}")

#--------------------Vergleich mit Literaturwerten--------------------#

Eisen = 0.606
Aluminium = 0.211
Blei = 1.419
Messing = 0.683
Delrin = 0.121

print(f"\n")
print(f"Abweichung zu Eisen: \n {Abk_W4 - Eisen} \n")
print(f"Abweichung zu Aluminium: \n {Abk_W4 - Aluminium} \n")
print(f"Abweichung zu Blei: \n {Abk_W4 - Blei} \n")
print(f"Abweichung zu Messing: \n {Abk_W4 - Messing} \n")
print(f"Abweichung zu Delrin: \n {Abk_W4 - Delrin} \n")