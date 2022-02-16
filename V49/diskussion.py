import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import matplotlib.pyplot as plt
import scipy.constants as const

import csv

gamma = const.value('proton gyromag. ratio') #[1/sT], gyromagnetisches Verhältnis
print(gamma)

T_1_lit = ufloat(3.09, 0.15)
T_2_lit = ufloat(1.52, 0.093)
T_lit = T_1_lit/T_2_lit
print(T_lit)


T_1_exp = ufloat(2.39, 0.08)
T_2_exp = ufloat(1.37, 0.07)
T_exp = T_1_exp/T_2_exp
print(T_exp)
print('{:.4f}'.format(T_2_exp**3))

abweichung = (T_exp - 1.85)/1.85 * 100
print('{:.4f}'.format(abweichung))

D = ufloat(1.71, 0.04)*1e-9
g = 0.0864
print('{:.10f}'.format(3 /2 /D /gamma**2 /g**2))
print('{:.10f}'.format(T_2_exp**3 /3 *2 *D *gamma**2 *g**2))
