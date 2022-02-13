import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
import matplotlib.pyplot as plt
import scipy.constants as const

import csv

gamma = const.value('proton gyromag. ratio') #[1/sT], gyromagnetisches Verhältnis
print(gamma)