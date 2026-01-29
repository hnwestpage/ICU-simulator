# Sensitivity Plots

import numpy as np
import math
import matplotlib.pyplot as plt

lmbda_range = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]   # Begen et al., (2024) {URBAN} Murthy et al. (2015) {low income/rural}
mu_1_range =  [2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]  # Moitra et al., (2017)
mu_2_range =  [1/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10, 1/11, 1/12, 1/13, 1/14, 1/15] # Essafi et al., 2022
rho_1_range =  [0.1, 0.08, 0.06, 0.04, 0.02, 0.01]            # Lo, (2001)
eta_range = [0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011]  #Lai et al., 2021
nu_range = [1/7, 1/8, 1/9, 1/10]                              #ACEP n.d. (accessed 2025)
ratio_range =  [1/8, 1/9, 1/10, 1/11]                         # Bhatla & Ryskina (2020)

lmbda_results = [150.76, 179.84, 213.8, 276.12, 333.44, 432.36, 559.44, 761.88, 1050.92, 1312.32, 1597.4, 1891.4]
mu_1_results = [273.4, 321.32,	411.28,	509.4,	736.6,	1116.88, 1376.44, 1586.88, 1732.56, 1848.84, 1925.16, 1996.52, 2077.44, 2128]


fig, ax = plt.subplots(nrows=2,ncols=2, sharex=True, layout='constrained')
ax[0,0].plot(lmbda_range, lmbda_results)
ax[0,1].plot(mu_1_range, mu_1_results)

plt.show()