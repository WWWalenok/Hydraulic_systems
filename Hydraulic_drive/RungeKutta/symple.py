
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import math

warnings.filterwarnings(action='once')


out_name = 'charracter'

large = 20; med = 14; small = 12
params = {'axes.titlesize': large,
          'legend.fontsize': med,
          'figure.figsize': (12, 6),
          'axes.labelsize': med,
          'axes.titlesize': med,
          'xtick.labelsize': med,
          'ytick.labelsize': med,
          'figure.titlesize': large}


plt.rcParams.update(params)
plt.style.use('seaborn-whitegrid')
sns.set_style("white")

# Import Data

dt = 0.0000001

def D2S(D):
	return ((D * 0.5)**2) * 3.1415926535

V1 = 0.01
S1 = D2S(100 * 1E-3)

V2 = 0.01
S2 = D2S(100 * 1E-3) - D2S(50 * 1E-3)


ro = 860
E = 1.56E9

mu = 0.62
S = D2S((14 - 6) * 1E-3)

min_x = 0
max_x = 630 * 1E-3
m = 1
b = 5
p_sliv = 134
p_input = 20 * 1E6

F = 0

U0=1

S0=S*U0

F_0=p_input*S1 - p_sliv*S2 + F

K = (ro * (S1 ** 3 + S2 ** 3)) / (2 * mu ** 2)

print(K)

v_r = (S0**2) / (2 * K) * ((b**2 + 4 * K / (S0**2) * abs(F_0)) ** 0.5 - b)

v_r_max = (S**2) / (2 * K) * ((b**2 + 4 * K / (S**2) * abs(F_0)) ** 0.5 - b)

print(v_r, v_r_max)

p1_r = p_input - v_r * abs(v_r) * ro * S1**2 / (S0**2 * 2 * mu**2)
p2_r = p_sliv + v_r * abs(v_r) * ro * S2**2 / (S0**2 * 2 * mu**2)

# Plot Line1 (Left Y Axis)

arr_F = []
arr_V_F_04 = []
arr_V_F_07 = []
arr_V_F_1 = []

arr_S = []
arr_V_S_04 = []
arr_V_S_07 = []
arr_V_S_1 = []



for i in range(1, 1000):
    TF = F_0 * i/ 1000
    arr_F.append(i / 1000)
    TS = 0.4 * S0
    arr_V_F_04.append((TS ** 2) / (2 * K) * ((b ** 2 + 4 * K / (TS ** 2) * abs(TF)) ** 0.5 - b))

    TS = 0.7 * S0
    arr_V_F_07.append((TS ** 2) / (2 * K) * ((b ** 2 + 4 * K / (TS ** 2) * abs(TF)) ** 0.5 - b))

    TS = 1 * S0
    arr_V_F_1.append((TS ** 2) / (2 * K) * ((b ** 2 + 4 * K / (TS ** 2) * abs(TF)) ** 0.5 - b))

for i in range(1, 1000):
    TS = S * i / 1000
    arr_S.append(i / 1000)
    TF = 0.4 * F_0
    arr_V_S_04.append((TS ** 2) / (2 * K) * ((b ** 2 + 4 * K / (TS ** 2) * abs(TF)) ** 0.5 - b))

    TF = 0.7 * F_0
    arr_V_S_07.append((TS ** 2) / (2 * K) * ((b ** 2 + 4 * K / (TS ** 2) * abs(TF)) ** 0.5 - b))

    TF = 1 * F_0
    arr_V_S_1.append((TS ** 2) / (2 * K) * ((b ** 2 + 4 * K / (TS ** 2) * abs(TF)) ** 0.5 - b))

ax1 = plt.subplot(121)
ax2 = plt.subplot(122)

V_F_04_line, = ax1.plot(arr_F, arr_V_F_04)
V_F_07_line, = ax1.plot(arr_F, arr_V_F_07)
V_F_1_line, = ax1.plot(arr_F, arr_V_F_1)

V_S_04_line, = ax2.plot(arr_S, arr_V_S_04)
V_S_07_line, = ax2.plot(arr_S, arr_V_S_07)
V_S_1_line, = ax2.plot(arr_S, arr_V_S_1)


ax1.set_xlabel('F / F0', fontsize=16)
ax1.tick_params(axis='x', rotation=0, labelsize=12)
ax1.set_ylabel('Скорость равновесия', fontsize=16)
ax1.tick_params(axis='y', rotation=0 )
ax1.grid(alpha=.2)

ax2.set_xlabel('S / S0', fontsize=16)
ax2.tick_params(axis='x', rotation=0, labelsize=12)
ax2.set_ylabel('Скорость равновесия', fontsize=16)
ax2.tick_params(axis='y', rotation=0 )
ax2.grid(alpha=.2)

ax1.legend(handles = [V_F_04_line, V_F_07_line, V_F_1_line], labels = ['S = S0 * 0.4', 'S = S0 * 0.7', 'S = S0'])
ax2.legend(handles = [V_S_04_line, V_S_07_line, V_S_1_line], labels = ['F = F0 * 0.4', 'F = F0 * 0.7', 'F = F0'])

plt.savefig("E:\\Repos\\Latex\\rsgc_1\\" + out_name + ".png", format = "png")

plt.show()



