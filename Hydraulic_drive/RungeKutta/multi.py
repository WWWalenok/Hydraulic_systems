print("Hi")


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import math

warnings.filterwarnings(action='once')
datas = []
for i in range(4):
    data = np.genfromtxt(
        "out_{a}.csv".format(a = i+1),
        delimiter=',',
        skip_header=1,
        skip_footer=10,
        names=
        [
            't',
            'v',
            'x',
            'p1',
            'p2',
            'f',
            'h',
            'n',
            'a',
            'cu',
            'y',
            'dt'
        ]
    )
    datas.append(data)

large = 20; med = 14; small = 12
params = {'axes.titlesize': large,
          'legend.fontsize': med,
          'figure.figsize': (12, 8),
          'axes.labelsize': med,
          'axes.titlesize': med,
          'xtick.labelsize': med,
          'ytick.labelsize': med,
          'figure.titlesize': large}


plt.rcParams.update(params)
plt.style.use('seaborn-whitegrid')
sns.set_style("white")

# Import Data

X_tag = 't'
Y1_tags = 'v'
Y2_tags = 'y'

# Plot Line1 (Left Y Axis)
fig, ax1 = plt.subplots(1, 1, figsize=(12,8))


I = 0

tags = []

for data in datas :
    print(str(1.56E9 * pow(10, I * -0.5)))
    ax1.plot(data[X_tag], data[Y1_tags])
    tags.append("E = {0:.3e}".format(1.56E9 * pow(10, I * -1)))
    I += 1

v = [0] * 7 * 2

plt.legend(tags)

ax1_d = [min(data[Y1_tags]), max(data[Y1_tags])]
ax2_d = [min(data[Y2_tags]), max(data[Y2_tags])]

plt.xlim([min(data[X_tag]), max(data[X_tag])])

# Plot Line2 (Right Y Axis)

ax1.set_title("Харрактер колебательного процесса в зависимости от\nмодуля обемной упругости рабочего тела", fontsize=22)

# Decorations
# ax1 (left Y axis)
ax1.set_xlabel(X_tag, fontsize=20)
ax1.tick_params(axis='x', rotation=0, labelsize=12)
ax1.set_ylabel(Y1_tags, color='tab:blue', fontsize=20)
ax1.tick_params(axis='y', rotation=0, labelcolor='tab:blue' )
ax1.grid(alpha=.2)


#ax1.set_ylim([ax1_d[0] * 1.01 - ax1_d[1] * 0.01, ax1_d[1] * 1.01 - ax1_d[0] * 0.01])

if False :

    ax2 = ax1.twinx()
    for data in datas :
        ax2.plot(data[X_tag], data[Y2_tags])

    ax2.set_ylabel(Y2_tags, color='tab:red', fontsize=20)

    ax2.tick_params(axis='y', labelcolor='tab:red')

    ax2.set_ylim([min(data[Y2_tags]), max(data[Y2_tags])])

    ax2.set_ylim([ax2_d[0] * 1.01 - ax2_d[1] * 0.01, ax2_d[1] * 1.01 - ax2_d[0] * 0.01])







#fig.tight_layout()

plt.tight_layout()

plt.show()



