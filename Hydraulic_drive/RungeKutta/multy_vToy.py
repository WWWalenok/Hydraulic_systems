
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import math
import sys



print(sys.argv)

files = sys.argv[2:]

if len(files) == 0:
    exit(10)

warnings.filterwarnings(action='once')
data = 0
for name in files:
    tdata = np.genfromtxt(
        "out_" + name + ".csv",
        delimiter=',',
        skip_header=0,
        names=True
    )
    data = tdata

large = 20; med = 14; small = 12
params = {'axes.titlesize': large,
          'legend.fontsize': med,
          'figure.figsize': (10, 4),
          'axes.labelsize': med,
          'axes.titlesize': med,
          'xtick.labelsize': med,
          'ytick.labelsize': med,
          'figure.titlesize': large}


plt.rcParams.update(params)
plt.style.use('seaborn-whitegrid')
sns.set_style("white")

# Import Data

X_tag = 'U'
Y1_tags = 'V'
Y2_tags = 'U'

# Plot Line1 (Left Y Axis)
fig, ax1 = plt.subplots(1, 1, figsize=(10,5))


I = 0

tags = []

ax1.plot(data[X_tag], data[Y1_tags])
tags.append(files[I])

v = [0] * 7 * 2

plt.legend(tags)


# Plot Line2 (Right Y Axis)

ax1.set_title("", fontsize=0)

# Decorations
# ax1 (left Y axis)
ax1.set_xlabel('Управление.', fontsize=20)
ax1.tick_params(axis='x', rotation=0, labelsize=12)
ax1.set_ylabel('Скорость: V  m/s.', fontsize=20)
ax1.tick_params(axis='y', rotation=0)
ax1.grid(alpha=.2)


#ax1.set_ylim([ax1_d[0] * 1.01 - ax1_d[1] * 0.01, ax1_d[1] * 1.01 - ax1_d[0] * 0.01])

if 1 :

    ax2 = ax1.twinx()
    ax2.plot(data[X_tag], data[Y2_tags], color = 'tab:red')

    ax2.set_ylabel("Управление", color='tab:red', fontsize=20)

    ax2.tick_params(axis='y', labelcolor='tab:red')







#fig.tight_layout()

plt.tight_layout()

plt.savefig("E:\\Repos\\Latex\\rsgc_1\\" + sys.argv[1] + ".png", format = "png")

plt.show()



