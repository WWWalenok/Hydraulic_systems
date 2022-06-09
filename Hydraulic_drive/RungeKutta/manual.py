print("Hi")


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import math

warnings.filterwarnings(action='once')

data = np.genfromtxt(
    "out_real.csv",
    delimiter=',',
    skip_header=0,
    names=True
)

print(data)


large = 20; med = 14; small = 12
params = {'axes.titlesize': large,
          'legend.fontsize': med,
          'figure.figsize': (12, 10),
          'axes.labelsize': med,
          'axes.titlesize': med,
          'xtick.labelsize': med,
          'ytick.labelsize': med,
          'figure.titlesize': large}


plt.rcParams.update(params)
plt.style.use('seaborn-whitegrid')
sns.set_style("white")

# Import Data

X_tag = 'T'
Ys_tags = [['QI', 'QS'], ['Q1', 'Q2'], ['U', 'Y']]
Zs_tags = [[], [], ['V']]
axs = []

# Plot Line1 (Left Y Axis)

axs.append(plt.subplot(311))
axs.append(plt.subplot(312))
axs.append(plt.subplot(313))

for i in range(len(axs)):
    ax = axs[i]
    Y_tags = Ys_tags[i]
    Z_tags = Zs_tags[i]

    lines = []


    for Y_tag in Y_tags:
        line, = ax.plot(data[X_tag], data[Y_tag], label=Y_tag)
        lines.append(line)

    ax.set_xlabel(X_tag, fontsize=20)
    ax.tick_params(axis='x', rotation=0, labelsize=12)
    ax.set_ylabel(Y_tags, color='tab:blue', fontsize=20)
    ax.tick_params(axis='y', rotation=0, labelcolor='tab:blue' )
    ax.grid(alpha=.2)

    tags = Y_tags

    if len(Z_tags) > 0:
        axZ = ax.twinx()
        axZ.set_ylabel(Z_tags, color='tab:blue', fontsize=20)
        for Z_tag in Z_tags:
            line, = axZ.plot(data[X_tag], data[Z_tag], color = 'C{a}'.format(a = len(Y_tags)), label=Z_tag)
            lines.append(line)

    ax.legend(handles = lines, labels = Y_tags + Z_tags)

plt.tight_layout()

plt.show()



