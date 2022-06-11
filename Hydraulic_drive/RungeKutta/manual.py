
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import math

warnings.filterwarnings(action='once')

data = np.genfromtxt(
    "out_real_stab.csv",
    delimiter=',',
    skip_header=0,
    names=True
)

out_name = 'phase_real'

large = 20; med = 14; small = 12
params = {'axes.titlesize': large,
          'legend.fontsize': med,
          'figure.figsize': (12, 4),
          'axes.labelsize': med,
          'axes.titlesize': med,
          'xtick.labelsize': med,
          'ytick.labelsize': med,
          'figure.titlesize': large}


plt.rcParams.update(params)
plt.style.use('seaborn-whitegrid')
sns.set_style("white")

# Import Data

X_tag = 'V'
Ys_tags = [['A'], ['P1'], ['P2'], [], []]
Zs_tags = [[], [], [], [], []]

tag_name = {'A': 'Ускорение', 'P1': 'Давление в поршневой', 'P2': 'Давление в штоковой', 'V' : 'Скорость'}

print(tag_name['A'])

# Plot Line1 (Left Y Axis)

axs = []

axs.append(plt.subplot(131))
axs.append(plt.subplot(132))
axs.append(plt.subplot(133))

for i in range(len(axs)):
    ax = axs[i]
    Y_tags = Ys_tags[i]
    Z_tags = Zs_tags[i]

    lines = []
    tags = []

    for Y_tag in Y_tags:
        line, = ax.plot(data[X_tag][3:], data[Y_tag][3:], label=Y_tag)
        lines.append(line)
        if Y_tag in tag_name:
            tags.append(tag_name[Y_tag])
        else:
            tags.append(Y_tag)

    ax.set_xlabel(tag_name[X_tag], fontsize=16)
    ax.tick_params(axis='x', rotation=0, labelsize=12)
    ax.set_ylabel(tag_name[Y_tags[0]], fontsize=16)
    ax.tick_params(axis='y', rotation=0, labelcolor='tab:blue' )
    ax.grid(alpha=.2)



    if len(Z_tags) > 0:
        axZ = ax.twinx()
        axZ.set_ylabel(Z_tags, fontsize=16)
        for Z_tag in Z_tags:
            line, = axZ.plot(data[X_tag][3:], data[Z_tag][3:], color = 'C{a}'.format(a = len(Y_tags)), label=Z_tag)
            lines.append(line)
            tags.append(Z_tag)

    #ax.legend(handles = lines, labels = tags)

plt.tight_layout()

plt.savefig("E:\\Repos\\Latex\\rsgc_1\\" + out_name + ".png", format = "png")

plt.show()



