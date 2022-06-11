
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import math

warnings.filterwarnings(action='once')
out_name = 'base'

names = ['real', 'linear 1']
datas = []
for name in names:
    datas.append(np.genfromtxt(
        "out_" + name + ".csv",
        delimiter=',',
        skip_header=0,
        names=True
    ))
names = ['diff model', 'linear']


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

X_tag = 'T'
Ys_tags = [['V'],['U'],[],[],[],[],[]]
Zs_tags = [[],[],[],[],[],[],[],[]]
axs = []

# Plot Line1 (Left Y Axis)

axs.append(plt.subplot(211))
axs.append(plt.subplot(212))
#axs.append(plt.subplot(313))
for i in range(len(axs)):
    tags = []
    lines = []
    for j in range(len(datas)):
        data = datas[j]
        name = names[j]
        ax = axs[i]
        Y_tags = Ys_tags[i]
        Z_tags = Zs_tags[i]

        for Y_tag in Y_tags:
            line, = ax.plot(data[X_tag], data[Y_tag], label=Y_tag, color = 'C{a}'.format(a = len(lines)))
            lines.append(line)
            tags.append(name + ": " + Y_tag)

        ax.set_xlabel(X_tag, fontsize=20)
        ax.tick_params(axis='x', rotation=0, labelsize=12)
        ax.set_ylabel(Y_tags, fontsize=20)
        ax.tick_params(axis='y', rotation=0)
        ax.grid(alpha=.2)

        if len(Z_tags) > 0:
            axZ = ax.twinx()
            axZ.set_ylabel(Z_tags, color='tab:blue', fontsize=20)
            for Z_tag in Z_tags:
                line, = axZ.plot(data[X_tag], data[Z_tag], color = 'C{a}'.format(a = len(lines)), label=Z_tag)
                lines.append(line)
                tags.append(name + ": " + Z_tag)

    ax.legend(handles = lines, labels = tags)

plt.tight_layout()

plt.savefig("E:\\Repos\\Latex\\rsgc_1\\" + out_name + ".png", format = "png")

plt.show()



