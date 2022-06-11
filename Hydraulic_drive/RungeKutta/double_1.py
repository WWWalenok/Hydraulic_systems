
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
out_name = ''
if(len(sys.argv) > 1):
    out_name = sys.argv[1]
else:
    out_name = 'base'
if len(files) == 0:

    files = ['real', 'linear']

warnings.filterwarnings(action='once')
datas = []
for name in files:
    data = np.genfromtxt(
        "out_" + name + ".csv",
        delimiter=',',
        skip_header=0,
        names=True
)
    datas.append(data)

large = 20; med = 14; small = 12
params = {'axes.titlesize': large,
          'legend.fontsize': med,
          'figure.figsize': (10, 6),
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
Y1_tags = 'V'
Y2_tags = 'P1'
Y3_tags = 'P2'

# Plot Line1 (Left Y Axis)
ax1 = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)
I = 0

tags = []

for data in datas:
    ax1.plot(data[X_tag], data[Y1_tags])
    tags.append(files[I])
    I += 1

v = [0] * 7 * 2

plt.legend(tags, loc = 'lower right')


# Plot Line2 (Right Y Axis)

ax1.set_title("", fontsize=0)

# Decorations
# ax1 (left Y axis)
ax1.set_xlabel('Время: T s.', fontsize=20)
ax1.tick_params(axis='x', rotation=0, labelsize=12)
ax1.set_ylabel(Y1_tags, fontsize=20)
ax1.tick_params(axis='y', rotation=0)

ax1.legend(tags)

for data in datas :
    ax2.plot(data[X_tag], data[Y2_tags])

ax2.set_ylabel(Y2_tags, fontsize=20)

ax2.tick_params(axis='y')

ax2.legend(tags)


for data in datas :
    ax3.plot(data[X_tag], data[Y3_tags])

ax3.set_ylabel(Y3_tags, fontsize=20)

ax3.tick_params(axis='y')

ax3.legend(tags)


ax1.grid(alpha=.2)
ax2.grid(alpha=.2)
ax3.grid(alpha=.2)

#fig.tight_layout()

plt.tight_layout()

plt.savefig("E:\\Repos\\Latex\\rsgc_1\\" + out_name + ".png", format = "png")

plt.show()



