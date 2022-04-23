print("Hi")


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.filterwarnings(action='once')
import math

datas = dict()

var_names = [
    't',
    'v',
    'x',
    'p1',
    'p2',
    'a'
]

for a in range(+5, 5 + 1):
    for b in range(-5, 5 + 1):
        for c in range(-5, 5 + 1):
            print("smple_{a}_{b}_{c}.csv".format(a=a, b=b, c=c))
            data = np.genfromtxt(
                "smple_{a}_{b}_{c}.csv".format(a = a, b = b, c = c),
                delimiter=',',
                skip_header=1,
                skip_footer=10,
                names=var_names
            )

            var_list = dict()
            for name in var_names:
                var_list[name] = data[name]
            datas["{a}_{b}_{c}".format(a=a, b=b, c=c)] = var_list


#print(datas)

large = 22; med = 16; small = 12
params = {'axes.titlesize': large,
          'legend.fontsize': med,
          'figure.figsize': (16, 10),
          'axes.labelsize': med,
          'axes.titlesize': med,
          'xtick.labelsize': med,
          'ytick.labelsize': med,
          'figure.titlesize': large}
plt.rcParams.update(params)
plt.style.use('seaborn-whitegrid')
sns.set_style("white")

X_tag = 't'
Y1_tag = 'v'
Y2_tag = 't'

# Plot Line1 (Left Y Axis)
fig, ax1 = plt.subplots(1, 1, figsize=(16,9))

for data in datas:
    #print(datas[data])
    ax1.plot(datas[data][X_tag], datas[data][Y1_tag], 'b')




# Decorations
# ax1 (left Y axis)
ax1.set_xlabel(X_tag, fontsize=20)
ax1.tick_params(axis='x', rotation=0, labelsize=12)
ax1.set_ylabel(Y1_tag, color='tab:blue', fontsize=20)
ax1.tick_params(axis='y', rotation=0, labelcolor='tab:blue' )
ax1.grid(alpha=.2)



plt.tight_layout()

plt.show()