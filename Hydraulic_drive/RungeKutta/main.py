print("Hi")


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings; warnings.filterwarnings(action='once')
import math



data = np.genfromtxt("out.txt", delimiter=',', skip_header=1,
                     skip_footer=10, names=['t', 'v', 'x', 'p1', 'p2', 'f', 'h', 'n', 'a', 'u', 'cu', 'dt'])

print(data)


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



# Import Data

X_tag = 't'
Y1_tags = 'p1'
Y2_tags = 'x'

# Plot Line1 (Left Y Axis)
fig, ax1 = plt.subplots(1, 1, figsize=(16,9))

ax1.plot(data[X_tag], data[Y1_tags])

v = [0] * 7 * 2

ax2 = ax1.twinx()
ax2.plot(data[X_tag], data[Y2_tags], color='red')

# Plot Line2 (Right Y Axis)

ax1.set_title("Зависимость", fontsize=22)

# Decorations
# ax1 (left Y axis)
ax1.set_xlabel(X_tag, fontsize=20)
ax1.tick_params(axis='x', rotation=0, labelsize=12)
ax1.set_ylabel(Y1_tags, color='tab:blue', fontsize=20)
ax1.tick_params(axis='y', rotation=0, labelcolor='tab:blue' )
ax1.grid(alpha=.2)

ax2.set_ylabel(Y2_tags, color='tab:red', fontsize=20)

ax2.tick_params(axis='y', labelcolor='tab:red')


#fig.tight_layout()

plt.tight_layout()

plt.show()



