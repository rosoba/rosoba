'''
Created on Jan 9, 2014

@author: rch
'''

import numpy as np
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots()
t = np.arange(0.01, 10.0, 0.01)
s1 = np.exp(t)
max_s1 = np.max(s1)
ax1.plot(t, s1, 'b-')
ax1.set_xlabel('time (s)')
ax1.set_ylim(0, max_s1 * 1.1)
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('exp', color='b')


ax2 = ax1.twinx()
max_s2 = max_s1 * 10.0
ax2.set_ylim(0, max_s2 * 1.1)
ax2.set_ylabel('sin', color='black')
plt.show()
