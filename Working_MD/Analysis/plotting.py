#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

# USAGE:

# PREAMBLE:

import numpy as np
import matplotlib.pyplot as plt
import sys

# SUBROUTINES

# MAIN PROGRAM:

time = np.loadtxt(time_file)
datalist = np.loadtxt(data_file)

plt.plot(time[:], datalist[:])

plt.title('TITLE')
plt.xlabel('XLABEL')
plt.ylabel('YLABEL')

plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')

plt.savefig('figure.eps')

