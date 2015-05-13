#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

# USAGE:

# PREAMBLE:

import numpy as np
import matplotlib.pyplot as plt
import sys

time_file = sys.argv[1]
temp_file = sys.argv[2]

# SUBROUTINES

# MAIN PROGRAM:

time = np.loadtxt(time_file)
datalist = np.loadtxt(temp_file)

plt.plot(time[:], datalist[:])

plt.title('Temperature of Liquid Argon System')
plt.xlabel(r'Time $(ps)$')
plt.ylabel(r'Temperature $(K)$')

plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')

plt.savefig('tempvstime.eps')

