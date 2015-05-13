#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python

# USAGE:

# PREAMBLE:

import numpy as np
import sys

log_file = sys.argv[1]
system = sys.argv[2]

# SUBROUTINES

# MAIN PROGRAM:

lf = open('%s' %(log_file), 'r')
tf = open('%s.time.dat' %(system), 'w')
ef = open('%s.tot_en.dat' %(system), 'w')
pf = open('%s.pot_en.dat' %(system), 'w')
kf = open('%s.kin_en.dat' %(system), 'w')
tempf = open('%s.temp.dat' %(system), 'w')

for line in lf:  	# for statement that flows through all lines in the .log file

	if line.startswith('Step'):   	# collect value for time conversion unit (1 step per n fs) 
		t = line.split() 	# creates a list that contains all components of the line
		step = float(t[-1])
		tf.write('%15.6f \n' %(step*0.002)) 

	if line.startswith('Total Energy'): 	 
		a = line.split() 	
		en = float(a[-1])
		ef.write('%15.6f \n' %(en))
	
	if line.startswith('Total Potential'):
		b = line.split()
		pot = float(b[-1])
		pf.write('%15.6f \n' %(pot))
	
	if line.startswith('Total Kinetic'):
		c = line.split()
		kin = float(c[-1])
		kf.write('%15.6f \n' %(kin))

	if line.startswith('Instantaneous'):
		d = line.split()
		temp = float(d[-1])
		tempf.write('%15.6f \n' %(temp))

lf.close()
tf.close()
ef.close()
pf.close()
kf.close()
tempf.close()

