# LogisticBifn.py:
#	Plot a bifurcation diagram for the Logistic map.

# Import modules
from numpy import *
from pylab import *

# Define the Logistic map's function 
def LogisticMap(r,x):
	return r * x * (1.0 - x)

# Setup parameter range
rlow  = 2.85
rhigh = 4.0

# Setup the plot
# It's numbered 1 and is 8 x 6 inches.
figure(1,(8,6))
# Stuff parameter range into a string via the string formating commands.
TitleString  = 'Logistic map, f(x) = r x (1 - x), '
TitleString += 'bifurcation diagram for r in [%g,%g]' % (rlow,rhigh)
title(TitleString)
# Label axes
xlabel('Control parameter r')
ylabel('{X(n)}')

# This command turns off matplotlib's interactive mode.
# If you want to see the graphics commands as they're created,
#	comment out this line and the ion() command at the end.
#ioff()

# We tell plot()'s autoscaling the range in (r,x) we want by
#	putting dots at the corners of the desired data window.
#	Otherwise, the plot constantly rescales with new data.
plot([rhigh], [1.0], 'k,')
plot([rhigh], [0.0], 'k,')
plot([rlow], [0.0], 'k,')
plot([rlow], [1.0], 'k,')
# Note that plot()'s autoscaling is erratic is how it sets the x scale,
#	sometimes adding extra space at low or high x.
# Should override the default autoscaling.

# Set the initial condition used across the different parameters
ic = 0.2
# Establish the arrays to hold the set of iterates at each parameter value
# The iterates we'll throw away
nTransients = 200
# This sets how much the attractor is filled in
nIterates = 250
# This sets how dense the bifurcation diagram will be
nSteps = 400
# Sweep the control parameter over the desired range
rInc = (rhigh-rlow)/float(nSteps)
for r in arange(rlow,rhigh,rInc):
	# Set the initial condition to the reference value
	state = ic
	# Throw away the transient iterations
	for i in xrange(nTransients):
		state = LogisticMap(r,state)
	# Now store the next batch of iterates
	rsweep = [ ]   # The parameter value
	x = [ ]        # The iterates
	for i in xrange(nIterates):
		state = LogisticMap(r,state)
		rsweep.append(r)
		x.append( state )
	plot(rsweep, x, 'k,') # Plot the list of (r,x) pairs as pixels

# Use this to save figure as a bitmap png file
#savefig('LogisticBifn', dpi=600)

# Turn on matplotlib's interactive mode.
#ion()

# Display plot in window
show()