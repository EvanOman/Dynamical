from numpy import *
import scipy.ndimage
import matplotlib.image
from PIL import Image
import sage.all
from sage.plot.contour_plot import ContourPlot
from sage.doctest.util import Timer
import time
from pylab import *

def mandel(n, m, itermax, xmin, xmax, ymin, ymax):

	start = time.time()

	'''
	Fast mandelbrot computation using numpy.
	
	(n, m) are the output image dimensions
	itermax is the maximum number of iterations to do
	xmin, xmax, ymin, ymax specify the region of the
	set to compute.
	'''
	# The point of ix and iy is that they are 2D arrays
	# giving the x-coord and y-coord at each point in
	# the array. The reason for doing this will become
	# clear below...
	ix, iy = mgrid[0:n, 0:m]
	# Now x and y are the x-values and y-values at each
	# point in the array, linspace(start, end, n)
	# is an array of n linearly spaced points between
	# start and end, and we then index this array using
	# numpy fancy indexing. If A is an array and I is
	# an array of indices, then A[I] has the same shape
	# as I and at each place i in I has the value A[i].
	x = linspace(xmin, xmax, n)[ix]
	y = linspace(ymin, ymax, m)[iy]
	# c is the complex number with the given x, y coords
	c = x+complex(0,1)*y
	del x, y # save a bit of memory, we only need z
	# the output image coloured according to the number
	# of iterations it takes to get to the boundary
	# abs(z)>2
	img = zeros(c.shape, dtype=int)
	# Here is where the improvement over the standard
	# algorithm for drawing fractals in numpy comes in.
	# We flatten all the arrays ix, iy and c. This
	# flattening doesn't use any more memory because
	# we are just changing the shape of the array, the
	# data in memory stays the same. It also affects
	# each array in the same way, so that index i in
	# array c has x, y coords ix[i], iy[i]. The way the
	# algorithm works is that whenever abs(z)>2 we
	# remove the corresponding index from each of the
	# arrays ix, iy and c. Since we do the same thing
	# to each array, the correspondence between c and
	# the x, y coords stored in ix and iy is kept.
	ix.shape = n*m
	iy.shape = n*m
	c.shape = n*m
	# we iterate z->z^2+c with z starting at 0, but the
	# first iteration makes z=c so we just start there.
	# We need to copy c because otherwise the operation
	# z->z^2 will send c->c^2.
	z = copy(c)
	
	#print "z: ", z + 2
	#print "c: ", c 

	for i in xrange(itermax):
		if not len(z): break # all points have escaped
		# equivalent to z = z*z+c but quicker and uses
		# less memory
		multiply(z, z, z)
		add(z, c, z)
		# these are the points that have escaped
		rem = abs(z)>2.0
		# colour them with the iteration number, we
		# add one so that points which haven't
		# escaped have 0 as their iteration number,
		# this is why we keep the arrays ix and iy
		# because we need to know which point in img
		# to colour
		img[ix[rem], iy[rem]] = i+1
		# -rem is the array of points which haven't
		# escaped, in numpy -A for a boolean array A
		# is the NOT operation.
		rem = -rem
		# So we select out the points in
		# z, ix, iy and c which are still to be
		# iterated on in the next step
		z = z[rem]
		ix, iy = ix[rem], iy[rem]
		c = c[rem]
	print 'Time taken:', time.time()-start
	
	img[img==0] = 101

	print img.T

	image = imshow(img.T, origin='lower left')
	image.write_png('mandel.png', noscale=True)

#Creates a 2 dimensional PARAMETER image of a singular perturbation of the complex quadratic map
##z = z**n1 + c2 + (beta)/(z.conjugate()**d)
def singPertParam(n1,d,beta, angle, n, m, itermax, xmin, xmax, ymin, ymax,filename,colorMap):
	start = time.time()

	#First we need to calculate the critical point that we are going to iterate on:
	r = ((d/n1)*abs(beta))**(1/(n1+d))

	#The set of critical points is a circle with the radius given by above, for now we will take the right-most real value on this circle:
	critPoint = r*complex(sage.all.cos(angle), sage.all.sin(angle))
	ix, iy = mgrid[0:n, 0:m]
	
	x = linspace(xmin, xmax, n)[ix]
	y = linspace(ymin, ymax, m)[iy]
	
	c = x+complex(0,1)*y
	del x, y # save a bit of memory, we only need z
	
	img = zeros(c.shape, dtype=int)
	
	ix.shape = n*m
	iy.shape = n*m
	c.shape = n*m
	
	#Sets all of the z_0 values
	if critPoint == 0:
		z = critPoint**n1 + copy(c)
	else:
		z = critPoint**n1 + copy(c) + (beta)/(critPoint.conjugate()**d)


	for i in xrange(itermax):
		if not len(z): break # all points have escaped
		
		#multiply(z, z, z)
		#add(z, c, z)

		z = z**n1 + c + (beta)/(z.conjugate()**d)

		# these are the points that have escaped
		rem = abs(z)>2.0
		
		img[ix[rem], iy[rem]] = i+1
		
		rem = -rem
		
		z = z[rem]
		ix, iy = ix[rem], iy[rem]
		c = c[rem]
	print 'Time taken:', time.time()-start

	img[img==0] = 101
	
	image = imshow(img.T, origin='lower left')
	if not colorMap == "":
		image.set_cmap(colorMap)
		image.write_png(filename+'.png', noscale=True)


#Creates a 2 dimensional PHASE image of a singular perturbation of the complex quadratic map
##z = z**n1 + c2 + (beta)/(z.conjugate()**d)
def singPertPhase(n1,d,beta,c2,n, m, itermax, xmin, xmax, ymin, ymax,filename,colorMap):
	start = time.time()

	ix, iy = mgrid[0:n, 0:m]

	x = linspace(xmin, xmax, n)[ix]
	y = linspace(ymin, ymax, m)[iy]

	z = x+complex(0,1)*y
	del x, y # save a bit of memory, we only need z

	img = zeros(z.shape, dtype=int)

	ix.shape = n*m
	iy.shape = n*m
	z.shape = n*m

	#Sets all of the z_0 values
	#z = z**n1 + c2 + (beta)/(z.conjugate()**d)

	for i in xrange(itermax):
		if not len(z): break # all points have escaped

		#multiply(z, z, z)
		#add(z, c, z)

		z = z**n1 + c2 + (beta)/(z.conjugate()**d)

		# these are the points that have escaped
		rem = abs(z)>2.0

		img[ix[rem], iy[rem]] = i+1

		rem = -rem

		z = z[rem]
		ix, iy = ix[rem], iy[rem]
		#c = c[rem]
	print 'Time taken:', time.time()-start

	img[img==0] = 101

	image = imshow(img.T, origin='lower left')
	image.set_cmap(colorMap)
	return image
	image.write_png(filename+'.png', noscale=True)


def burningShip(n, m, itermax, xmin, xmax, ymin, ymax,filename):

	start = time.time()

	'''
	Fast mandelbrot computation using numpy.
	
	(n, m) are the output image dimensions
	itermax is the maximum number of iterations to do
	xmin, xmax, ymin, ymax specify the region of the
	set to compute.
	'''
	
	ix, iy = mgrid[0:n, 0:m]
	
	x = linspace(xmin, xmax, n)[ix]
	y = linspace(ymin, ymax, m)[iy]
	
	c = x+complex(0,1)*y
	del x, y # save a bit of memory, we only need z
	
	img = zeros(c.shape, dtype=int)
	
	ix.shape = n*m
	iy.shape = n*m
	c.shape = n*m
	
	z = copy(c)

	for i in xrange(itermax):
		if not len(z): break # all points have escaped
		
		z = (abs(real(z)) + abs(imag(z))*complex(0,1))**2 + c

		# these are the points that have escaped
		rem = abs(z)>2.0
		
		img[ix[rem], iy[rem]] = i+1
	
		rem = -rem
		
		z = z[rem]
		ix, iy = ix[rem], iy[rem]
		c = c[rem]
	print 'Time taken:', time.time()-start
	img[img==0] = 101
	image = imshow(img.T, origin='lower left')
	image.write_png(filename+'.png', noscale=True)

"""
	Creates an image of the specified type, size, and parameters

	@param imgType-     Type of the image to be created. Two options for now:
	                         1. "param": creates a parameter space image(for now just on the c value)
	                         2. "phase": creates a phase space image(for now just on the c value)
	@param xSize-       number of x pixels
	@param ySize-       number of y pixels
	@param xMin/xMax-   min/max x values to check
	@param yMin/yMax-   min/max y values to check
	@param itermax-     max number of iterations to try
	@param n/d/beta/c-  the values of the paramer to be used in the function: z = z**n + c + (beta)/(z.conjugate()**d)
	@param angle-       the angle of the value along the crtical unit circle to use
	@param colorMap     the color map with which to render the image (optional)

"""
def makeImage(imgType, xSize, ySize, xMin, xMax, yMin, yMax, itermax, n, d, beta, c, angle, colorMap):

	#First we will make a consistent file name for the output image:
	filename = imgType + "--n=" + str(n) + "--d=" + str(d) + "--beta=" + str(beta) + "--angle=" + str(angle)

	#Now we call the appropriate function based on the user input
	if imgType == "param":
		singPertParam(n, d, beta, angle, xSize, ySize, itermax, xMin, xMax, yMin, yMax, filename, "spectral")
	elif imgType == "phase":
		singPertPhase(n, d, beta, c, xSize, ySize, itermax, xMin, xMax, yMin, yMax, filename, "spectral")
	else:
		print("Not a valid image type, doing nothing")


#Would be nice to have a single function to call in which I can specify what image I want to create(param or phase), which parameter to vary over(might cover the first item, changing varying over c rather than z is the difference between phase and parameter pictures), and maybe even which function to use to iterate. 

#Accomplish ^^ partially with some kind of GUI?

#Best bet for ^^ would be to use web, would be really slow unless it ran python in the background(possible?)

#Need a good image naming convention which automatically stores the parameters used and the type of picture made as the file name

#Use image metadata for ^^?

#Would be even better to not only make a filename with that info but also to push the new image to a latex document which could contain an acutal text description of the image as the figure caption.

#Parameter space is like the original Mandelbrot set: we choose all the parameters but one, let that parameter vary, and then iterate a critical point to check if it stays bounded. For the z^2 + c case we only have one critical point: z=0. In the singular perturbations case we have a set of critical points which are the points along the circle given by: r = ((d/n1)*abs(beta))**(1/(n1+d))

#Phase space is like the Julia sets: we choose all the parameters and then iterate different z values.