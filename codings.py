import sys

def singpert(x, n, d, c, beta):
	return x**n + c + beta/(x**d)

def makeEncodings(n ,d, beta, minX, maxX, iterNum):
	critR = .17728
	scaleIter = (float(maxX) - float(minX))/iterNum

	for i in xrange(iterNum):
		c = maxX - scaleIter * i
		codingStr = ""
		myMap = lambda x: singpert(x, n, d, c, beta)
		xVal = critR
		for i in xrange(50):
			if xVal > critR:
				codingStr += " R "
				if xVal > 100:
					continue
			elif xVal == critR:
				codingStr += " C "
			elif xVal > 0 and xVal < critR:
				codingStr += " L "
			elif xVal < 0 and xVal > -critR:
				codingStr += " r "
			elif xVal < 0 and xVal < -critR:
				codingStr += " l "
			xVal = myMap(xVal)
		print(codingStr + ":\t\t" + str(c))



if __name__ == "__main__":
	makeEncodings(2, 2, .001, float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]))
