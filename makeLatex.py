from glob import *
import re

def format_number(num):
	tup = dec.as_tuple()
	delta = len(tup.digits) + tup.exponent
	digits = ''.join(str(d) for d in tup.digits)
	if delta <= 0:
		zeros = abs(tup.exponent) - len(tup.digits)
		val = '0.' + ('0'*zeros) + digits
	else:
		val = digits[:delta] + ('0'*tup.exponent) + '.' + digits[delta:]
	val = val.rstrip('0')
	if val[-1] == '.':
		val = val[:-1]
	if tup.sign:
		return '-' + val
	return val

def main():
	print("Here are the png's in this directory:")
	for filename in glob('./*.png'):
		#Gets rid of the forward slash and the file extension
		newstr = filename[2:len(filename) - 4]

		imgType = newstr[0:5]

		values = re.findall(r"[-+]?\d*\.\d+|\d+", filename)

		print("This is a " + imgType + " space image with values: $n=" + values[0] + "$, $d=" + values[1] + "$, and $\\beta=" + values[2] + "$")

if __name__ == "__main__":
    main()