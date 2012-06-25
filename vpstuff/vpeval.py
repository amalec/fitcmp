# based on http://lybniz2.sourceforge.net/safeeval.html

from math import *
from string import strip
import os, sys, re
from vpstuff.ncol import ColumnDensity as N

# make a list of safe functions
safe_list = ['math','acos', 'asin', 'atan', 'atan2', 'ceil', 
	'cos', 'cosh', 'de grees', 'e', 'exp', 'fabs', 'floor', 
	'fmod', 'frexp', 'hypot', 'ldexp', 'log', 'log10', 
	'modf', 'pi', 'pow', 'radians', 'sin', 'sinh', 'sqrt', 
	'tan', 'tanh']

# use the list to filter the local namespace
safe_dict = dict([ (k, locals().get(k, None)) for k in safe_list ])

# add any needed builtins, non-builtins
safe_dict['abs'] = abs
safe_dict['bin'] = bin
safe_dict['bool'] = bool
safe_dict['c'] = 299792458.0
safe_dict['chr'] = chr
safe_dict['cmp'] = cmp
safe_dict['complex'] = complex
safe_dict['enumerate'] = enumerate
safe_dict['False'] = False
safe_dict['filter'] = filter
safe_dict['float'] = float
safe_dict['format'] = format
safe_dict['hex'] = hex
safe_dict['int'] = int
safe_dict['len'] = len
safe_dict['list'] = list
safe_dict['map'] = map
safe_dict['max'] = max
safe_dict['min'] = min
safe_dict['N'] = N
safe_dict['oct'] = oct
safe_dict['range'] = range
safe_dict['reduce'] = reduce
safe_dict['reversed'] = reversed
safe_dict['round'] = round
safe_dict['slice'] = slice
safe_dict['sorted'] = sorted
safe_dict['str'] = str
safe_dict['sum'] = sum
safe_dict['True'] = True
safe_dict['tuple'] = tuple
safe_dict['zip'] = zip

def vpeval(eval_str):
	if eval_str == None or strip(eval_str) == '':
		return
	try:
		eval_str = re.sub("{[\s^\n]*([\.\w\* <>\?]+?)[\s^\n]+([0-9]*\.?[0-9]+)[\s^\n]*}", lookup_line, eval_str)
		# eval_str = re.sub("([0-9]*\.?[0-9]+)[nN]", calculate_N, eval_str)
		result =  eval(eval_str, {"__builtins__":None}, safe_dict)
		if type(result) == type(list()):
			for r in result:
				print r
		elif result == None or result == '':
			return
		else:
			print result
	except Exception as err:
		print err

def lookup_line(line_match):
	if not os.getenv('ATOMDIR'):
		raise Exception("Could not open file $ATOMDIR")
	afn = os.getenv('ATOMDIR')
	if not os.path.isfile(afn):
		raise Exception("Could not open file $ATOMDIR = %s" % afn)
	try:
		atom = open(afn).readlines()
	except:
		raise Exception("Could not open file $ATOMDIR = %s" % afn)
	
	species = line_match.group(1)
	rest_wave = line_match.group(2)
	
	for a in atom:
		sr = re.findall("^[\s^\n]*(%s)[\s^\n]+([0-9]*\.?[0-9]+)" % re.escape(species), a)
		if len(sr) != 0:
			float(sr[0][1])
			if sr[0][1].find(rest_wave) != -1:
				return sr[0][1]
	raise Exception("Could not find line matching '%s %s' in %s" % (species, rest_wave, afn))

# def calculate_N(N_match):
# 	return "N(%s)" % N_match.group(1)