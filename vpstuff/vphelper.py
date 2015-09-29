from math import log
import os, re, traceback, sys
from vpstuff.termcon import TerminalController, ProgressWrapper
from time import gmtime, strftime
term = TerminalController()

ERROR_LOG_FN = 'error.fitcmp.log'

def termBold(msg):
	return term.render('${BOLD}%s${NORMAL}') %msg
	
def termWarn(msg):
	return term.render('${YELLOW}%s${NORMAL}') %msg

def show_error(msg, exit = False):
	print term.render('${RED}${BOLD}Error: ${NORMAL} %s' % msg)
	exc_type, exc_value, exc_traceback = sys.exc_info()
	infomsg = ''
	if exc_type:
		infomsg = 'Additional debugging information is available - written to %s' % ERROR_LOG_FN
		lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
		debugmsg = ''.join('   ' + line for line in lines) 
	if exit:
		infomsg = 'Exiting. ' + infomsg
	print infomsg
	
	error_log = open(ERROR_LOG_FN, 'a+')
	error_log.write('[%s] Error: %s\n\n' % (strftime("%Y-%m-%d %H:%M:%S", gmtime()), msg))
	if exc_type:
		error_log.write(debugmsg+'\n')
	if exit:
		error_log.write('Exiting.\n\n')
		sys.exit()

def isNotString(thestring, stringlist):
	for thisstring in stringlist:
		if (thisstring == thestring):
			return False
	return True

def readAtom(filename):
	ionlist = []
	ignoremarkers = ['#Q', 'end', '!', '??', '<<', '>>', '<>', '__']
	atomfile = open(filename, 'r')
	lines = atomfile.readlines()
	for iondat in lines:
		iondat = iondat.split()
		# ignore start, comment and end markers
		if (isNotString(iondat[0], ignoremarkers)):
			# check for space in ion name
			try:
				ion = [iondat[0], float(iondat[1]), float(iondat[2])]
			except ValueError:
				ion = [iondat[0] + ' ' + iondat[1], float(iondat[2]), float(iondat[3])]
			ionlist.append(ion)
	return ionlist

def nameList(ionlist):
	"""http://mail.python.org/pipermail/python-list/2000-February/025219.html"""
	uniquedict = {}
	for i in ionlist:
		uniquedict[i[0]] = 0
	nlist = uniquedict.keys()
	nlist.sort()
	return nlist
	
def filterIons(ionlist, includenames):
	filterlist = []
	for ion in ionlist:
		if (ion[0] in includenames): filterlist.append(ion)
	return filterlist

def labs2lobs(labs, z):
	return labs*(z+1.0)

def lobs2labs(lobs, z):
	return lobs/(z+1.0)

def zobsList(ionlist, zlist):
	"""Returns [[zvalue1, ionlist1], ...]"""
	lobslist = []
	for z in zlist:
		tmplist = []
		for ion in ionlist:
			tmplist.append([ion[0], labs2lobs(ion[1], z), ion[2]])
		lobslist.append([z, tmplist])
	return lobslist

def truncWL(obslist, minw, maxw):
	tlist = []
	for zi in obslist:
		tmplist = []
		for ion in zi[1]:
			if (ion[1] > minw and ion[1] < maxw):
				tmplist.append(ion)
		if (tmplist != []): tlist.append([zi[0], tmplist])
	return tlist

def minfF(obslist, minF):
	mlist = []
	for zi in obslist:
		tmplist = []
		for ion in zi[1]:
			if (ion[2] > minF):
				tmplist.append(ion)
		if (tmplist != []): mlist.append([zi[0], tmplist])
	return mlist

def printIL(name, wl, f, z, flag):
	print '%s  %s (%4.2f)    %4.3f    f = %.5e' % (flag, name.ljust(9), lobs2labs(wl, z), wl, f)

printIon = lambda ion, z, flag = ' ': printIL(ion[0], ion[1], ion[2], z, flag)

def printSuspects(atomfilename, namelist, zlist, minw, maxw, wl, tol = 0.2):
	obslist = truncWL(zobsList(filterIons(readAtom(atomfilename), namelist), zlist), minw, maxw)
	suspects = []
	susnames = []
	susz = []
	suszn = []
	for zi in obslist:
		for ion in zi[1]:
			if (ion[1] > wl-tol and ion[1] < wl+tol):
				suspects.append(ion)
				suszn.append([zi[0], ion[0]])
				if (ion[0] not in susnames): susnames.append(ion[0])
				if (zi[0] not in susz): susz.append(zi[0])
	# now print the results
	# print suspects
	# print suszn
	# print susz
	for zi in obslist:
		if (zi[0] in susz):
			print '\nz = %.5f\n-----------' % zi[0]
			for ion in zi[1]:
				if ([zi[0], ion[0]] in suszn):
					if (ion in suspects):
						printIon(ion, zi[0], '*')
					else:
						printIon(ion, zi[0])

sign = lambda x: cmp(x,0)
""" From: http://mail.python.org/pipermail/python-list/2008-August/504465.html """

def Nsum(nlist):
	thesum = 0.0
	for nitem in nlist:
		thesum += sign(nitem)*10**abs(nitem)
	return log(thesum, 10)

def TnTn(n1, n2, n3, n4):
	an1 = Nsum([n1, -n2])
	an3 = Nsum([n3, -n4])
	an2 = n2
	an4 = n4
	print Nsum([an1, an2, an3, an4])
	print an2
	print an3
	print an4

def open_atom():
	if not os.getenv('ATOMDIR'):
		raise Exception("Could not open file $ATOMDIR")
	afn = os.getenv('ATOMDIR')
	if not os.path.isfile(afn):
		raise Exception("Could not open file $ATOMDIR = %s" % afn)
	try:
		atom = open(afn).readlines()
	except:
		raise Exception("Could not open file $ATOMDIR = %s" % afn)
	return atom

def parse_atom():
	atom = open_atom()
	
	s0p, s1p = "[\s^\n]*", "[\s^\n]+" # white space
	ion = "([\.\w\* <>\?]+?)"         # species
	elt = "([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)" # E float
	# the format is ion lambda f gamma mass q comments
	
	#                 ion       wv        f         gamma     mass (optional)    q (optional)
	reptrn = "^"+s0p  +ion+s1p  +elt+s1p  +elt+s1p  +elt+s1p  +'('+elt+s1p+')?'  +'('+elt+s1p+')?'
	
	atomlst = []
	for a in atom:
		p = re.findall(reptrn, a)
		if p:
			p0 = p[0]
			pai = {'ion': p0[0], 'wv': p0[1], 'f': p0[3], 'gamma': p0[5], 'mass': p0[8], 'q': p0[11], 'linecopy': a}
			# print pai
			atomlst.append(pai)
			
	return atomlst

def find_line(species, approx_rest_wave):
	atom_sub_sp = find_lines_byspecies(species)
	deltawv = 9999999.0
	dwv_i = None
	for i, a in enumerate(atom_sub_sp):
		cdwv = abs(float(a['wv']) - float(approx_rest_wave))
		if deltawv > cdwv:
			deltawv = cdwv
			dwv_i = i
	
	return atom_sub_sp[dwv_i]

def find_lines_byspecies(species):
	atom = parse_atom()
	atom_sub_sp = [a for a in atom if a['ion'] == species]
	if not atom_sub_sp:
		raise Exception("Could not find line matching '%s %s' (no ion match) in atom.dat file" % (species, approx_rest_wave))
	
	return atom_sub_sp

def RepresentsInt(s):
	try: 
		int(s)
		return True
	except ValueError:
		return False