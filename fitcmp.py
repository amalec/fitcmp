# fitcmp by Adrian Malec
# Displays vpfit-like plot for two fort.13 files for comparison.
# WARNING: All the fitting regions in each file must correspond.
# WARNING: I started writing this code when I commenced my PhD, it is terrible. Good luck. (I've been slowly refactoring the code though...)

import sys, os, readline, re
from optparse import OptionParser

from matplotlib import use as useGUI
useGUI('TkAgg')

from vpstuff.vpcmp import *
from vpstuff.vpeval import vpeval
from vpstuff.termcon import TerminalController, ProgressWrapper
from string import ljust, strip, lower
from matplotlib.pyplot import ion, close


VERSION = 0.85

term = TerminalController()

def main():
	scriptpath = os.path.dirname(sys.argv[0])
	ion() # turn on matplotlib interactive mode
	
	usage = "Usage: %prog [options] <vpfit executable path>"
	ver = "%prog " + str(VERSION)
	parser = OptionParser(usage=usage, version=ver)
	parser.add_option("-o", "--oldf13", dest="oldf", help="Old fort.13 file (optional)", default=False)
	parser.add_option("-n", "--newf13", dest="newf", help="New fort.13 file (fort.13 by default)", default="fort.13")
    # parser.add_option("--vpfix1", dest="vpfix1", help="Fixes input issue with some versions of VPFIT", action="store_true", default=False)
	
	(options, args) = parser.parse_args()
	if args == []:
		print "No VPFIT executable path argument.\nNote: -v, --vpfit-path flag has been removed.\nA valid VPFIT executable path must now always be specified as an argument."
		parser.print_help()
		sys.exit()
	
    # vppath, oldf, newf, vpfix1 = args[0], options.oldf, options.newf, options.vpfix1
	vppath, oldf, newf, vpfix1 = args[0], options.oldf, options.newf, False
	
	def vpfix(n): # strange way of doing it, but may be handy later on, or alternatively I'll rewrite the fort.17 functions to automatically detect the vpfit version
		if n == 1:
			return vpfix1
	
	if not os.path.isfile(vppath):
		print "VPFIT executable: %s NOT FOUND" % vppath
		sys.exit()
	if not os.path.isfile(newf):
		print "New fort.13 file: %s NOT FOUND" % newf
		sys.exit()
	if not oldf:
		oldf = newf
	if not os.path.isfile(oldf):
		print "Old fort.13 file: %s NOT FOUND" % newf
		sys.exit()
	
	rftList_old, c_old, rftList_new, c_new = readAllf17(oldf, newf, vppath, vpfix)
	
	colour_config, tick_config = readColourConfig(scriptpath + '/colours.cfg')
	
	settings = {'plot_type': 1, 'flux_bottom': 1, 'tick_type': 1, 'group_delta': 0.07,
	            'tick_label': 1, 'plot_width': 650, 'plot_height': 490}
	
	# from: http://docs.python.org/library/readline.html
	histfile = os.path.join(os.path.expanduser("~"), ".f13cmphist")
	try:
	    readline.read_history_file(histfile)
	except IOError:
	    pass
	import atexit
	atexit.register(readline.write_history_file, histfile)
	#
	
	print "Available fitting regions:"
	
	listRegions(rftList_new)
	print
	selHelp()
	print
	fids = [[1, 1]] # Plot figure ID and fit ID
	
	showPlot(rftList_old[fids[-1][1]-1], rftList_new[fids[-1][1]-1], c_old, c_new, colour_config, tick_config, settings, figureid = fids[-1][0])
	
	while True:
		try:
			sin = raw_input("[" + termBold(str(fids[-1][1])) + "] ") # removed lower()
			sin = sin.split(';')
			for s in sin:
				s = strip(s)
				if s == '':
					showPlot(rftList_old[fids[-1][1]-1], rftList_new[fids[-1][1]-1], c_old, c_new, colour_config, tick_config, settings, figureid = fids[-1][0])
				elif s[0] == 'a':
					if len(s) == 1:
						autoold, autonew = True, True
					elif s == 'an':
						autoold, autonew = False, True
					elif s == 'ao':
						autoold, autonew = True, False
					else:
						print "Invalid 'a' command"
						continue
					print termWarn("[Entering autoreload mode, Ctrl+C, Enter to exit]")
					while 1:
						try:
							if autonew:
								rftList_new, c_new = reloadAllFits(newf, rftList_new, vppath, vpfix, silent = True)
							if autoold:
								rftList_old, c_old = reloadAllFits(oldf, rftList_old, vppath, vpfix, silent = True)
							showMultiPlots(fids, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings)
						except KeyboardInterrupt:
							print termWarn("[Exiting autoreload mode]")
							break
				elif s[0] == 'v':
					vcmd = re.findall("v(n|o)[\s^\n]*/[\s^\n]*(([a-zA-Z]{1,2})([\s^\n]*-[\s^\n]*([a-zA-Z0-9\s^\n,]+))?)([\s^\n]*/[\s^\n]*([+-]?[0-9]*\.?[0-9]+)[\s^\n]*,[\s^\n]*([+-]?[0-9]*\.?[0-9]+))?[\s^\n]*", s)
					if len(vcmd) == 0:
						print "Invalid 'v' command"
						continue
					else:
						# print vcmd[0]
						rftList_old, c_old, rftList_new, c_new = readAllf17(oldf, newf, vppath, vpfix)
						if vcmd[0][0] == 'n':
							rftList, comp = rftList_new, c_new
						else:
							rftList, comp = rftList_old, c_old
						tiedz_lbl = vcmd[0][2]
						if vcmd[0][6] != '':
							vlow = float(vcmd[0][6])
						else:
							vlow = None
						if vcmd[0][7] != '':
							vhigh = float(vcmd[0][7])
						else:
							vhigh = None
						if vcmd[0][4].strip() != '':
							splist = vcmd[0][4].strip().split(',')
							splist = [spl.strip() for spl in splist]
						else:
							splist = []
							# print splist
				
					showStackPlot(tiedz_lbl, rftList, comp, colour_config, tick_config, settings, splist, vlow = vlow, vhigh = vhigh, cursor_on = True)
				elif s[0:2] == 'e/':
					vpeval(s[2:])
				elif s == 'n':
					fids[-1][1] = fids[-1][1] % len(rftList_new) + 1
					showPlot(rftList_old[fids[-1][1]-1], rftList_new[fids[-1][1]-1], c_old, c_new, colour_config, tick_config, settings, figureid = fids[-1][0])
				elif s == 'p':
					fids[-1][1] = (fids[-1][1] - 2) % len(rftList_new) + 1
					showPlot(rftList_old[fids[-1][1]-1], rftList_new[fids[-1][1]-1], c_old, c_new, colour_config, tick_config, settings, figureid = fids[-1][0])
				elif s == 'q' or s == 'quit' or s == 'exit':
					print "Quitting fitcmp."
					sys.exit()
				elif s == '?' or s == 'h' or s == 'help':
					selHelp()
				elif s == 'l':
					listRegions(rftList_new)
				elif s == '!':
					selLines(oldf, newf, rftList_old[fids[-1][1]-1], rftList_new[fids[-1][1]-1], c_old, c_new)
				elif s[0] == '^':
					if len(s) == 1:
						fids.append([fids[-1][0]+1, fids[-1][1]])
					elif s == '^c':
						for fid in fids:
							close(fid[0])
						fids = [[1, fids[-1][1]]]
					else:
						print "Invalid figure window command"
						continue
					showPlot(rftList_old[fids[-1][1]-1], rftList_new[fids[-1][1]-1], c_old, c_new, colour_config, tick_config, settings, figureid = fids[-1][0])
				elif s[:3] == 'ps/':
					exportAll(s[3:], rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings, ext = 'ps')
				elif s[:4] == 'eps/':
					exportAll(s[4:], rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings, ext = 'eps')
				elif s[:4] == 'pdf/':
					exportAll(s[4:], rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings, ext = 'pdf')
				elif s == 'f13':
					printf13Paths(oldf, newf)
				elif s[0] == 'r':
					reloadmultifigs = True
					if len(s) == 1:
						rftList_old, c_old, rftList_new, c_new = readAllf17(oldf, newf, vppath, vpfix)
					elif s == 'rc':
						rftList_old, c_old, rftList_new, c_new = reloadCFits(fids[-1][1], oldf, rftList_old, newf, rftList_new, vppath, vpfix)
						reloadmultifigs = False
					elif s == 'rn':
						rftList_new, c_new = reloadAllFits(newf, rftList_new, vppath, vpfix)
					elif s == 'ro':
						rftList_old, c_old = reloadAllFits(oldf, rftList_old, vppath, vpfix)
					elif s == 'rcn':
						rftList_new, c_new = reloadFit(fids[-1][1], newf, rftList_new, vppath, vpfix)
						reloadmultifigs = False
					elif s == 'rco':
						rftList_old, c_old = reloadFit(fids[-1][1], oldf, rftList_old, vppath, vpfix)
						reloadmultifigs = False
					# else:
					# 	try:
					# 		rfid = int(s[1:])
					# 		if (rfid > len(rftList_new) or rfid < 1):
					# 			raise ValueError
					# 		rftList_old, c_old, rftList_new, c_new = reloadFit(rfid, oldf, rftList_old, c_old, newf, rftList_new, c_new, vppath, vpfix)
					# 	except ValueError:
					# 		print "Enter a valid fit id to reload"
					if reloadmultifigs:
						showMultiPlots(fids, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings)
					else:
						showPlot(rftList_old[fids[-1][1]-1], rftList_new[fids[-1][1]-1], c_old, c_new, colour_config, tick_config, settings, figureid = fids[-1][0])
				elif s[0] == '{':
					opt_match = re.match("{([\w ]+)?[, ]?([\w .]+)?}", s)
					if opt_match:
						if not opt_match.group(1) or not opt_match.group(2): 
							selOptionHelp()
							continue
					
						opt_name  = strip(opt_match.group(1))
						opt_value = strip(opt_match.group(2))
					
						if opt_name == 'plot_type' or opt_name == 'pt':
							try:
								opt_value = int(opt_value)
								if opt_value < 1 or opt_value > 5 or opt_value == None:
									raise ValueError
							except ValueError:
								print "Invalid setting for plot_type option"
								continue
							settings['plot_type'] = opt_value
							showMultiPlots(fids, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings)
						elif opt_name == 'plot_width' or opt_name == 'pw':
							try:
								opt_value = int(opt_value)
								if opt_value < 1 or opt_value == None:
									raise ValueError
							except ValueError:
								print "Invalid setting for plot_width option"
								continue
							settings['plot_width'] = opt_value
							for fid in fids:
								close(fid[0])
							showMultiPlots(fids, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings)
						elif opt_name == 'plot_height' or opt_name == 'ph':
							try:
								opt_value = int(opt_value)
								if opt_value < 1 or opt_value == None:
									raise ValueError
							except ValueError:
								print "Invalid setting for plot_height option"
								continue
							settings['plot_height'] = opt_value
							for fid in fids:
								close(fid[0])
							showMultiPlots(fids, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings)
						elif opt_name == 'flux_bottom' or opt_name == 'fb':
							try:
								opt_value = int(opt_value)
								if opt_value < 1 or opt_value > 3 or opt_value == None:
									raise ValueError
							except ValueError:
								print "Invalid setting for flux_bottom option"
								continue
							settings['flux_bottom'] = opt_value
							showMultiPlots(fids, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings)
						elif opt_name == 'tick_type' or opt_name == 'tt':
							try:
								opt_value = int(opt_value)
								if opt_value < 1 or opt_value > 2 or opt_value == None:
									raise ValueError
							except ValueError:
								print "Invalid setting for tick_type option"
								continue
							settings['tick_type'] = opt_value
							showMultiPlots(fids, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings)
						elif opt_name == 'group_delta' or opt_name == 'gd':
							try:
								opt_value = float(opt_value)
								if opt_value <= 0.0 or opt_value > 100 or opt_value == None:
									raise ValueError
							except ValueError:
								print "Invalid setting for group_delta option"
								continue
							settings['group_delta'] = opt_value
							showMultiPlots(fids, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings)
						elif opt_name == 'tick_label' or opt_name == 'tl':
							try:
								opt_value = int(opt_value)
								if opt_value < 1 or opt_value > 2 or opt_value == None:
									raise ValueError
							except ValueError:
								print "Invalid setting for group_delta option"
								continue
							settings['tick_label'] = opt_value
							showMultiPlots(fids, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings)
						else:
							selOptionHelp()
					else:
						print "Invalid option setting command"
						continue
				else:
					if re.match(".*open.*pod.*bay.*doors.*", lower(s)):
						print "I'm sorry, %s. I'm afraid I can't do that." % os.getenv('USER')
					try:
						nid = int(s)
						if (nid > len(rftList_new) or nid < 1):
							raise ValueError
						fids[-1][1] = nid
						showPlot(rftList_old[fids[-1][1]-1], rftList_new[fids[-1][1]-1], c_old, c_new, colour_config, tick_config, settings, figureid = fids[-1][0])
					except ValueError:
						print "Enter a valid fit id or command"
		except KeyboardInterrupt:
			while 1:
				s = lower(strip(raw_input("Are you sure you want to quit? y/[n]: ")))
				if s == 'y':
					print "Quitting fitcmp."
					sys.exit()
				elif s == '' or s == 'n':
					break

def showMultiPlots(fidlist, rft_old, rft_new, comp_old, comp_new, colour_config, tick_config, settings, show = True, saveFile = ''):
	for fid in fidlist:
		showPlot(rft_old[fid[1]-1], rft_new[fid[1]-1], comp_old, comp_new, colour_config, tick_config, settings, fid[0], show, saveFile)

def exportAll(fname, rftList_old, rftList_new, c_old, c_new, colour_config, tick_config, settings, ext = 'ps'):
	print "--- Exporting %s files" % ext
	for i in range(len(rftList_new)):
		fn = '%s_%s.%s' % (fname, str(i+1).zfill(len(str(len(rftList_new)))), ext)
		showPlot(rftList_old[i], rftList_new[i], c_old, c_new, colour_config, tick_config, settings, show = False, saveFile = fn)

def printf13Paths(oldf, newf):
	print "Old fort.13 path: %s" % abspath(oldf)
	print "New fort.13 path: %s" % abspath(newf)

def reloadOneFit(rfid, f13, rftList, vppath, vpfix, pr = None):
	#print "--- Reloading region data"
	#print "---"
	tmpRList, c = readFort13(f13, pr)
	for i in range(len(tmpRList)):
		rftList[i][RFT_R] = tmpRList[i]
	#print "--- Reloading fit data"
	#print "---"
	rftList[rfid-1][RFT_F], rftList[rfid-1][RFT_T] = vpf17ReadAll(f13, vppath, len(tmpRList), vpfix, pr, rfid)[0]
	return rftList, c

def reloadFit(rfid, f13, rftList, vppath, vpfix, silent = False):
	if silent:
		pr = None
	else:
		pr = ProgressWrapper(term, 'Reloading data for fit %d' % (rfid), PR_F13READ+PR_F17GENERATE+PR_F17READ)
	rftList, c = reloadOneFit(rfid, f13, rftList, vppath, vpfix, pr)
	if not silent:
		pr.clear()
	return rftList, c

def reloadCFits(rfid, oldf, rftList_old, newf, rftList_new, vppath, vpfix):
	pr = ProgressWrapper(term, 'Reloading old and new data for fit %d' % (rfid), (PR_F13READ+PR_F17GENERATE+PR_F17READ)*2)
	rftList_old, c_old = reloadOneFit(rfid, oldf, rftList_old, vppath, vpfix, pr)
	rftList_new, c_new = reloadOneFit(rfid, newf, rftList_new, vppath, vpfix, pr)
	pr.clear()
	return rftList_old, c_old, rftList_new, c_new

def reloadAllFits(f, rftList, vppath, vpfix, silent = False):
	if silent:
		pr = None
	else:
		pr = ProgressWrapper(term, 'Retrieving data for all fits', (PR_F13READ+PR_F17GENERATE+PR_F17READ))
	rftList, c = vpf17AutomateAll(f, vppath, vpfix, pr)
	if not silent:
		pr.clear()
	return rftList, c

def selLines(oldf, newf, rft_old, rft_new, c_old, c_new):
	print "Old file (%s)" % oldf
	getLines(rft_old, c_old)
	print
	print "New file (%s)" % newf
	getLines(rft_new, c_new)

def getLines(rft, c):
	cids = tExtract(rft, T_COM)
	sl = [[strip(c[i][C_L]), i] for i in range(len(c)) if (i+1 in cids)]
	for i in range(len(sl)):
		print ljust("<%d>" % (sl[i][1]+1), 9) + sl[i][0]

def readAllf17(oldf, newf, vppath, vpfix):
	pr = ProgressWrapper(term, 'Retrieving data for all fits', (PR_F13READ+PR_F17GENERATE+PR_F17READ)*2)
	
	#print "--- Reading %s and retrieving fort.17 data for all fits..." % oldf
	#print "---"
	rftList_old, c_old = vpf17AutomateAll(oldf, vppath, vpfix, pr)
	
	#print "--- Reading %s and retrieving fort.17 data for all fits..." % newf
	#print "---"
	rftList_new, c_new = vpf17AutomateAll(newf, vppath, vpfix, pr)
	
	#print "--- Done."
	
	pr.clear()
	
	return rftList_old, c_old, rftList_new, c_new

def termBold(msg):
	return term.render('${BOLD}%s${NORMAL}') %msg
	
def termWarn(msg):
	return term.render('${YELLOW}%s${NORMAL}') %msg

def selOptionHelp():
	print "Available options and settings:"
	print termBold(" option (shrtver)   setting        description")
	print          " plot_type (pt)     1*             data hist, fit line"
	print          "                    2              data hist, fit hist"
	print          "                    3              data line, fit line"
	print          "                    4              data point w/ err, fit line"
	print          "                    5              data err contour, fit line"
	print          " plot_width (pw)    650*           plot width, in pixels"
	print          " plot_height (ph)   490*           plot height, in pixels"
	print          " flux_bottom (fb)   1*             bottom flux limit is 0"
	print          "                    2              bottom flux limit is min(data or fit) rounded"
	print          "                    3              bottom flux limit is min(data or fit)"
	print          " tick_type (tt)     1*             normal tick marks"
	print          "                    2              group fine/isotopic structure lines"
	print          " tick_label (tl)    1*             smart z label (if available), IDs"
	print          "                    2              IDs only"
	print          " group_delta (gd)   0.07*          grouping delta for the above option in Angstroms"
	print          " * = default setting"

def selHelp():
	print "Commands: (use ; to separate multiple entries)"
	print termBold(" <n>") + "                        - Select <n> fit for display"
	print termBold(" n") + "                          - Display next region"
	print termBold(" p") + "                          - Display previous region"
	print termBold(" l") + "                          - List all regions"
	print termBold(" !") + "                          - Print all lines within region, for new and old fort.13 files"
	print termBold(" ^[c]") + "                       - Create a new figure window or [c]lose all figure windows"
	print termBold(" r") + "                          - Reload all fort.17 data for new and old fort.13 files"
	# print termBold(" r<n>") + "                       - Reload all fort.17 data for <n> fit"
	print termBold(" rc") + "                         - Reload all fort.17 data for current region"
	print termBold(" r[n,o]") + "                     - Reload all fort.17 data for new/old fort.13 file"
	print termBold(" rc[n,o]") + "                    - Reload all fort.17 data for new/old fort.13, current region only"
	print termBold(" a[n,o]") + "                     - Autoreload mode for new/old fort.13, current region"
	print termBold(" v[n,o]/lbl[-sp,]/vmin,vmax") + " - Velocity stack plot around the lines with label lbl,"
	print                "                              vmin/max are opt, sp is a comma-delimited subset of species"
	print termBold(" e/expr") + "                     - Evaluates expression, atom.dat lookup using e.g. {H I 1215}"
	print termBold(" (e)ps/<f>") + "                  - Output all displays to sequentially named <f>_n.[e]ps files"
	print termBold(" pdf/<f>") + "                    - Output all displays to sequentially named <f>_n.pdf files"
	print termBold(" f13") + "                        - Prints the paths of the fort.13 files in use"
	print termBold(" {option,setting}") + "           - Set option. Use {} to see available options"
	print termBold(" ?") + "                          - Show this message"
	print termBold(" q") + "                          - Quit immediately"

def listRegions(rftList):
	fidl = 1
	for rft in rftList:
		print ljust("(%d)" % fidl, 6) + strip(rft[RFT_R][R_L])
		fidl += 1

if __name__ == "__main__":
    main()
