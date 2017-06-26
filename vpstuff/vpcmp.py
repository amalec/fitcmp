from numpy import *
import pylab as p
import os, re
from math import log, sqrt
from string import strip
from matplotlib.ticker import *
from os.path import dirname, abspath, isfile
from glob import glob
from ConfigParser import RawConfigParser
from vpstuff.constants import C
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.widgets import MultiCursor
from matplotlib import rcParams
from vpstuff.vphelper import find_line, find_lines_byspecies, show_error, termBold, termWarn, term, RepresentsInt
from vpstuff.pyvpfit import *
from tempfile import NamedTemporaryFile
from operator import itemgetter
import vpstuff.fcs as fcs
import asciitable, sys

#TODO: Add a feature that counts the number of lines in both new and old regions, and prints a warning if there is a differecne

COMMENT_MARKER = '!'
SEP_MARKER = '*'
F17SEP_MARKER = '-'
NO_MODE, REGION_MODE, COMPONENT_MODE = 'NO_MODE', 'REGION_MODE', 'COMPONENT_MODE'
FIT_MODE, TICK_MODE = 'FIT_MODE', 'TICK_MODE'
FIT_THRESH = 1.0e-99

# index labels
R_WL, R_WH, R_L = 0, 1, 2 # Fit region: low wavelength bound, high wavelength bound, line copy
C_NAME, C_L = 0, 1 # Fit components: ion name, line copy
F_WL, F_DATA, F_ERR, F_FIT = 0, 1, 2, 3 # Output fit: wavelength, data, error, fit data
T_WL, T_SPEC, T_COM = 0, 1, 2 # Tick marks: wavelength, species no. (not reliable/random in vpfit10), component no.
RFT_R, RFT_F, RFT_T = 0, 1, 2 # RFT: Region, Fit, Tick marks

# relative estimated progress
PR_F13READ            = 1
PR_F17GENERATE        = 5
PR_F17READ            = 1

class fort13FormatError(Exception):
	def __init__(self):
		return
	def __str__(self):
		return 'Incorrectly formatted fort.13 file!'

class fort17FormatError(Exception):
	def __init__(self):
		return
	def __str__(self):
		return 'Incorrectly formatted fort.17 file!'

class fort17FileError(Exception):
	def __init__(self):
		return
	def __str__(self):
		return 'Could not find fort.17 or vpfit_chunk*.txt files!'

def sanitiseLine(line):
	return line[0:line.find(COMMENT_MARKER)]

def readFort13(filename, pr = None):
	if pr: pr.setMessage('Reading and parsing %s' % filename)
	fort13f = open(filename, 'r')
	lines = fort13f.readlines()
	mode = NO_MODE
	regionList = []; cmpList = []
	for line in lines:
		dat = line.split()
		if line.strip() == '':
			raise fort13FormatError()
		if (dat[0][0] == COMMENT_MARKER): continue # ! at start of line
		if (dat[0] == SEP_MARKER):
			if (mode == NO_MODE):
				mode = REGION_MODE
				continue
			elif (mode == REGION_MODE):
				mode = COMPONENT_MODE
				continue
			elif (mode == COMPONENT_MODE):
				raise fort13FormatError
		
		if (mode == REGION_MODE):
			try:
				regionList.append([float(dat[2]), float(dat[3]), line]) # R_WL, R_WH, R_L
			except:
				raise fort13FormatError
		
		if (mode == COMPONENT_MODE):
			try:
				cleanDat = sanitiseLine(line).split()
				if (len(cleanDat) == 9):
					# ideally this would handle any no. of spaces in ion name
					cmpList.append([cleanDat[-9] + ' ' + cleanDat[-8], line])
				elif (len(cleanDat) == 8):
					cmpList.append([cleanDat[-8], line]) # C_NAME, C_L
			except IndexError:
				raise fort13FormatError
	
	if pr: pr.addProgress(PR_F13READ)
	return regionList, cmpList
######
def readFort13z(filename, pr = None):
	cleanz = lambda z: ''.join([c for c in z if c in '0123456789.']) # sanitise redshift
	if pr: pr.setMessage('Reading and parsing %s' % filename)
	fort13f = open(filename, 'r')
	lines = fort13f.readlines()
	mode = NO_MODE
	regionList = []; cmpList = []
	for line, i in enumerate(lines):
		dat = line.split()
		if (dat[0][0] == COMMENT_MARKER): continue # ! at start of line
		if dat[0][0] == COMMENT_MARKER:
			dat = dat[1:]
		if (dat[0] == SEP_MARKER):
			if (mode == NO_MODE):
				mode = REGION_MODE
				continue
			elif (mode == REGION_MODE):
				mode = COMPONENT_MODE
				continue
			elif (mode == COMPONENT_MODE):
				raise fort13FormatError
		
		if (mode == REGION_MODE):
			try:
				regionList.append([float(dat[2]), float(dat[3]), line]) # R_WL, R_WH, R_L
			except:
				raise fort13FormatError
		
		if (mode == COMPONENT_MODE):
			try:
				cleanDat = sanitiseLine(line).split()
				if (len(cleanDat) == 9):
					# ideally this would handle any no. of spaces in ion name
					cmpList.append([cleanDat[-9] + ' ' + cleanDat[-8], 
						float(cleanz(cleanDat[-6])), line])
				elif (len(cleanDat) == 8):
					cmpList.append([cleanDat[-8], 
						float(cleanz(cleanDat[-6])), line]) # C_NAME, C_L
			except IndexError:
				raise fort13FormatError
	
	if pr: pr.addProgress(PR_F13READ)
	return regionList, cmpList

####

def readFort17All(pr = None):
	vpchunkf = glob('vpfit_chunk*.txt')
	
	if isfile('fort.17'): # for VPFIT 9.5
		if pr: pr.setMessage('Reading and parsing fort.17 file')
		fort17f = open('fort.17', 'r')
		lines = fort17f.readlines()
		fort17f.close()
	elif vpchunkf: # for VPFIT 10
		if pr: pr.setMessage('Reading and parsing vpfit_chunk*.txt files')
		vpchunkf.sort() # rather important to sort them
		lines = []
		for vpcf in vpchunkf:
			fort17f = open(vpcf, 'r')
			lines.extend(fort17f.readlines())
			fort17f.close()
	else:
		raise fort17FileError
	
	#lines.append(F17SEP_MARKER)
	mode = NO_MODE
	allList = []
	fitList = []; tickList = []
	fitDone, tickDone = False, False
	for i, line in enumerate(lines):
		dat = line.split()
		if dat[0][0] == '!':
			dat = dat[1:]
		#if (dat[0][0] == COMMENT_MARKER): continue # ! at start of line
		if (dat[0][0] == F17SEP_MARKER): # - at start of line
			if (mode == NO_MODE):
				mode = FIT_MODE
				continue
			elif (mode == FIT_MODE):
				mode = TICK_MODE
				fitDone = True
				continue
			elif (mode == TICK_MODE):
				mode = FIT_MODE
				tickDone = True
				continue
		
		if (mode == FIT_MODE):
			try:
				if (abs(float(dat[3])) > FIT_THRESH):
					# F_WL, F_DATA, F_ERR, F_FIT
					fitList.append([float(dat[0]), float(dat[1]), float(dat[2]), float(dat[3])])
			except:
				raise fort17FormatError
		
		if (mode == TICK_MODE):
			try:
				# T_WL, T_SPEC, T_COM
				tickList.append([float(dat[0]), int(dat[1]), int(dat[2])])
			except:
				raise fort17FormatError
		if fitDone and tickDone:
			allList.append([fitList, tickList])
			#print fitList[0][0]
			fitDone, tickDone = False, False
			fitList = []; tickList = []
	allList.append([fitList, tickList]) # don't forget the last ones!
	
	if pr: pr.addProgress(PR_F17READ)
	return allList

def fort17Clean():
	vpchunkf = glob('vpfit_chunk*.txt')
	
	try:
		if isfile('fort.17'): # for VPFIT 9.5
			os.remove('fort.17')
		elif vpchunkf: # for VPFIT 10
			for vpcf in vpchunkf:
				os.remove(vpcf)
	except:
		pass

def vpf17ReadAll(fort13path, vppath, nreg, vpfix, pr = None, regionid = None):
	# fid is not 0 indexed, i.e. first fit id is 1!
	if pr: pr.setMessage('Generating fort.17 file')
	vpin =  '%s\\n' % 'd'
    # if not vpfix(1):
	vpin += '\\n'
	vpin += '\\n'
	vpin += '%s\\n' % fort13path
	
	if regionid:
		regs = [regionid]
	else:
		regs = range(1, nreg+1)
	
	for i in regs:
		if nreg != 1: # only one region...
			vpin += '%s\\n' % 'y' # plot
		vpin += '%d\\n' % i
		if i == 1 or regionid:
			vpin += 'as\\n'
			vpin += '\\n'
			vpin += '/NULL\\n'
		else:
			vpin += '\\n'
	vpin += '%s\\n' % 'n' # plot
	vpin += '%s\\n' % 'n'
	
	cmd = 'printf "%s" | %s > /dev/null' % (vpin, vppath)
	# print cmd
	os.system(cmd)
	f13path = dirname(fort13path)
	if pr: pr.addProgress(PR_F17GENERATE)
	vpf17data = readFort17All(pr)
	fort17Clean()
	
	return vpf17data

def vpf17AutomateAll(fort13path, vppath, vpfix, pr = None):
	r, c = readFort13(fort13path, pr)
	rftList = []
	ft = vpf17ReadAll(fort13path, vppath, len(r), vpfix, pr)
	# print len(r), len(ft)
	assert(len(r) == len(ft))
	for i in range(1,len(r)+1):
		rftList.append([r[i-1], ft[i-1][0], ft[i-1][1]]) # RFT_R, RFT_F, RFT_T
	return rftList, c

def fExtract(rft, F_INDEX):
	""" F_INDEX is either of F_WL, F_DATA, F_ERR, F_FIT """
	somedat  = [dat[F_INDEX] for dat in rft[RFT_F] if (dat[F_WL] >= rft[RFT_R][R_WL] and dat[F_WL] <= rft[RFT_R][R_WH])]
	return somedat

def tExtract(rft, T_INDEX):
	"""T_INDEX is either of T_WL, T_SPEC, T_COM """
	tdat  = [dat[T_INDEX] for dat in rft[RFT_T] if (dat[T_WL] >= rft[RFT_R][R_WL] and dat[T_WL] <= rft[RFT_R][R_WH])]
	return tdat

def rft2Data(rft, wl_flag, data_flag, hist_step = True):
	# extracts a wavelength array from the fort.17 data file that is within the specified region bounds
	wldat  = fExtract(rft, wl_flag)
	datdat = fExtract(rft, data_flag)
	
	if hist_step:
		# make bin arrays
		wlbin = [wldat[0]-(wldat[1]-wldat[0])/2.0, wldat[0]+(wldat[1]-wldat[0])/2.0] # leftcap
	
		for i in range(1,len(wldat)-1): # body
			wlbin.append(wldat[i]-(wldat[i]-wldat[i-1])/2.0)
			wlbin.append(wldat[i]+(wldat[i+1]-wldat[i])/2.0)
		wlbin.append(wldat[len(wldat)-1]-(wldat[len(wldat)-1]-wldat[len(wldat)-2])/2.0)
		wlbin.append(wldat[len(wldat)-1]+(wldat[len(wldat)-1]-wldat[len(wldat)-2])/2.0) # rightcap
	
		datbin = []
		for i in range(0,len(datdat)):
			datbin.extend([datdat[i], datdat[i]])
	
		return wlbin, datbin
	else:
		return wldat, datdat

def readColourConfig(config_file):
	# TODO: check config file for errors
	config = RawConfigParser()
	if not isfile(config_file):
		show_error("Could not open colour config file = %s" % config_file, True)
	config.read(config_file)
	colour_config = {}
	tick_config = []
	for section in config.sections():
		if section != 'plot':
			tick_config.append({'match': config.get(section, 'match'), 'colour': config.get(section, 'colour')})
	for col in config.items('plot'):
		colour_config[col[0]] = col[1]
	return colour_config, tick_config

def Nb_estimate(norm_flux, delta_v, f, lam_0):
	# From J. Whitmore
	# delta_v is between the two successive clicks
	x = 1.52e-5
	N = log( (-log(norm_flux)*delta_v) / (f*lam_0*x), 10)
	b = delta_v / (2.0*sqrt(2.0*log(2.0)))
	return N, b

# TODO: Maybe eventually support stacked common-velocity display
def showStackPlot(tiedz_lbl, rft_all, comp, colour_config, tick_config, settings, splist, vlow = None, vhigh = None, cursor_on = False, saveFile = None):
	if settings['crs_display'] == 1:
		addplot = 1
		crsplus = 'CRS + '
	else:
		addplot = 0
		crsplus = ''
		
	pc = parseComps(comp)
	sel_comps = []
	
	if tiedz_lbl[:3].lower() == 'red':
		redshift = float(tiedz_lbl[3:])
		for sp in splist:
			sel_comps.append([sp, redshift])
	else:
	
		for i, pci in enumerate(pc):
			lbl = pci[5]
			species = pci[0]
			if lbl.lower() == tiedz_lbl.lower() and (species in splist or splist == []):
				redshift = float(pci[4])
				sel_comps.append([species, redshift])
			if RepresentsInt(tiedz_lbl):
				if int(tiedz_lbl) == i+1 and (species in splist or splist == []):
					redshift = float(pci[4])
					sel_comps.append([species, redshift])
			# print species
			
	
	found_lines = []
	if sel_comps:
		for sc in sel_comps:
			sp_lines = find_lines_byspecies(sc[0])
			for sp in sp_lines:
				found_lines.append([sc[0], sc[1], sp]) # species, redshift, rest wavelength - replaced with fl result
		
		# group 
		if settings['tick_type'] != 1:
			ingroups = []
			notgroups = []
			groupdelta = settings['group_delta']
			
			for sc in sel_comps: # in failed attempts at implementing this I was looping over found_lines, which obviously generated many duplicates
				newgroup = False
				gfl = filter(lambda g: g[0] == sc[0], found_lines) # 1. find other lines with matching species
				gfl = sorted(gfl, key=lambda g: float(g[2]['wv'])) # 2. sort by wavelength
				if len(gfl) > 1:
					septotal = 0.0
					group = []
					for i in range(1, len(gfl)):
						separation = float(gfl[i][2]['wv'])-float(gfl[i-1][2]['wv']) # abs is omitted here to make the algorithm exclude filling the result arrays with duplicates
						
						if separation + septotal <= groupdelta:
							if not newgroup: newgroup = True
							septotal += separation
							group.append(gfl[i-1])
							if i == len(gfl)-1:
								group.append(gfl[i])
								ingroups.append(group)
						else:
							if newgroup:
								group.append(gfl[i-1])
								ingroups.append(group)
								group = []
								newgroup = False
								septotal = 0.0
							else:
								notgroups.append(gfl[i-1])
							if i == len(gfl)-1:
								notgroups.append(gfl[i])
				else:
					notgroups.append(gfl[0])
			
			# now weigh the grouped lines if there are any
			
			if ingroups:
				for ig in ingroups:
					# there's one more layer.. the group within the list of groups
					f_sum = sum([float(fl[2]['f']) for fl in ig])
					g_weighted_wl = sum([float(fl[2]['f'])*float(fl[2]['wv']) for fl in ig]) / f_sum
					weighted_fl = ig[0] # assuming that grouped lines share the same properties
					weighted_fl[2]['f'] = str(f_sum) # str for consistency
					weighted_fl[2]['wv'] = str(g_weighted_wl)
					notgroups.append(weighted_fl)
			
			found_lines = notgroups
		
		temp_fl = []
		# effectively remove lines that don't fall into the fort.13 file regions, and identify regions where they do
		for rindx, rft in enumerate(rft_all):
			for fl in found_lines:
				wv_obs = float(fl[2]['wv'])*(fl[1]+1.0)
				if wv_obs <= rft[RFT_R][R_WH] and wv_obs >= rft[RFT_R][R_WL]:
				# if wv_obs <= rft[RFT_R][R_WH]+10.0 and wv_obs >= rft[RFT_R][R_WL]-10.0:
					temp_fl.append([fl[0], fl[1], fl[2], rindx]) # 0 species, 1 redshift, 2 rest wavelength, 3 region index
		found_lines = sorted(temp_fl, key=itemgetter(3)) # used to be 2, 3 works better
		
		
		if len(found_lines) == 0:
			print "No lines with label '%s' found in current regions" % tiedz_lbl
			return
		
		#######################################################
		# extract data and plot it
		
		if not saveFile:
			p.close('all')
		pdpi = rcParams['figure.dpi']
		
		fig, ax = p.subplots(len(found_lines)+addplot, 1, sharex=True, figsize=(float(settings['vplot_width'])/pdpi, float(settings['vplot_height'])/pdpi)) #default: 8*80, 13*80
		fig.canvas.set_window_title(crsplus + 'Velocity stack')
		
		class TwoClick:
			def __init__(self, figure, settings):
				self.second_click = False
				self.figure = figure
				self.cid_onclick = self.figure.canvas.mpl_connect('button_press_event', self.onclick)
				self.settings = settings
				self.clicks = {'x1': None, 'y1': None, 'x2': None, 'y2': None}
			
			def onclick(self, event):
				try:
					if event.button == 3: # right click
						fl_indx = event.inaxes.get_geometry()[2]-2 # based on subplot id
						if self.settings['crs_display'] != 1:
							fl_indx += 1
						
						if event.ydata <= 0.0 or event.ydata >= 1.0:
							return
						
						if not self.second_click:
							self.clicks['x1'], self.clicks['y1'] = event.xdata, event.ydata
						else:
							self.clicks['x2'], self.clicks['y2'] = event.xdata, event.ydata
						
						if self.second_click:
							# do the calc
							norm_flux = self.clicks['y1']
							delta_v = abs(self.clicks['x1']-self.clicks['x2'])
							if not delta_v <= 0.0:
								f = float(found_lines[fl_indx][2]['f'])
								lam_0 = float(found_lines[fl_indx][2]['wv'])*1.0e-10 # in metres
								wv_obs0 = float(found_lines[fl_indx][2]['wv'])*(float(found_lines[fl_indx][1])+1.0)
								wv_obs = self.clicks['x1']/C*1000.0*wv_obs0+wv_obs0
								z = wv_obs/float(found_lines[fl_indx][2]['wv'])-1.0
								N, b = Nb_estimate(norm_flux, delta_v, f, lam_0)
								print "   %s %s   %s   %s       0.000SE      0.00   0.00E+00  0  !     " % (found_lines[fl_indx][0].ljust(8), ("%.5f" % N).rjust(8), ("%.7f" % z).rjust(11), ("%.4f" % b).rjust(9))
							else:
								self.second_click = False
								return
								
						self.second_click = not self.second_click
						# 
				except TypeError:
					pass
		
		class KeyAction:
			def __init__(self, figure):
				self.figure = figure
				self.cid_onrelease = self.figure.canvas.mpl_connect('key_release_event', self.onrelease)
			
			def onrelease(self, event):
				if event.key == ',':
					print "Setting left bound to %.2f" % event.xdata
					self.figure.get_axes()[-1].set_xlim(left=event.xdata)
					self.figure.canvas.draw()
				if event.key == '.':
					print "Setting right bound to %.2f" % event.xdata
					self.figure.get_axes()[-1].set_xlim(right=event.xdata)
					self.figure.canvas.draw()
		
		tc = TwoClick(fig, settings)
		ka = KeyAction(fig)
		
		vall = []
		
		for i, fl in enumerate(found_lines):
			rft = rft_all[fl[3]]
			wl_raw, dat_raw = rft2Data(rft, F_WL, F_DATA, False)
			wlbin, datbin = rft2Data(rft, F_WL, F_DATA, settings['plot_type'] <= 2)
			wldat, fitdat = rft2Data(rft, F_WL, F_FIT, settings['plot_type'] == 2)
			daterr = fExtract(rft, F_ERR)
			# wldat = fExtract(rft, F_WL)
			fitdat_raw = fExtract(rft, F_FIT)
			twl = tExtract(rft, T_WL)
			tcom = tExtract(rft, T_COM)
			tsp = tExtract(rft, T_SPEC)
			
			wv_obs = float(fl[2]['wv'])*(fl[1]+1.0)
			vel_raw = [(w-wv_obs)/wv_obs*C/1000.0 for w in wl_raw]
			vbin = [(w-wv_obs)/wv_obs*C/1000.0 for w in wlbin] # in km/s
			vdat = [(w-wv_obs)/wv_obs*C/1000.0 for w in wldat] # in km/s
			tv = [(w-wv_obs)/wv_obs*C/1000.0 for w in twl] # in km/s
			vdelta = abs(vel_raw[0]-vel_raw[1])
			
			vdata = {
				'vel_raw': vel_raw, 
				'dat_raw': dat_raw, 
				'fitdat_raw': fitdat_raw, 
				'vbin': vbin, 
				'vdat': vdat, 
				'tv': tv, 
				'datbin': datbin, 
				'daterr': daterr, 
				'fitdat': fitdat, 
				'twl': twl, 
				'tcom': tcom, 
				'tsp': tsp, 
				'vdelta': vdelta
				# 'ylabel': "%s_%i" % (fl[0], int(float(fl[2]['wv']))),
			}
			vall.append(vdata)
			velocityPlot(ax[i+addplot], vdata, pc, fl, colour_config, tick_config, settings) ###### velocityPlot
			ax[i+addplot].set_ylabel("%s %i" % (fl[0], int(float(fl[2]['wv']))), stretch='extra-condensed')
		
		minvel = min([min(vi['vel_raw']) for vi in vall])
		maxvel = max([max(vi['vel_raw']) for vi in vall])
		
		fig.subplots_adjust(hspace = 0.00, top = 0.98, bottom = 0.05, right = 0.97, left = 0.10)
		if cursor_on:
			multi = MultiCursor(fig.canvas, ax, color='r', lw=1)
		for a in ax[:1]:
			hideXLabels(a)
		ax[-1].set_xlabel('Velocity (km/s) [z = %.7f]' % found_lines[0][1])
		packet = {
			'z': [found_lines[0][1]]
		}
		dd_writer(packet, 'vel.z')
		
		if vlow != None and vhigh != None:
			ax[-1].set_xlim(vlow, vhigh)
		else:
			ax[-1].set_xlim(minvel, maxvel)
		
		if settings['crs_display'] == 1:
			# calculate the CRS, use the bluemost pixel as the reference value for the grid
		
			# zerovel_indxs = [findClosest(0.0, vi['vel_raw']) for vi in vall]
			# zerovel = [vall[i]['vel_raw'][zi] for i, zi in enumerate(zerovel_indxs)]
			# print zerovel
			vdel = sum([vi['vdelta'] for vi in vall])/len(vall)
		
			numbins = int((maxvel-minvel)/vdel)
			crs_vel = [minvel+i*vdel for i in range(numbins)]
			# crs_mask = []
			crs_res = []
			avgNoise = lambda l: sum(l)/sqrt(len(l))
			# avgRes = lambda l: sum(l)/len(l)
		
			for cv in crs_vel: 
				vel_indxs = [findClosest(cv, vi['vel_raw']) for vi in vall] # find closest bin crs_vel <-> vel_raw stacks
			
				res_vals = []
				for j, vel_i in enumerate(vel_indxs):
					if abs(cv - vall[j]['vel_raw'][vel_i]) < vdel and vall[j]['daterr'][vel_i] > 0: # if within a bin of cv and a valid pixel
						residual = calcResidual(vall[j]['dat_raw'][vel_i], vall[j]['fitdat_raw'][vel_i], vall[j]['daterr'][vel_i])
						res_vals.append(residual)
			
				if len(res_vals) > 0:
					# crs_mask.append(1)
					crs_res.append(avgNoise(res_vals))
				else:
					# crs_mask.append(0)
					crs_res.append(NaN)
		
			ax[0].axhline(1.0, c=colour_config['res_zero_one'])
			ax[0].axhline(-1.0, c=colour_config['res_zero_one'])

			# ax[0].plot(wldat_old, res_old, c=colour_config['res_zero_one'], linestyle='--')
			ax[0].plot(crs_vel, crs_res, c=colour_config['residual'])
			packet = {
				'vel': crs_vel,
				'res': crs_res
			}
			dd_writer(packet, 'vel.crs')
		
		if saveFile:
			print "Saving velocity stack plot to %s" % (saveFile)
			p.savefig(saveFile)
			p.close()
		else:
			p.ioff()
			print termWarn("[Close display window to continue]")
			print "Select left and right bounds using <,> and <.> keys"
			p.show()
			p.ion()
	else:
		print "No lines found with label '%s'" % tiedz_lbl
		return

def findClosest(targetVal, valList):
	""" Searches valList for closest match to targetVal and returns the
	    corresponding valList index"""
	diffs = [abs(x-targetVal) for x in valList]
	return diffs.index(min(diffs))

def dd_writer(packet, fpart):
	if fcs.dd:
		dd_fname = "%s.%s.dat" % (fcs.dd_prefix, fpart)
		if fcs.dd_test:
			print "@FILENAME: " + dd_fname
			asciitable.write(packet, sys.stdout)
		else:
			print "DD: " + dd_fname
			asciitable.write(packet, dd_fname)


def velocityPlot(ax, data, pc, fline, colour_config, tick_config, settings):
	ax.axhline(1.0, c=colour_config['zero_one'], linestyle = ':')
	ax.axhline(0.0, c=colour_config['zero_one'], linestyle = ':')
	
	dd_fpart = "%s_%.4f" % (''.join(ff for ff in fline[0] if ff.isalpha()), float(fline[2]['wv']))
	
	if settings['plot_type'] == 4:
		ax.errorbar(data['vbin'], data['datbin'], yerr=data['daterr'], color=colour_config['data'], fmt=None)
		packet = {
			'vbin': data['vbin'],
			'datbin': data['datbin'],
			'daterr': data['daterr']
		}
		dd_writer(packet, dd_fpart+'.vel.data.errbar')

			
	elif settings['plot_type'] == 5:
		filla = [d-e for d, e in zip(data['datbin'], data['daterr'])]
		fillb = [d+e for d, e in zip(data['datbin'], data['daterr'])]
		ax.fill_between(data['vbin'], filla, y2 = fillb, lw=0.0, color=colour_config['data_contour'])
		packet = {
			'vbin': data['vbin'],
			'filla': filla,
			'fillb': fillb
		}
		dd_writer(packet, dd_fpart+'.vel.data.contour')
	else:
		ax.plot(data['vbin'], data['datbin'], color=colour_config['data'])
		packet = {
			'vbin': data['vbin'],
			'datbin': data['datbin']
		}
		dd_writer(packet, dd_fpart+'.vel.data')
		
	ax.plot(data['vdat'], data['fitdat'], color=colour_config['fit_new'])
	packet = {
		'vdat': data['vdat'],
		'fitdat': data['fitdat']
	}
	dd_writer(packet, dd_fpart+'.vel.fit')
	
	
	# vdata = {
	# 	'vel_raw': vel_raw,
	# 	'dat_raw': dat_raw,
	# 	'fitdat_raw': fitdat_raw,
	# 	'vbin': vbin,
	# 	'vdat': vdat,
	# 	'tv': tv,
	# 	'datbin': datbin,
	# 	'daterr': daterr,
	# 	'fitdat': fitdat,
	# 	'twl': twl,
	# 	'tcom': tcom,
	# 	'tsp': tsp,
	# 	'vdelta': vdelta
	# 	# 'ylabel': "%s_%i" % (fl[0], int(float(fl[2]['wv']))),
	# }
	# vall.append(vdata)
	# velocityPlot(ax[i+addplot], vdata, pc, fl, colour_config, tick_config, settings) ###### velocityPlot
	# ax[i+addplot].set_ylabel("%s %i" % (fl[0], int(float(fl[2]['wv']))), stretch='extra-condensed')
		
	# set view bounds
	dax = ax.axis() # xmin, xmax, ymin, ymax
	# ymax = dax[3]
	
	if settings['flux_bottom'] == 1:
		if (dax[2] >= 0.0):
			ymin = 0.0 - 0.1
		else:
			ymin = dax[2]
		ymin -= 0.05
		ax.yaxis.set_major_locator(FixedLocator([0.0, 1.0]))
	elif settings['flux_bottom'] == 2:
		ymin = round(min(data['datbin']+data['fitdat'])-0.06, 1)
		ymin -= 0.05
		ax.yaxis.set_major_locator(FixedLocator([ymin, 1.0]))
	elif settings['flux_bottom'] == 3:
		tmpmin = min(data['datbin']+data['fitdat'])
		tmpmax = max(data['datbin']+data['fitdat'])
		ymin = tmpmin
		ymax = tmpmax
		ax.yaxis.set_major_locator(FixedLocator([ymin, ymax]))
	
	if settings['vel_res'] == 2:
		YMAX_OFFSET = 0.55
		RESIDUAL_OFFSET = 0.3
		RESIDUAL_SCALE = 0.05
	else:
		YMAX_OFFSET = 0.2
	
	TICK_OFFSET = 0.05
	TICK_SCALE = 0.05
	
	ymax = 1.00 + YMAX_OFFSET*(1.0-ymin)
	
	ax.axis([dax[0], dax[1], ymin, ymax])
	
	if settings['vel_res'] == 2:
		resy = 1.00 + RESIDUAL_OFFSET*(1.0-ymin)
	
	##
	
	# tick positioning
	y0y1 = [1.00 + TICK_SCALE*(ymax-ymin), 1.00 - TICK_SCALE*(ymax-ymin)]
	ytxt = None # don't plot text labels on velocity stack plot	
	
	if settings['tick_type'] == 1:
		drawTicks(data['twl'], data['tcom'], data['tsp'], pc, ax, tick_config, settings, y0y1, ytxt, fline = fline, fpart = dd_fpart)
	elif settings['tick_type'] == 2:
		drawGroupedTicks(data['twl'], data['tcom'], data['tsp'], pc, ax, tick_config, settings, y0y1, ytxt, fline = fline, fpart = dd_fpart)
	elif settings['tick_type'] == 3:
		drawGroupedTicks(data['twl'], data['tcom'], data['tsp'], pc, ax, tick_config, settings, y0y1, ytxt, fline = fline, weighted = True, fpart = dd_fpart)
	
	if settings['vel_res'] == 2:
		ax.axhline(resy+RESIDUAL_SCALE*(ymax-ymin), c=colour_config['res_zero_one'])
		ax.axhline(resy-RESIDUAL_SCALE*(ymax-ymin), c=colour_config['res_zero_one'])
		residual = [calcResidual(data['dat_raw'][vel_i], data['fitdat_raw'][vel_i], data['daterr'][vel_i]) for vel_i in range(len(data['daterr']))]
		proj_res = [r*RESIDUAL_SCALE*(ymax-ymin)+resy for r in residual] # plot coordinates projected residual
		ax.plot(data['vdat'], proj_res, c=colour_config['residual'])
		packet = {
			'vdat': data['vdat'],
			'residual': residual, #included for reference
			'proj_res': proj_res
		}
		dd_writer(packet, dd_fpart+'vel.res')
	
	ax.yaxis.set_minor_locator(FixedLocator([0.25, 0.5, 0.75]))
	if settings['flux_bottom'] == 3:
		ax.yaxis.set_major_formatter(FormatStrFormatter("%.1e"))
	# else:
	# 	ax.yaxis.set_major_formatter(FormatStrFormatter("%0.1f"))
	ax.yaxis.set_minor_formatter(NullFormatter())
	

def setPlotBounds(axes, settings, fitdat_old, fitdat_new, datbin, wlbin):
	# set view bounds
	daxis = axes.axis() # xmin, xmax, ymin, ymax
	ymax = daxis[3]
	
	if settings['flux_bottom'] == 1:
		if (daxis[2] >= 0.0):
			ymin = 0.0 - 0.1
		else:
			ymin = daxis[2]
		ymin -= 0.05
		axes.yaxis.set_major_locator(FixedLocator([0.0, 1.0]))
	elif settings['flux_bottom'] == 2:
		ymin = round(min(fitdat_old+fitdat_new+datbin)-0.06, 1)
		ymin -= 0.05
		axes.yaxis.set_major_locator(FixedLocator([ymin, 1.0]))
	elif settings['flux_bottom'] == 3:
		tmpmin = min(fitdat_old+fitdat_new+datbin)
		tmpmax = max(fitdat_old+fitdat_new+datbin)
		ymax = tmpmax
		ymin = tmpmin
		axes.yaxis.set_major_locator(FixedLocator([ymin, ymax]))
	
	if settings['flux_bottom'] == 3:
		axes.yaxis.set_major_formatter(FormatStrFormatter("%.2e"))
	# else:
	# 	ax.yaxis.set_major_formatter(FormatStrFormatter("%0.1f"))
	
	xmax = max(wlbin)
	xmin = min(wlbin)

	axes.axis([xmin, xmax, ymin, ymax])

def dumpFits(fname, rftList, comps):
	pc = parseComps(comps)
	# print pc
	#['CrII', '13.00818', '', False, '2.3090528', 'ac', True, '3.9817', 'ac', True, '-0.069', '', False, 'QA', '0.00', 'E+00', 1, '   CrII     13.00818    2.3090528ac   3.9817ac   -0.069QA     0.00  1.00E+00  0 ! 2\n']
	# species0, N1,lbl2,small3, z4,lbl5,small6, b7,lbl8,small9, e10,lbl11,small12, bturb13, temp14, rgnflg15, id16, linecopy17
	for i, rft in enumerate(rftList):
		wldat = fExtract(rft, F_WL)
		fitdat = fExtract(rft, F_FIT)
		errdat = fExtract(rft, F_ERR)
		datdat = fExtract(rft, F_DATA)
		rlc = rft[RFT_R][R_L]
		
		assert(len(wldat) == len(fitdat) == len(errdat) == len(datdat))
		
		f = open("%s_%i.fit.dat" % (fname, i+1),'w')
		
		f.write("# " + rlc.strip() + "\n")
		f.write("# wl\tdata\tfit\terror\n")

		for j in range(len(wldat)):
			f.write("%s\t%s\t%s\t%s\n" % (wldat[j], fitdat[j], errdat[j], datdat[j]))
		
		f.close()
		
		
		twl = tExtract(rft, T_WL)
		tcom = tExtract(rft, T_COM)
		tsp = tExtract(rft, T_SPEC)
		
		assert(len(twl) == len(tcom) == len(tsp))
		
		f = open("%s_%i.comp.dat" % (fname, i+1),'w')
		
		f.write("# wl\tline no\tspecies\tN\tNlbl\tz\tzlbl\tb\tblbl\n")
		
		for j in range(len(twl)):
			ti = tcom[j]-1
			f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (twl[j], tcom[j], pc[ti][0], pc[ti][1], pc[ti][2], pc[ti][4], pc[ti][5], pc[ti][7], pc[ti][8]))
			
		f.close()
		
		# print twl, tcom, tsp
		# [6836.9024, 6837.04517, 6837.20158, 6837.31171, 6837.59291, 6837.94271] [1, 2, 3, 4, 5, 6] [1, 1, 1, 1, 1, 1]
		
		
		
		
		
		
		


def plotData(rft_old, rft_new, comp_old, comp_new, axes, settings, colour_config, live):
	wlbin, datbin = rft2Data(rft_new, F_WL, F_DATA, settings['plot_type'] <= 2)
	wldat_old, fitdat_old = rft2Data(rft_old, F_WL, F_FIT, settings['plot_type'] == 2)
	wldat_new, fitdat_new = rft2Data(rft_new, F_WL, F_FIT, settings['plot_type'] == 2)
	daterr = fExtract(rft_new, F_ERR)
	# wldat_old = fExtract(rft_old, F_WLL)
	# wldat_new = fExtract(rft_new, F_WL)
	# fitdat_old = fExtract(rft_old, F_FIT)
	# fitdat_new = fExtract(rft_new, F_FIT)
	
	axes.axhline(1.0, c=colour_config['zero_one'], linestyle = ':')
	axes.axhline(0.0, c=colour_config['zero_one'], linestyle = ':')
	
	if settings['plot_type'] == 4:
		axes.errorbar(wlbin, datbin, yerr=daterr, color=colour_config['data'], fmt=None)
		packet = {
			'wlbin': wlbin,
			'datbin': datbin,
			'daterr': daterr
		}
		dd_writer(packet, 'data.errbar')
	elif settings['plot_type'] == 5:
		filla = [d-e for d, e in zip(datbin, daterr)]
		fillb = [d+e for d, e in zip(datbin, daterr)]
		axes.fill_between(wlbin, filla, y2 = fillb, lw=0.0, color=colour_config['data_contour'])
		packet = {
			'wlbin': wlbin,
			'filla': filla,
			'fillb': fillb
		}
		dd_writer(packet, 'data.contour')
	else:
		axes.plot(wlbin, datbin, color=colour_config['data'])
		packet = {
			'wlbin': wlbin,
			'datbin': datbin
		}
		dd_writer(packet, 'data')
	
	new_color = colour_config['fit_new']
	if live: new_color = colour_config['fit_live']
	
	df = settings['decompose_fit']
	if df > 1:
		if df == 2:
			comp = comp_new
			rft = rft_new
		if df == 3:
			comp = comp_old
			rft = rft_old
		pc = parseComps(comp)
		rlc = rft[RFT_R][R_L]
		rwl = rft[RFT_R][R_WL]
		rwh = rft[RFT_R][R_WH]
		twl = tExtract(rft, T_WL)
		tcom = tExtract(rft, T_COM)
		tcom_unique = uniqueList(tcom)
		
		# space reserved for some component info pre-processing (e.g. for summed column density lines and tied parameters)
		# ...
		
		# remove lines not falling in the current region
		pc = [c for c in pc if c[16]+1 in tcom_unique]
		# remove misc lines 
		pc = [c for c in pc if c[0] not in ['<<', '>>', '<>', '__']]
		
		# retrieve fit data
		decomp_fits = []
		for pci in pc:
			tempf13 = NamedTemporaryFile(prefix = 'fort.13', dir = './')
			tempf13.write('   *\n')
			tempf13.write(rlc)
			tempf13.write('   *\n')
			tempf13.write(pci[17]) #linecopy of component
			tempf13.flush()
			rftf = vpf17ReadAll(os.path.basename(tempf13.name), settings['vppath'], 1, None, regionid = 1)[0][0] # essentially RFT_F
			fiti, wli = [di[F_FIT] for di in rftf if di[F_FIT] < settings['decomp_thresh']], [di[F_WL] for di in rftf if di[F_FIT] < settings['decomp_thresh']] # quick and dirty, might not work well for regions with continuum adjustments - though these have been pruned out, so it will be fine
			decomp_fits.append([wli, fiti]) # note swap
			tempf13.close()
		
		# 
		# R_WL, R_WH, R_L = 0, 1, 2 # Fit region: low wavelength bound, high wavelength bound, line copy
		# C_NAME, C_L = 0, 1 # Fit components: ion name, line copy
		# F_WL, F_DATA, F_ERR, F_FIT = 0, 1, 2, 3 # Output fit: wavelength, data, error, fit data
		# T_WL, T_SPEC, T_COM = 0, 1, 2 # Tick marks: wavelength, species no. (not reliable/random in vpfit10), component no.
		# RFT_R, RFT_F, RFT_T = 0, 1, 2 # RFT: Region, Fit, Tick marks
		# # species0, N1,lbl2,small3, z4,lbl5,small6, b7,lbl8,small9, e10,lbl11,small12, bturb13, temp14, rgnflg15, id16, linecopy17
		# pci = [p0[0], p0[1], p0[2]+p0[3], p0[2] != '', p0[4], p0[5]+p0[6], p0[5] != '', p0[7], p0[8]+p0[9], p0[8] != '', p0[10], p0[11]+p0[12], p0[11] != '', p0[13], p0[14], p0[16], i, c[C_L]]
		
	if df == 3:
		for counter, fi in enumerate(decomp_fits):
			axes.plot(fi[0], fi[1], color=colour_config['fit_old'], linestyle = '--')
			packet = {
				'wlbin': fi[0],
				'decomp_fit': fi[1]
			}
			dd_fpart = "%04i" % counter+1
			dd_writer(packet, dd_fpart+'.decompfit.old')
	else:
		axes.plot(wldat_old, fitdat_old, color=colour_config['fit_old'], linestyle = '--')
		packet = {
			'wldat': wldat_old,
			'fitdat': fitdat_old
		}
		dd_writer(packet, 'fit.old')
		
	if df == 2:
		for counter, fi in enumerate(decomp_fits):
			axes.plot(fi[0], fi[1], color=new_color)
			packet = {
				'wlbin': fi[0],
				'decomp_fit': fi[1]
			}
			dd_fpart = "%04i" % counter+1
			dd_writer(packet, dd_fpart+'.decompfit.new')
	else:
		axes.plot(wldat_new, fitdat_new, color=new_color)
		packet = {
			'wldat': wldat_new,
			'fitdat': fitdat_new
		}
		dd_writer(packet, 'fit.new')
	
	setPlotBounds(axes, settings, fitdat_old, fitdat_new, datbin, wlbin)
	
	if settings['flux_bottom'] == 2 or settings['flux_bottom'] == 3:
		labels = axes.get_yticklabels()
		for label in labels:
			label.set_rotation(90)

def showPlot(rft_old, rft_new, comp_old, comp_new, colour_config, tick_config, settings, figureid = 1, show = True, saveFile = '', live = False):
	assert(len(rft_old) == len(rft_new))
	pdpi = rcParams['figure.dpi']
	fig = p.figure(figureid, figsize=(float(settings['plot_width'])/pdpi, float(settings['plot_height'])/pdpi))
	fig.canvas.set_window_title('Fit comparison %i' % (figureid))
	p.clf()
	# fig.canvas.flush_events()
	
	# need to delete the old figure callbacks
	# print fig.canvas.callbacks.callbacks
	# print
	try:
		fig.canvas.callbacks.callbacks['pick_event'] = {}
		del fig.canvas.callbacks.callbacks['button_press_event'][max(fig.canvas.callbacks.callbacks['button_press_event'].keys())]
	except KeyError:
		pass
	# print fig.canvas.callbacks.callbacks
	
	def onclick(event):
		try:
			if event.button == 3: # right click
				print '\nwl=%f, y=%f' % (event.xdata, event.ydata)
		except TypeError:
			pass
	
	def onpick(event):
		if event.mouseevent.button == 1:
			thistext = event.artist
			tlabel = thistext.get_label()
			# ind = event.ind
			print '\n%s' % (tlabel)
	
	
	cid1 = fig.canvas.mpl_connect('pick_event', onpick)
	cid2 = fig.canvas.mpl_connect('button_press_event', onclick)
	
	# print cid1, cid2
	
	# print fig.get_figwidth() 8.125
	# print fig.get_figheight() 6.125
	# print fig.get_dpi() 80
	
	left, width = 0.05, 0.90 + 0.03
	rect1 = [left, 0.8, width, 0.10] # left, bottom, width, height
	rect2 = [left, 0.6, width, 0.2]
	rect3 = [left, 0.1, width, 0.5]
	
	axUpper      = p.axes(rect1)
	axMiddle     = p.axes(rect2, sharex=axUpper)
	axLower      = p.axes(rect3, sharex=axUpper)
	
	axUpper.xaxis.set_major_formatter(NullFormatter())
	hideXLabels(axUpper)
	title = axUpper.set_title(rft_new[RFT_R][R_L], fontsize='smaller')
	title.set_y(1.1)  # move it up a bit higher than the default
	title.set_x(0)  # align the title left, axes coords
	title.set_horizontalalignment('left')  # align the title left, axes coords
	
	axMiddle.xaxis.set_major_formatter(NullFormatter())
	axMiddle.yaxis.set_major_formatter(NullFormatter())
	axMiddle.yaxis.set_major_locator(NullLocator())
	hideXLabels(axMiddle)
	
	axLower.xaxis.set_major_locator(MaxNLocator(10))
	axLower.xaxis.set_minor_locator(MaxNLocator(100))
	axLower.xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))

	plotResidual(rft_old, rft_new, axUpper, colour_config)
	plotTicks(rft_old, rft_new, comp_old, comp_new, axMiddle, settings, tick_config)
	plotData(rft_old, rft_new, comp_old, comp_new, axLower, settings, colour_config, live)
	
	if (show): 
		# p.show()
		fig.canvas.draw()

	if (saveFile != ''): p.savefig(saveFile, orientation='landscape')

def calcResidual(data, fit, error):
	if error > 0.0:
		r = (data-fit)/error
	else:
		r = NaN
	return r

def residualDat(rft):
	datdat = fExtract(rft, F_DATA)
	fitdat = fExtract(rft, F_FIT)
	wldat = fExtract(rft, F_WL)
	errdat = fExtract(rft, F_ERR)	
	assert(len(datdat) == len(fitdat) == len(wldat) == len(errdat))
	#print [(datdat[i], fitdat[i], errdat[i]) for i in range(len(wldat))] 
	resdat = [calcResidual(datdat[i],fitdat[i],errdat[i]) for i in range(len(wldat))]
	return resdat

def plotResidual(rft_old, rft_new, axes, colour_config):
	wldat_old = fExtract(rft_old, F_WL)
	wldat_new = fExtract(rft_new, F_WL)
	res_old = residualDat(rft_old)
	res_new = residualDat(rft_new)
	
	axes.axhline(1.0, c=colour_config['res_zero_one'])
	axes.axhline(-1.0, c=colour_config['res_zero_one'])
	
	axes.plot(wldat_old, res_old, c=colour_config['res_zero_one'], linestyle='--')
	packet = {
		'wldat': wldat_old,
		'res': res_old
	}
	dd_writer(packet, 'res.old')
	axes.plot(wldat_new, res_new, c=colour_config['residual'])
	packet = {
		'wldat': wldat_new,
		'res': res_new
	}
	dd_writer(packet, 'res.new')

def uniqueList(alist):
	"""Based on http://mail.python.org/pipermail/python-list/2000-February/025219.html """
	uniquedict = {}
	for i in alist:
		uniquedict[i] = 0
	nlist = uniquedict.keys()
	nlist.sort()
	return nlist

def assignCompColor(pc, cNo, config):
	cList = [c[0] for c in pc]
	for c in config:
		try:
			if re.match(c['match'], cList[cNo]):
				return c['colour']
		except IndexError:
			# print 'cList'
			# print cList 
			# print "c['match']"
			# print c['match']
			# print 'cNo'
			# print cNo
			pass
	return 'black'
	
def groupStructure(twl, tcom, tsp, pc, groupdelta):
	# print twl
	# print tcom
	# print tsp
	cmp_z = [float(c[4]) for c in pc]
	
	gst = [[gwl, gtcom, gtsp, gwl/(cmp_z[gtcom-1]+1.0)] for gwl, gtcom, gtsp in zip(twl, tcom, tsp)]
	gst.sort(key=lambda g: g[0])
	
	tcoms = uniqueList(tcom)
	
	ingroups = []
	notgroups = []
	
	for com in tcoms:
		newgroup = False
		gst_com = filter(lambda g: g[1] == com, gst)
		if len(gst_com) > 1:
			septotal = 0.0
			group = []
			for i in range(1, len(gst_com)):
				separation = gst_com[i][3]-gst_com[i-1][3]
				if separation + septotal <= groupdelta:
					if not newgroup: newgroup = True
					septotal += separation
					group.append(gst_com[i-1])
					if i == len(gst_com)-1:
						group.append(gst_com[i])
						ingroups.append(group)
				else:
					if newgroup:
						group.append(gst_com[i-1])
						ingroups.append(group)
						group = []
						newgroup = False
						septotal = 0.0
					else:
						notgroups.append(gst_com[i-1])
					if i == len(gst_com)-1:
						notgroups.append(gst_com[i])
		else:
			notgroups.append(gst_com[0])
	
	# print "groups: ", ingroups
	# print "notgroups: ", notgroups
	return ingroups, notgroups

def parseComps(comps): # returns parsed components, still as strings.
	s0p, s1p = "[\s^\n]*", "[\s^\n]+"                # white space
	lin = "([\.\w\* <>\?]+?)"                        # species
	flt = "([0-9]*\.?[0-9]+)"                        # float
	pmf = "([-+]?[0-9]*\.?[0-9]+)"                   # +/- float
	elt = "([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)" # E float
	lbl = "([a-z]{0,2})?([A-Z%]{0,2})?"              # label
	nbl = "([a-z]{0,2})?([A-Z]{0,2}|%|&)?"           # N label
	ing = "(0|([1-9]+[0-9]*))"                       # integer
	ino = "(0|([1-9]+[0-9]*))?"                      # optional integer
	#                 species   N             z             b             e             bturb     temp      rgnflg
	reptrn = "^"+s0p  +lin+s1p  +pmf+nbl+s1p  +flt+lbl+s1p  +pmf+lbl+s1p  +elt+lbl+s1p  +flt+s1p  +elt+s1p  +ino
	
	rt = []
	i = 0
	for c in comps:
		p = re.findall(reptrn, c[C_L])
		if p:
			p0 = p[0]
			# species0, N1,lbl2,small3, z4,lbl5,small6, b7,lbl8,small9, e10,lbl11,small12, bturb13, temp14, rgnflg15, id16, linecopy17
			pci = [p0[0], p0[1], p0[2]+p0[3], p0[2] != '', p0[4], p0[5]+p0[6], p0[5] != '', p0[7], p0[8]+p0[9], p0[8] != '', p0[10], p0[11]+p0[12], p0[11] != '', p0[13], p0[14], p0[16], i, c[C_L]]
			rt.append(pci)
			i += 1
			# print pci
		else:
			if c[C_L].strip()[0] == '!': # commented line
				continue
			else:
				raise fort13FormatError
				# rt.append(None) # something went wrong so throw a spanner
		
	return rt


def plotTicks(rft_old, rft_new, cmpList_old, cmpList_new, axes, settings, config):
	twl_new = tExtract(rft_new, T_WL)
	tcom_new = tExtract(rft_new, T_COM)
	tsp_new = tExtract(rft_new, T_SPEC)
	twl_old = tExtract(rft_old, T_WL)
	tcom_old = tExtract(rft_old, T_COM)
	tsp_old = tExtract(rft_old, T_SPEC)
	
	pcn = parseComps(cmpList_new)
	pco = parseComps(cmpList_old)
	
	# tick positioning
	y0y1_new = [1.0, 2.0]
	ytxt_new = 0.5
	y0y1_old = [3.0, 4.0]
	ytxt_old = 2.4
	
	if settings['tick_type'] == 1:
		drawTicks(twl_new, tcom_new, tsp_new, pcn, axes, config, settings, y0y1_new, ytxt_new, fpart = 'new')
		drawTicks(twl_old, tcom_old, tsp_old, pco, axes, config, settings, y0y1_old, ytxt_old, fpart = 'old')
	elif settings['tick_type'] == 2:
		drawGroupedTicks(twl_new, tcom_new, tsp_new, pcn, axes, config, settings, y0y1_new, ytxt_new, fpart = 'new')
		drawGroupedTicks(twl_old, tcom_old, tsp_old, pco, axes, config, settings, y0y1_old, ytxt_old, fpart = 'old')
	elif settings['tick_type'] == 3:
		drawGroupedTicks(twl_new, tcom_new, tsp_new, pcn, axes, config, settings, y0y1_new, ytxt_new, weighted = True, fpart = 'new')
		drawGroupedTicks(twl_old, tcom_old, tsp_old, pco, axes, config, settings, y0y1_old, ytxt_old, weighted = True, fpart = 'old')
	
	taxis = axes.axis() # xmin, xmax, ymin, ymax
	axes.axis([taxis[0], taxis[1], 0.0, 4.5])

def tickText(compid, pc, settings):
	txt = ''
	if settings['tick_label'] == 1:
		if pc[compid][5] != '': # component has a label
			species_pc = [c for c in pc if c[0] == pc[compid][0]]
			if len(species_pc) == 1: # only one component of this species
				txt = pc[compid][5]
			else: # multiple components of this species
				# find all components with this exact label, if there are more than 1, set txt to compid with label as subscript (fixed redshift)
				if len([c for c in species_pc if c[5] == pc[compid][5]]) > 1:
					txt = "%i$_\mathrm{\mathsf{%s}}$" % (compid+1, pc[compid][5])
				else:
					if len(pc[compid][5]) == 2:
						txt = "$_\mathrm{\mathsf{%s}}$%s" % (pc[compid][5][0], pc[compid][5][1])
					else:
						txt = pc[compid][5]
					# Alternate:
					# #check if all labels are of length 2 and have common first or second character. If so, don't include it.
					# albl = [c for c in species_pc if len(c[5]) == 2]
					# if len(albl) > 1:
					# 	clbl = [c for c in albl if c[5][0].lower() == albl[0][5][0].lower()]
					# 	dlbl = [c for c in albl if c[5][1].lower() == albl[0][5][1].lower()]
					# 	if len(albl) == len(clbl):
					# 		txt = pc[compid][5][1]
					# 	elif len(albl) == len(dlbl):
					# 		txt = pc[compid][5][0]
					# 	else:
					# 		txt = pc[compid][5]
					# else:
					# 	txt = pc[compid][5]
		else:
			txt = str(compid+1)
	if settings['tick_label'] == 2:
		txt = str(compid+1)
	return txt

def drawTicks(twl, tcom, tsp, pc, axes, config, settings, y0y1, ytxt, fline = None, fpart = ''):
	dd_fpart = fpart+'.ticks'
	packet_color = []
	packet_tx = []
	packet_label = []
	packet_text = []
	# packet_zspecies = []
	
	for i, wl in enumerate(twl):
		tcol = assignCompColor(pc, tcom[i]-1, config)
		# packet_zspecies.append(pc[tcom[i]-1][0])
		packet_color.append(tcol)
		if fline:
			wv_obs = float(fline[2]['wv'])*(fline[1]+1.0)
			tx_i = (wl-wv_obs)/wv_obs*C/1000.0
		else:
			tx_i = wl
		packet_tx.append(tx_i)
		axes.plot([tx_i, tx_i], y0y1, color=tcol)
		tx_lbl = pc[tcom[i]-1][17].strip('\n')
		packet_label.append('"'+tx_lbl+'"')
		tx_txt = tickText(tcom[i]-1, pc, settings)
		packet_text.append(tx_txt)
		if ytxt:
			axes.text(tx_i, ytxt, tx_txt, horizontalalignment='center', size = 'smaller', color=tcol, picker=2, label=tx_lbl)
		
	packet = {
		'tx': packet_tx,
		'zz_copy': packet_label,
		'text': packet_text,
		'color': packet_color
		# 'zspecies': packet_zspecies
	}
	dd_writer(packet, dd_fpart)

def drawGroupedTicks(twl, tcom, tsp, pc, axes, config, settings, y0y1, ytxt, fline = None, weighted = False, fpart = ''):
	groupdelta = settings['group_delta']
	ingroups, notgroups = groupStructure(twl, tcom, tsp, pc, groupdelta)
	if notgroups:
		ng_wl = [ng[0] for ng in notgroups]
		ng_tcom = [ng[1] for ng in notgroups]
		ng_tsp = [ng[2] for ng in notgroups]
		drawTicks(ng_wl, ng_tcom, ng_tsp, pc, axes, config, settings, y0y1, ytxt, fline = fline, fpart = fpart+'.notgrp')
	if ingroups:
			drawGroups(ingroups, pc, axes, config, settings, y0y1, ytxt, fline = fline, weighted = weighted, fpart = fpart+'.ingrp')

def drawGroups(ingroups, pc, axes, config, settings, y0y1, ytxt, fline = None, weighted = False, fpart = ''):
	dd_fpart = fpart+'.ticks'
	
	if not weighted:
		y0y1 = [y0y1[0], y0y1[0], y0y1[1], y0y1[1]]
	
	packet_color = []
	packet_txa = []
	packet_txb = []
	packet_label = []
	packet_text = []
	# packet_zspecies = []
	
	for group in ingroups:
		g_wl = [g[0] for g in group]
		g_tcom = [g[1] for g in group]
		g_tsp = [g[2] for g in group]
		tcol = assignCompColor(pc, g_tcom[0]-1, config)
		packet_color.append(tcol)
		# packet_zspecies.append(pc[tcom[i]-1][0])
		if fline: wv_obs = float(fline[2]['wv'])*(fline[1]+1.0)
		if weighted:
			#                         species               approx_rest_wave
			g_f = [float(find_line(   pc[g_tcom_i-1][0], g_wl_i/(float(pc[g_tcom_i-1][4])+1.0)   )['f']) for g_tcom_i, g_wl_i in zip(g_tcom, g_wl)]
			g_weighted_wl = sum([g_f_i * g_wl_i for g_f_i, g_wl_i in zip(g_f, g_wl)]) / sum(g_f)
			if fline:
				# calculate velocity of weighted component
				g_weighted_x = (g_weighted_wl-wv_obs)/wv_obs*C/1000.0
			else:
				g_weighted_x = g_weighted_wl
			axes.plot([g_weighted_x, g_weighted_x], y0y1, color=tcol)
			
			tx_lbl = pc[g_tcom[0]-1][17].strip('\n')
			packet_label.append('"'+tx_lbl+'"')
			tx_txt = tickText(g_tcom[0]-1, pc, settings)
			packet_text.append(tx_txt)
			
			packet_txa.append(g_weighted_x)
			packet_txb.append(g_weighted_x)
			
			if ytxt:
				axes.text(g_weighted_x, ytxt, tx_txt, horizontalalignment='center', size = 'smaller', color=tcol, picker=2, label=tx_lbl)
		else:
			gwl_min, gwl_max = min(g_wl), max(g_wl)
			if fline:
				gx_min, gx_max = (gwl_min-wv_obs)/wv_obs*C/1000.0, (gwl_max-wv_obs)/wv_obs*C/1000.0
			else:
				gx_min, gx_max = gwl_min, gwl_max
			axes.fill([gx_min, gx_max, gx_max, gx_min], y0y1, color=tcol, alpha=0.6)
			
			tx_lbl = pc[g_tcom[0]-1][17].strip('\n')
			packet_label.append('"'+tx_lbl+'"')
			tx_txt = tickText(g_tcom[0]-1, pc, settings)
			packet_text.append(tx_txt)
			
			packet_txa.append(gx_min)
			packet_txb.append(gx_max)
			
			if ytxt:
				axes.text((gx_min+gx_max)/2.0, ytxt, tx_txt, horizontalalignment='center', size = 'smaller', color=tcol, picker=2, label=tx_lbl)
	
	packet = {
		'txa': packet_ta,
		'txb': packet_tb,
		'zz_copy': packet_label,
		'text': packet_text,
		'color': packet_color
		# 'zspecies': packet_zspecies
	}
	dd_writer(packet, dd_fpart)

def hideXLabels(axes):
	""" http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg06932.html """
	for label in axes.get_xticklabels():
		label.set_visible(False)

def onlySpaces(astring): # not very brilliant but works - ... seriously, what does this do again? - appears to replace whitespace with spaces
	line = astring.split()
	r = ''
	for l in line:	r += l + ' '
	return r

########

def deltacolor(i, fit, labelkey, valkey = None):
	iiter = len(fit['fitprogress'])-1
	fl = fit['fitprogress'][-1]['lines'][i]
	if fl[labelkey].upper() == fl[labelkey] and fl[labelkey] != '':
		return '${BOLD}${BLACK}'
	else:
		if not valkey or iiter == 0:
			return '${WHITE}'
		else:
			if fit['fitprogress'][-1]['try'] == fit['fitprogress'][-2]['try']:
				delta = float(fit['fitprogress'][-1]['lines'][i][valkey])-float(fit['fitprogress'][-2]['lines'][i][valkey])
				if delta < 0.0:
					return '${CYAN}'
				elif delta > 0.0:
					return '${YELLOW}'
				else:
					return '${WHITE}'
			else:
				return '${WHITE}'
	

def iterprint(fit, enditer):
	iiter = len(fit['fitprogress'])-1
	
	if iiter == 0:
		print term.render('Fitting file ${BOLD}%s${NORMAL} with ${BOLD}%s${NORMAL} region%s and ${BOLD}%s${NORMAL} ion%s...' 
			% (fit['f13path'], fit['nregions'], ('s','')[int(fit['nregions']) == 1], fit['nions'], ('s','')[int(fit['nions']) == 1]),)
		print
		print term.render('${BG_BLACK}   ion          N            z             b         x4           bturb   temp             ${NORMAL}')
		print
		chi2delta = ''
	else:
		chi2delta = ' (%.4f)' % (float(fit['fitprogress'][-1]['chi2'])-float(fit['fitprogress'][-2]['chi2']))
	
	print term.render(' iter ${BOLD}%s${NORMAL} | try %s${BOLD}%s${NORMAL} | total ${BOLD}%i${NORMAL}                chi^2/v = ${BOLD}%s${NORMAL} = ${BOLD}%s${NORMAL}%s / ${BOLD}%s${NORMAL}'
		% (fit['fitprogress'][-1]['iter'], ('', '${YELLOW}')[int(fit['fitprogress'][-1]['try']) > 1], fit['fitprogress'][-1]['try'], iiter,
		fit['fitprogress'][-1]['chi2v'], fit['fitprogress'][-1]['chi2'], chi2delta, fit['fitprogress'][-1]['v']))
	
	for i, l in enumerate(fit['fitprogress'][-1]['lines']):
		print term.render('${BG_BLACK}   %s %s%s${NORMAL}${BG_BLACK}%s%s${NORMAL}${BG_BLACK} %s%s${NORMAL}${BG_BLACK}%s%s${NORMAL}${BG_BLACK} %s%s${NORMAL}${BG_BLACK}%s%s${NORMAL}${BG_BLACK} %s%s${NORMAL}${BG_BLACK}%s%s${NORMAL}${BG_BLACK} %s%s${NORMAL}${BG_BLACK} %s%s${NORMAL}${BG_BLACK} %s ${BOLD}${BLACK}%s${NORMAL}'
			% (l['ion'].ljust(8), deltacolor(i, fit, 'Nlbl', 'N'), l['N'].rjust(8), deltacolor(i, fit, 'Nlbl'), l['Nlbl'].ljust(2),
			deltacolor(i, fit, 'zlbl', 'z'), l['z'].rjust(11), deltacolor(i, fit, 'zlbl'), l['zlbl'].ljust(2),
			deltacolor(i, fit, 'blbl', 'b'), l['b'].rjust(9), deltacolor(i, fit, 'blbl'), l['blbl'].ljust(2),
			deltacolor(i, fit, 'x4lbl', 'x4'), l['x4'].rjust(9), deltacolor(i, fit, 'x4lbl'), l['x4lbl'].ljust(2),
			deltacolor(i, fit, 'blbl', 'bturb'), l['bturb'].rjust(9),
			deltacolor(i, fit, 'blbl', 'temp'), l['temp'].rjust(10),
			l['region'].rjust(2),
			l['cmid']))
	print
	
	if fit['fitprogress'][-1]['msg'] != '':
		for ml in fit['fitprogress'][-1]['msg'].split('\n'):
			if ml.strip() != '': print term.render('${BG_BLACK}${MAGENTA}%s${NORMAL}' % (ml.ljust(92)))
		print
	
	if enditer:
		print term.render('${BG_GREEN}${BLACK}'+'-- Fitting complete! --'.center(92)+'${NORMAL}')
		print
		print ' Parameter errors:'
		for il in fit['finalstats']['err'].split('\n'):
			if il.strip() != '': print term.render('${BG_BLACK}${WHITE}%s${NORMAL}') % (il.ljust(92))
		print
		print ' Fit statistics:'
		for il in fit['finalstats']['fitstats'].split('\n'):
			print term.render('${BG_BLACK}${WHITE}%s${NORMAL}') % (il.ljust(92))
		print
		print ' Region statistics:'
		for il in fit['finalstats']['regstats'].split('\n'):
			print term.render('${BG_BLACK}${WHITE}%s${NORMAL}') % (il.ljust(92))
		print

def livefit(vpfit_path, fort13, nvar, bvar, zvar, x4var, csvar, itercb):
	try:
		vpf = VPFIT(vpfit_path)
		vpf.fit(f13path = fort13,
			n = nvar, b = bvar, z = zvar, x4 = x4var, cs = csvar,
			itercallback = itercb)
	except KeyboardInterrupt:
		print term.render('${BG_BLACK}${RED}'+'-- Terminating --'.center(92)+'${NORMAL}')
		del vpf
	
	return