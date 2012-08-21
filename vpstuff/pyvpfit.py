from vpstuff.pexpect import *
import sys, time

RE_UFLOAT = lambda name: '(?P<%s>[0-9]*\.?[0-9]+)' % name
RE_SFLOAT = lambda name: '(?P<%s>[-+]?[0-9]*\.?[0-9]+)' % name
RE_EFLOAT = lambda name: '(?P<%s>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)' % name
RE_UINT   = lambda name: '(?P<%s>0|([1-9]+[0-9]*))' % name
RE_UINTOP = lambda name: '(?P<%s>(0|([1-9]+[0-9]*))?)' % name
RE_LCCHAR = lambda name: '(?P<%s>[a-z])' % name
RE_SPEC   = lambda name: '(?P<%s>[\.\w\* <>\?]+?)' % name
RE_NLBL   = lambda name: '(?P<'+name+'>([a-z]{0,2})?([A-Z]{0,2}|%|&)?)'
RE_LBL    = lambda name: '(?P<'+name+'>([a-z]{0,2})?([A-Z%]{0,2})?)'
RE_ALL    = lambda name: '(?P<%s>.*?)' % name
RE_CMID   = lambda name: '(?P<%s>[! 0-9]+)' % name
RE_W0     = '[ \t\f\v]*'
RE_W1     = '[ \t\f\v]+'
RE_NL     = '\r\n' # VPFIT prints these for newlines

# EXPECTS:
# novars 4
# nopchan 2 ?
# pcvals
# 'linear' n vals

class VPFIT(object):
	ITER_TIMEOUT = 1*60*60 # seconds
	
	def __init__(self, vpf_path):
		self.vpf_version = None
		self.setup = {}
		self.fit_data = {'fitprogress': [], 'finalstats': {}}
		
		self.vpf_path = vpf_path
		self.vpf = spawn(vpf_path)
		# self.vpf.logfile = sys.stdout #DEBUG
		
		self.vpf.expect('VPFIT %s' % (RE_UFLOAT('vpf_version')), timeout=5)
		self.vpf_version = self.vpf.match.groupdict()['vpf_version']
		if float(self.vpf_version) != 10.0:
			raise VersionError("Unexpected VPFIT version; expected 10.0")
		
		try:
			self.vpf.expect(RE_W0+'Total column density after:'+RE_W0+RE_LCCHAR('lastchtied'), timeout = 1)
			self.setup.update(self.vpf.match.groupdict())
		except TIMEOUT:
			self.setup['lastchtied'] = None
		
		self.vpf.expect(RE_W0+'Using data from : (?P<atom>.+?)')
		self.setup.update(self.vpf.match.groupdict())
		
	def fit(self, f13path = 'fort.13', itercallback = None, **kwargs):
		self.fit_data['f13path'] = f13path
		self.vpf.expect(RE_W0+'option \(key\) \(key\)\.\.\.'+RE_W0+RE_NL+RE_W0+'>')
		self.vpf.sendline('f')
		# print kwargs
		for key in kwargs:
			if kwargs[key] and str(kwargs[key]).strip() != '':
				self.vpf.expect(RE_W0+'setup: \? print values, <CR> OK, n,z,b,x4,cs,sf,il,w,me,p,d,v to change'+RE_W0+RE_NL+RE_W0+'>')
				self.vpf.sendline(key)
				self.vpf.expect(RE_W0+'New value.*?'+RE_W0+RE_NL+RE_W0+'>')
				self.vpf.sendline(str(kwargs[key]))
			# print key
		self.vpf.sendline()
		self.vpf.expect(RE_W0+'Column density \(n\), logN \(l\) or emission \(e\), scalefactor \[l, 1\.0\]'+RE_W0+RE_NL+RE_W0+'>')
		self.vpf.sendline()
		self.vpf.expect(RE_W0+'Parameter input file, # entries\? \[fort\.13,1\]'+RE_W0+RE_NL+RE_W0+'>')
		self.vpf.sendline(self.fit_data['f13path'])
		self.vpf.expect(RE_W0+RE_UINT('nregions')+RE_W0+'regions fitted'+RE_W0)
		self.fit_data.update(self.vpf.match.groupdict())
		self.vpf.expect(RE_W0+'no\. of ions for fitting is'+RE_W0+RE_UINT('nions')+RE_W0)
		self.fit_data.update(self.vpf.match.groupdict())
		# nions = int(self.fit_data['nions'])
		
		reptrn = RE_W0+RE_SPEC('ion') \
			+RE_W1+RE_SFLOAT('N')+RE_NLBL('Nlbl')+RE_W1 \
			+RE_UFLOAT('z')+RE_LBL('zlbl')+RE_W1 \
			+RE_SFLOAT('b')+RE_LBL('blbl')+RE_W1 \
			+RE_EFLOAT('x4')+RE_LBL('x4lbl')+RE_W1 \
			+RE_EFLOAT('bturb')+RE_W1 \
			+RE_EFLOAT('temp')+RE_W1 \
			+RE_UINTOP('region') \
			+RE_CMID('cmid')+RE_NL
		
		while 1:
			self.vpf.expect(RE_W0+'iteration'+RE_W0+':'+RE_W0+RE_UINT('iter')+RE_W0+'\('+RE_W0+RE_UINT('try')+RE_W0+'\)'+RE_W0+RE_NL
				+RE_W0+'chi-squared'+RE_W0+':'+RE_W0+RE_UFLOAT('chi2v')+RE_W0+'\('+RE_W0+RE_UFLOAT('chi2')+RE_W0+','+RE_W0+RE_UFLOAT('v')+RE_W0+'\)'+RE_W0+RE_NL,
				timeout = self.ITER_TIMEOUT)
			self.fit_data['fitprogress'].append(self.vpf.match.groupdict())
			self.vpf.expect(RE_NL)
			
			
			self.fit_data['fitprogress'][-1]['lines'] = []
			self.fit_data['fitprogress'][-1]['msg'] = ''
			
			ended = False
			earlystop = False
			while 1:
				iterdone = False
				e = self.vpf.expect([reptrn, RE_W0+'(?=iteration'+RE_W0+':)', '('+RE_W0+'iterations stopped early'+RE_W0+RE_NL+')', RE_W0+'Parameter errors:'+RE_W0+RE_NL,RE_ALL('msg')+RE_NL, TIMEOUT], timeout = (5, self.ITER_TIMEOUT)[earlystop])
				# print e
				if e == 0: self.fit_data['fitprogress'][-1]['lines'].append(self.vpf.match.groupdict())
				elif e == 1: break
				elif e == 2:
					earlystop = True
				elif e == 3:
					self.fit_data['finalstats'] = {'err': '', 'fitstats': '', 'regstats': ''}
					# print 'at err'
					while 1:
						e = self.vpf.expect([RE_W0+'statistics for whole fit:'+RE_W0+RE_NL, RE_ALL('err')+RE_NL], timeout = 5)
						if e == 1: self.fit_data['finalstats']['err'] += self.vpf.match.groupdict()['err'] + '\n'
						elif e == 0: break
					while 1:
						e = self.vpf.expect([RE_W0+'Statistics for each region :'+RE_W0+RE_NL, RE_ALL('fitstats')+RE_NL], timeout = 5)
						if e == 1: self.fit_data['finalstats']['fitstats'] += self.vpf.match.groupdict()['fitstats'] + '\n'
						elif e == 0: break
					while 1:
						e = self.vpf.expect([RE_W0+'Plot\? y,n, c=change device, t=change ticks \[y\]'+RE_W0+RE_NL+'>', RE_ALL('regstats')+RE_NL], timeout = 5)
						if e == 1: self.fit_data['finalstats']['regstats'] += self.vpf.match.groupdict()['regstats'] + '\n'
						elif e == 0: break
					# print 'finished err'
					self.vpf.sendline('n')
					self.vpf.expect(RE_W0+'Fit more lines\? \[n\]'+RE_W0+RE_NL+'>')
					self.vpf.sendline('n')
					ended = True
					iterdone = True
				elif e == 4:
					pass
					# print 'at 3'
					# try: # catch any messages
					while 1:
						# print 'waiting'
						e = self.vpf.expect([RE_W0+'(?=(iteration'+RE_W0+':))', RE_W0+'(?=(Parameter errors)|(iterations stopped early))', RE_ALL('msg')+RE_NL, TIMEOUT], timeout = 5)
						# print 'emsg', e
						if e == 2: 
							if self.vpf.match.groupdict()['msg'].strip() != 'ion          N            z             b         bturb   temp' and self.vpf.match.groupdict()['msg'].strip() != '':
								self.fit_data['fitprogress'][-1]['msg'] += self.vpf.match.groupdict()['msg']+'\n'
						elif e == 0 or e == 3:
							iterdone = True
							break
						elif e == 1:
							iterdone = False
							break
				elif e == 5:
					iterdone = True
				if iterdone:
					if itercallback:
							itercallback(self.fit_data, enditer = ended)
					break
				if earlystop:
					if itercallback:
							itercallback(self.fit_data, enditer = ended)
			# print 'out2'

			if ended: break
			# print self.fit_data['fitprogress'][-1]['msg']
		# print self.fit_data

class VPFITError(Exception):
	pass
	
class VersionError(VPFITError):
	def __init__(self, msg):
		self.msg = msg