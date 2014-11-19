#!/usr/bin/python2
import frag
import fio
import os.path
import argparse


# TODO collect part

supportedMethod = ["SFM", "MTA", "XO1", "CFM", "XF"]
supportedLevel = range(1, 5)
flagLayer = 'HMLN'


class Interpreter(object):
	method = ""
	level = 3
	xoutput = True
	verbose = False
	forceMatrix = []
	fragmentor = None
	def __init__(self, arg):
		filename = arg.COORD
		
		self.method = arg.method	
		self.level = arg.level
		self.xoutput = arg.xoutput
		self.verbose = arg.verbose
		self.layer = arg.layer
		self.collect = arg.collect
		if filename[-3:] in fio.typeConverter:
			if self.verbose: print "  readin molecule from file: " + filename
			self.molecule = fio.typeConverter[filename[-3:]](filename)
			if self.verbose: print self.molecule
		#if self.verbose:
		#	print self.molecule
		if self.method == "XO1":
			if self.verbose: print "  readin force matrix from file: " + filename[:-3] + "log"
			self.forceMatrix = fio.forceMatrixFromGaussOutput(filename[:-3] + "log")
			if self.verbose: print self.forceMatrix

	def XO1(self):
		self.fragmentor = frag.iFrag(self.molecule, 3.5 + (self.level - 3) * 0, self.level * 8 / 5)
		self.fragmentor.loadForceMatrix(self.forceMatrix)
		self.fragmentor.Run()

	def MTA(self):
		self.fragmentor = frag.MTA(self.molecule, self.level * 1., self.level * 5)
		self.fragmentor.Run()


	def XF(self):
		self.fragmentor = frag.FragForOpt(self.molecule, 3.5 + (self.level - 3) * 0 , self.level *  8 / 5.)
		self.fragmentor.Run()


	def CFM(self):
		self.fragmentor = frag.CFM(self.molecule, self.level - 1)
		self.fragmentor.Run()

	def SFM(self):
		self.fragmentor = frag.SFM(self.molecule, self.level)
		self.fragmentor.Run()

	def Frag(self):
		if self.method == "XO1":
			self.XO1()
		elif self.method == "MTA":
			self.MTA()
		elif self.method == "CFM":
			self.CFM()
		elif self.method == "SFM":
			self.SFM()
		elif self.method == "XF":
			self.XF()
		if self.verbose:
			print "Number of %s fragments: %d" % (self.method, len(self.fragmentor.frags))

	def output(self):
		if self.xoutput:
			eggs =  self.fragmentor.xout()
			print "[eggs]"
			print eggs
			eggs = eggs.split('\n')
			print "[levels]"
			for i in eggs[:-1]:
				print i[:i.index(':')] + ': ' + self.layer



def main(argv = None):
	parser = argparse.ArgumentParser(
			description = "Fragmentation program for SFM, MTA, CFM, XO1.")
	parser.add_argument("COORD", 
			help = "molecule file, including coordinates and linkage list, gau/gjf type only")
	parser.add_argument("-f", 
			dest = "method", 
			choices = supportedMethod, 
			required = True, 
			help = "Supporded fragmentation method, XO1 needs force(-error) list in gaussian-log file")
	parser.add_argument("-n", 
			dest = "level", 
			type = int, 
			default = 3, 
			help = "precision level of fragmentation, the higher it is, the bigger each fragment is, default 2")
	parser.add_argument("-x", "--xo", 
			dest = "xoutput", 
			action = "store_const", 
			const = True, 
			default = False, 
			help = "use xo-schemefile output, default=False, giving off gjf-files")
	parser.add_argument("-v", 
			"--verbose", 
			dest = "verbose",
			action = "store_const", 
			const = True, 
			default = False, 
			help = "run in vervose mode")
	parser.add_argument("-l", 
			dest = "layer", 
			choices = flagLayer, 
			default = 'H', 
			help = "layer of this fragmentation, default is H")
	parser.add_argument("-c", 
			"--collect", 
			dest = "collect",
			action = "store_const", 
			const = True, 
			default = False, 
			help = "collect calculation results")
	arg = parser.parse_args()
	ip = Interpreter(arg)
	if arg.method in supportedMethod:
		ip.Frag()
		if arg.xoutput: ip.output()
		if arg.collect: 
			pass
	return 0

if __name__  == "__main__":
	reselt = main()
