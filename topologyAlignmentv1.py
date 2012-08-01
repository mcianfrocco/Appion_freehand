#!/usr/bin/env python
#
import os,sys,re
import time
import math
import shutil
import glob
import cPickle
import tarfile
import subprocess
#appion
from appionlib import appionScript
from appionlib import apDisplay
from appionlib import apFile
from appionlib import apStack
from appionlib import apParam
from appionlib import apEMAN
from appionlib import appiondata
from appionlib import apProject
from appionlib.apSpider import operations
from appionlib import apIMAGIC
from appionlib import apImagicFile
from appionlib.apImagic import imagicFilters
from appionlib.apImagic import imagicAlignment

#=====================
#=====================
class TopologyRepScript(appionScript.AppionScript):

	#=====================
	def setupParserOptions(self):
		self.parser.set_usage("Usage: %prog --stack=ID --start=# --end=# [options]")
		self.parser.add_option("-N", "--num-part", dest="numpart", type="int",
			help="Number of particles to use", metavar="#")
		self.parser.add_option("-s", "--stack", dest="stackid", type="int",
			help="Stack database id", metavar="ID#")

		self.parser.add_option("--msaproc", dest="msaproc", type="int", default=1,
			help="Number of processor to use for CAN", metavar="#")

		self.parser.add_option("--lowpass", "--lp", dest="lowpass", type="int",
			help="Low pass filter radius (in Angstroms)", metavar="#")
		self.parser.add_option("--highpass", "--hp", dest="highpass", type="int",
			help="High pass filter radius (in Angstroms)", metavar="#")
		self.parser.add_option("--bin", dest="bin", type="int", default=1,
			help="Bin images by factor", metavar="#")
		self.parser.add_option("-i", "--iter", dest="iter", type="int", default=20,
			help="Number of iterations", metavar="#")
		self.parser.add_option("--start", dest="start", type="int",
			help="number of classes to create in first iteration")
		self.parser.add_option("--end", dest="end", type="int",
			help="number of classes to create in last iteration")
		self.parser.add_option("--mask", dest="mask", type='int', metavar="#",
			help="radius of circular mask to apply in pixels (default=(boxsize/2)-2)")
		self.parser.add_option("--itermult", dest="itermult", type="float", metavar="FLOAT", default=10.0,
			help="multiplier for determining number of times data will be presented to the network. Number of particles in your stack will by multiplied by this value to determine # of iterations")
		self.parser.add_option("--learn", dest="learn", type="float", metavar="FLOAT", default=0.01,
			help="direct learning rate - fraction that closest unit image will be moved toward presented data, 0.01 suggested for cryo, higher for neg stain")
		self.parser.add_option("--ilearn", dest="ilearn", type="float", metavar="FLOAT", default=0.0005,
			help="indirect learning rate - fraction that connection unit images will be moved should be lower than direct rate")
		self.parser.add_option("--age", dest="maxage", type="int", metavar="INT", default=25,
			help="number of iterations an edge connecting two units can be unused before it's discarded")

		### IMAGIC MSA options
		self.parser.add_option("--msaiter", dest="msaiter", type="int", default=50,
			help="number of MSA iterations")
		self.parser.add_option("--numeigen", dest="numeigen", type="int", default=20,
			help="total number of eigen images to calculate")
		self.parser.add_option("--overcorrection", dest="overcorrection", type="float", default=0.8,
			help="overcorrection facter (0-1)")
		self.parser.add_option("--activeeigen", dest="activeeigen", type="int", default=10,
			help="number of active eigen images to use for classification")

		### true/false
		self.parser.add_option("--keep-all", dest="keepall", default=False,
			action="store_true", help="Keep all intermediate node images")
		self.parser.add_option("--premask", dest="premask", default=False,
			action="store_true", help="Mask raw particles before processing")
		self.parser.add_option("--no-mask", dest="nomask", default=False,
			action="store_true", help="Do not apply a mask to the class averages")
		self.parser.add_option("--no-center", dest="nocenter", default=False,
			action="store_true", help="Do not center particles after each iteration")
		self.parser.add_option("--classiter", dest="classiter", default=False,
			action="store_true", help="Perform iterative averaging of class averages")
		self.parser.add_option("--uploadonly", dest="uploadonly", default=False,
			action="store_true", help="Just upload results of completed run")
		self.parser.add_option("--invert", dest="invert", default=False,
			action="store_true", help="Invert before alignment")

		### choices
		self.mramethods = ("eman","imagic")
		self.parser.add_option("--mramethod", dest="mramethod",
			help="Method for multi-reference alignment", metavar="PACKAGE",
			type="choice", choices=self.mramethods, default="eman")
		self.msamethods = ("can","imagic")
		self.parser.add_option("--msamethod", dest="msamethod",
			help="Method for MSA", metavar="PACKAGE",
			type="choice", choices=self.msamethods, default="can")
		self.cluster = ("barcelona","himem")
		self.parser.add_option("--cluster", dest="cluster", 
			help="If using cluster, specify queue (defaults to barcelona)", metavar="QUEUE",
			type="choice", choices=self.cluster, default="barcelona")

	def checkConflicts(self):
		if self.params['stackid'] is None:
			apDisplay.printError("stack id was not defined")
		if self.params['start'] is None:
			apDisplay.printError("a number of starting classes was not provided")
		if self.params['end'] is None:
			apDisplay.printError("a number of ending classes was not provided")
		if self.params['runname'] is None:
			apDisplay.printError("run name was not defined")
		stackdata = apStack.getOnlyStackData(self.params['stackid'], msg=False)
		stackfile = os.path.join(stackdata['path']['path'], stackdata['name'])
		if self.params['numpart'] > apFile.numImagesInStack(stackfile):
			apDisplay.printError("trying to use more particles "+str(self.params['numpart'])
				+" than available "+str(apFile.numImagesInStack(stackfile)))

		self.boxsize = apStack.getStackBoxsize(self.params['stackid'])
		self.workingboxsize = math.floor(self.boxsize/self.params['bin'])
		if self.params['numpart'] is None:
			self.params['numpart'] = apFile.numImagesInStack(stackfile)
		if not self.params['mask']:
			self.params['mask'] = (self.boxsize/2)-2
		self.workingmask = math.floor(self.params['mask']/self.params['bin'])
		if self.params['mramethod'] == 'imagic':
			self.imagicroot = apIMAGIC.checkImagicExecutablePath()
			self.imagicversion = apIMAGIC.getImagicVersion(self.imagicroot)
		## check if running on cluster:
		if apParam.getExecPath("qsub",die=False) is None:
			self.params['cluster']=None


	#=====================
	def insertTopolRepJob(self):
		print 'inside insertTopolRepJob'
		topoljobq = appiondata.ApTopolRepJobData()
		print topoljobq
		topoljobq['runname'] = self.params['runname']
		print "Runname = %s" %(topoljobq['runname'])
		topoljobq['path'] = appiondata.ApPathData(path=os.path.abspath(self.params['rundir']))
		print self.params['rundir']
		print os.path.abspath(self.params['rundir'])
		print topoljobq['path']
		topoljobdatas = topoljobq.query(results=1)
		print "topoljobdatas = %s" %(topoljobdatas)
		if topoljobdatas:
			alignrunq = appiondata.ApAlignRunData()
			alignrunq['runname'] = self.params['runname']
			alignrunq['path'] = appiondata.ApPathData(path=os.path.abspath(self.params['rundir']))
			alignrundata = alignrunq.query(results=1)
			if topoljobdatas[0]['finished'] is True or alignrundata:
				apDisplay.printError("This run name already exists as finished in the database, please change the runname")
		topoljobq['REF|projectdata|projects|project'] = apProject.getProjectIdFromStackId(self.params['stackid'])
		topoljobq['timestamp'] = self.timestamp
		topoljobq['finished'] = False
		topoljobq['hidden'] = False
		if self.params['commit'] is True:
			topoljobq.insert()
		self.params['topoljob'] = topoljobq
		return


	#=====================
	def start(self):
		self.insertTopolRepJob()
		self.stack = {}


#=====================
if __name__ == "__main__":
	topRep = TopologyRepScript()
	topRep.start()
	topRep.close()



