#!/usr/bin/env python

import optparse
import os,sys
#from optparse import OptionParser
import glob
import subprocess
import linecache
import struct
import shutil

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -f <stack>")
	parser.add_option("-f",dest="stack",type="string",metavar="FILE",
		help="IMAGIC particle stack")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
		help="debug")
	options,args = parser.parse_args()

	if len(args) > 0:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) <= 0:
		parser.print_help()
		sys.exit()
	params={}
	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params

def checkConflicts(params):
        if not params['stack']:
                print "\nWarning: no stack specified\n"

def getEMANPath():
        ### get the imagicroot directory        
        emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()

        if emanpath:
                emanpath = emanpath.replace("EMAN2DIR=","")
        if os.path.exists(emanpath):
                return emanpath
        print "EMAN2 was not found, make sure it is in your path"
        sys.exit()

def run(params):

	stack = params['stack']
	
        # get box size
        im=EMData.read_images(stack,[0])
        nx = im[0].get_xsize()
        del im
	nimg = EMUtil.get_image_count(stack)

        img = EMData(nx,nx,nimg)
        img.write_image(stack[:-4]+'.mrc')

	i = 1

        while i < nimg:
                d = EMData()
                d.read_image(stack, i)
                region = Region(0, 0, i, nx, nx, 1)
                d.write_image(stack[:-4]+".mrc",0,EMUtil.get_image_ext_type("mrc"), False, region, EMUtil.EMDataType.EM_FLOAT, True)
        	i = i + 1

#=========================
#=========================
if __name__ == "__main__":
        params=setupParserOptions()

        getEMANPath()
        from EMAN2 import *
        from sparx  import *

        checkConflicts(params)
        run(params)

