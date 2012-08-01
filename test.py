#!/usr/bin/env python

import re
import subprocess
import os
import sys
def getOPENMPIPath():
        ### get the openmpi directory        
        openpath = subprocess.Popen("env | grep MPIHOME", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
	print openpath
        test = openpath.find('imagic')
	print test
	if test >= 0:
                print "only imagic openmpi is detected, load openmpi module"
                sys.exit()

	if test is None:

	        if openpath:
        	        openpath = openpath.replace("MPIHOME=","")
	        if os.path.exists(openpath):
        	        return openpath
	        print "openmpi is not loaded, make sure it is in your path"
        	sys.exit()

def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

getOPENMPIPath()


