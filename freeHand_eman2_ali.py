#!/usr/bin/env python

import optparse
from sys import *
import os,sys,re
from optparse import OptionParser
import glob
import subprocess
from os import system
import linecache

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog -u <untilted stack> -m <model> -p <parameter file>")
        parser.add_option("-u",dest="untilted",type="string",metavar="FILE",
                help="untilted stack (white particles in IMAGIC format)")
        parser.add_option("-m",dest="model",type="string",metavar="FILE",
                help="3D model(s) for alignment (Single SPIDER volume or multi-volume HDF file)")
        parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Parameter file with refinement info (free_param.par)")
        parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 2:
                parser.print_help()
		sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def getEMANPath():        
        ### get the imagicroot directory        
        emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        

        if emanpath:                
                emanpath = emanpath.replace("EMAN2DIR=","")                
        if os.path.exists(emanpath):                        
                return emanpath        
        print "EMAN2 was not found, make sure eman2/2.05 is in your path"        
        sys.exit()

def getOPENMPIPath():
        ### get the openmpi directory        
        openpath = subprocess.Popen("env | grep MPIHOME", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        test = openpath.find('imagic')
	if test >= 0:
                print "OPENMPI is not loaded, make sure it is in your path"
                sys.exit()

	if test is None:

	        if openpath:
        	        openpath = openpath.replace("MPIHOME=","")
	        if os.path.exists(openpath):
        	        return openpath
	        print "OPENMPI is not loaded, make sure it is in your path"
        	sys.exit()


def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

def main2(params):
	debug = params['debug']
	param = params['param']
	untilt = params['untilted']
	model = params['model']

	#Get current working directory
	script = sys.argv[0]
	cwd = '%s/lib' %(script[:-22])

	#Get parameter info: angular step
	p = open(param,'r')
	a = 'angular' 
	angl = grep(a,p)
	aL = angl.split()
	ang = aL[2]
	
	#Shift
	s = 'shift'
	p = open(param,'r')
	shi = grep(s,p)
	sh = shi.split()
	sx = sh[2]

	#Pixel size
	p13 = open(param,'r')
	pixel = 'pix'
	pixe = grep(pixel,p13)
	pi = pixe.split()
	pix = pi[2]

	#Radius
	r = 'radius'
	p = open(param,'r')
	radiu = grep(r,p)
	radi = radiu.split()
	rad = radi[2]

	p2 = open(param,'r')
	#SNR
	n = 'snr'
	nr = grep(n,p2)
	rn = nr.split()
	snr = rn[2]

	p3 = open(param,'r')
	#ts
	ts1 = 'ts'
	ts2 = grep(ts1,p3)
	ts3 = ts2.split()
	ts = ts3[2]

	#Box size
	p4 = open(param,'r')
	bxsz = 'boxSize'
	bxs = grep(bxsz,p4)
	bx = bxs.split()
	box = bx[2]

	p5 = open(param,'r')
	#Number of particles
	nmpts = 'num_part'
	nmpt = grep(nmpts,p5)
	nmp = nmpt.split()
	tot = nmp[2]

	#CS
	p6 = open(param,'r')
	cs1 = 'cs'
	cs2 = grep(cs1,p6)
	cs3 = cs2.split()
	cs = cs3[2]

	#Accelerating voltage
	p7 = open(param,'r')
	v1 = 'volt'
	v2 = grep(v1,p7)
	v3 = v2.split()
	volt = v3[2]

	#Free hand angular search
	p8 = open(param,'r')
	fs1 = 'freeHand_ang_search'
	fs2 = grep(fs1,p8)
	fs3 = fs2.split()
	angSearch = fs3[2]
	
	#Free hand Low resolution  
	p9 = open(param,'r')
	mr1 = 'min_res'
	mr2 = grep(mr1,p9)
	mr3 = mr2.split()
	min_res = mr3[2]

	#Free hand Max resolution  
	p10 = open(param,'r')
	mr4 = 'min_res'
	mr5 = grep(mr4,p10)
	mr6 = mr5.split()
	max_res = mr6[2]

	#Free hand first particle
	p11 = open(param,'r')
	fir1 = 'first'
	fir2 = grep(fir1,p11)
	fir3 = fir2.split()
	first = fir3[2]

	#Free hand last particle
	p12 = open(param,'r')
	ls1 = 'last'
	ls2 = grep(ls1,p12)
	ls3 = ls2.split()
	last = ls3[2]

	#Free hand Max resolution  
	p10 = open(param,'r')
	mr4 = 'max_res'
	mr5 = grep(mr4,p10)
	mr6 = mr5.split()
	max_res = mr6[2]

        #Free hand increment  
        p13 = open(param,'r')
        inc1 = 'incr'
        inc2 = grep(inc1,p13)
        inc3 = inc2.split()
        incr = inc3[2]

        #Free hand increment  
        p14 = open(param,'r')
        m1 = 'mag'
        m2 = grep(m1,p14)
        m3 = m2.split()
        mag = m3[2]

        #Free hand increment  
        p15 = open(param,'r')
        m1 = 'num_mod'
        m2 = grep(m1,p15)
        m3 = m2.split()
        num_mod = int(m3[2])

	p17 = open(param,'r')
        pp1 = 'cutoff'
        pp2 = grep(pp1,p17)
        pp3 = pp2.split()
        cutoff = pp3[2]

        p18 = open(param,'r')
        pp1 = 'calc'
        pp2 = grep(pp1,p17)
        pp3 = pp2.split()
        calc = pp3[2]

		
	#Prepare stack for EMAN2 refinement 
        print '\n'
        print 'Converting stack into EMAN2 format'
        print '\n'

	#Filter particles to specified resolution limits
                
	cmd = 'proc2d %s %s_prep.img apix=%s hp=%s lp=%s' %(untilt,untilt[:-4],pix,min_res,max_res)
        subprocess.Popen(cmd,shell=True).wait()

	if debug is True:
		print '%s/up_head.py %s_prep.img %s' %(cwd,untilt[:-4],pix)

       	cmd = '%s/up_head.py %s_prep.img %s' %(cwd,untilt[:-4],pix)
        subprocess.Popen(cmd,shell=True).wait()

        #Run refinement
        print '\n'                
	print 'Running EMAN2 refinement'                
	print '         Angular step = %s' %(ang)                
	print '         Shift range = %s' %(sx)                
	print '         Shift step size (ts)  = %s' %(ts)                
	print '         Pixel Size = %s' %(pix)                
	print '         Radius = %s' %(rad)                
	print '         SNR = %s' %(snr)
	print '	        CC_cut = %s' %(cutoff)                
	print '\n'                

	if num_mod == 1:
		
		if debug is True:
			print '%s/run_freehand.sh start.hdf %s %s %s %s %s %s %s %s' %(cwd,model,str(sx),str(ang),str(rad),str(snr),str(ts),str(cutoff),cwd)

		cmd = '%s/run_freehand.sh start.hdf %s %s %s %s %s %s %s %s' %(cwd,model,str(sx),str(ang),str(rad),str(snr),str(ts),str(cutoff),cwd)
               	subprocess.Popen(cmd,shell=True).wait()
	else:
		if debug is True:
			print '%s/run_freehand_sort.sh start.hdf %s %s %s %s %s %s %s %s' %(cwd,model,sx,ang,rad,snr,ts,cutoff,cwd)
	
		cmd = '%s/run_freehand_sort.sh start.hdf %s %s %s %s %s %s %s %s' %(cwd,model,sx,ang,rad,snr,ts,cutoff,cwd)
		subprocess.Popen(cmd,shell=True).wait()		

if __name__ == "__main__":     
	getEMANPath()             
	getOPENMPIPath()
	from EMAN2 import *     
	from sparx  import *     
	params=setupParserOptions()     
	main2(params)

