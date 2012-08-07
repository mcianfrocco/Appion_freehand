#!/usr/bin/env python

import optparse
from sys import *
import os,sys,re
from optparse import OptionParser
import glob
import subprocess
from os import system
import linecache

#=========================
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

#=========================
def checkConflicts(params):
	if not params['untilted']:
		print "\nWarning: no stack specified\n"
	elif not os.path.exists(params['untilted']):
		print "\nError: stack file '%s' does not exist\n" % params['untilted']
		sys.exit()
        if not params['model']:
                print "\nWarning: no model specified\n"
        elif not os.path.exists(params['model']):
                print "\nError: model file '%s' does not exist\n" % params['model']
                sys.exit()
	if not params['param']:
		print "\nError: no free_param.par file specified"
		sys.exit()
	if not os.path.isfile(params['param']):
		print "\nError: free_param.par file does not exist\n" 
		sys.exit()
	if os.path.exists('data'):
		print "\ndata/ folder already exists, please remove and restart\n"
		sys.exit()
	if os.path.isfile('data.sel'):
		print "\ndata.sel already exists, please remove and restart\n"
                sys.exit()
	if os.path.isfile('xmipp_alignment_info.doc'):
                print "\nxmipp_alignment_info.doc  already exists, please remove and restart\n"
                sys.exit()

#========================
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#========================
def getXMIPPPath():        
        ### get the imagicroot directory        
        xpath = subprocess.Popen("env | grep xmipp-2.4-mpi", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        

        if len(xpath) == 0:
		print "XMIPP was not found, make sure xmipp-2.4-mpi is in your path"
	        sys.exit()
        	
#========================
def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

def im_to_xmipp(stack,tot):

	spicmd = 'VM\n'
	spicmd += 'mkdir data\n'
	spicmd += 'do lb1 [part]=1,%s\n' %(tot)
	spicmd += '	CP\n'
	spicmd += '	%s@{******[part]}\n' %(stack)
	spicmd += '	data/data{******[part]}\n'
	spicmd += 'lb1\n'
	runSpider(spicmd)

def runSpider(lines):
       spifile = "currentSpiderScript.spi"
       if os.path.isfile(spifile):
               os.remove(spifile)
       spi=open(spifile,'w')
       spi.write("MD\n")
       spi.write("TR OFF\n")
       spi.write("MD\n")
       spi.write("VB OFF\n")
       spi.write("MD\n")
       spi.write("SET MP\n")
       spi.write("(0)\n")
       spi.write("\n")
       spi.write(lines)

       spi.write("\nEN D\n")
       spi.close()
       spicmd = "spider spi @currentSpiderScript"
       spiout = subprocess.Popen(spicmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr.read()
       output = spiout.strip().split()
       if "ERROR" in output:
               print "Spider Error, check 'currentSpiderScript.spi'\n"
               sys.exit()
       # clean up
       os.remove(spifile)
       if os.path.isfile("LOG.spi"):
               os.remove("LOG.spi")
       resultf = glob.glob("results.spi.*")
       if resultf:
               for f in resultf:
                       os.remove(f)

#========================
def align(params):
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
        print 'Converting stack into XMIPP format'
        print '\n'

	#Filter particles to specified resolution limits
                
	cmd = 'proc2d %s %s_prep.img apix=%s hp=%s lp=%s' %(untilt,untilt[:-4],pix,min_res,max_res)
        subprocess.Popen(cmd,shell=True).wait()

	cmd = 'proc2d %s_prep.img %s_prep.spi spider' %(untilt[:-4],untilt[:-4])
	subprocess.Popen(cmd,shell=True).wait()

        cmd = 'rm %s_prep.img %s_prep.hed' %(untilt[:-4],untilt[:-4])
	subprocess.Popen(cmd,shell=True).wait()
	
	im_to_xmipp('%s_prep' %(untilt[:-4]),tot)
	
	#Create select file
	cmd = 'xmipp_selfile_create "data/*.spi" > data.sel'
	subprocess.Popen(cmd,shell=True).wait()

	#Run refinement
        print '\n'                
	print 'Running XMIPP refinement'                
	print '         Angular step = %s' %(ang)                
	print '         Shift range = %s' %(sx)                
	print '         Shift step size (ts)  = %s' %(ts)                
	print '         Pixel Size = %s' %(pix)                
	print '         Radius = %s' %(rad)                
	print '         SNR = %s' %(snr)
	print '\n'                

	cmd = 'xmipp_header_extract -i data.sel -o original_angles.doc'
	subprocess.Popen(cmd,shell=True).wait()

	cmd = 'mkdir ReferenceLibrary ProjMatchClasses'
	subprocess.Popen(cmd,shell=True).wait()

	cmd = 'mpirun -np %s xmipp_mpi_angular_project_library  -i %s -experimental_images original_angles.doc -o ReferenceLibrary/ref -sampling_rate %s -sym c1h -compute_neighbors  -angular_distance -1' %(incr,model,ang)
	if debug is True:
		print cmd
	subprocess.Popen(cmd,shell=True).wait()

	cmd = 'mpirun -np %s xmipp_mpi_angular_projection_matching  -i original_angles.doc -o ProjMatchClasses/proj_match -ref ReferenceLibrary/ref -Ri 0 -Ro %s -max_shift 1000 -search5d_shift %s -search5d_step  %s -mem 2 -thr 1 -sym c1h' %(incr,rad,sx,ts)
	if debug is True:
                print cmd
        subprocess.Popen(cmd,shell=True).wait()

	#cmd = 'xmipp_angular_class_average  -i ProjMatchClasses/proj_match.doc -lib ReferenceLibrary/ref_angles.doc -dont_write_selfiles  -limit0 -1 -limitR 10 -o ProjMatchClasses/proj_match -split'
	#if debug is True:
        #        print cmd
        #subprocess.Popen(cmd,shell=True).wait()

	#cmd="mpirun -np %s xmipp_mpi_reconstruct_wbp -i ProjMatchClasses/proj_match_classes.sel -o %s_ali.spi -threshold 0.02 -sym c1  -use_each_image -weight" %(incr,model[:-4])
	#if debug is True:
        #        print cmd
        #subprocess.Popen(cmd,shell=True).wait()

	#cmd="xmipp_fourier_filter -i %s_ali.spi -o %s_ali_filt25A.spi -low_pass 25 -sampling %s" %(model[:-4],model[:-4],pix)
	#if debug is True:
        #        print cmd
        #subprocess.Popen(cmd,shell=True).wait()

	cmd = 'cp ProjMatchClasses/proj_match.doc xmipp_alignment_info.doc'
	subprocess.Popen(cmd,shell=True).wait()

	#Clean up
	if debug is False:
	
		cmd = 'rm -r ReferenceLibrary ProjMatchClasses data/ %s_prep.spi original_angles.doc data.sel' %(untilt[:-4])
		subprocess.Popen(cmd,shell=True).wait()

if __name__ == "__main__":     
	getXMIPPPath()             
	params=setupParserOptions()     
	checkConflicts(params)
	align(params)
