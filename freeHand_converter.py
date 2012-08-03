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
        parser.set_usage("%prog -t <untilted stack> -m <model> -p <parameter file> -c <ctf> -o <parms> --prog=NAME")
        parser.add_option("-t",dest="tilted",type="string",metavar="FILE",
                help="tilted stack (black, raw particles in IMAGIC format)")
        parser.add_option("-m",dest="model",type="string",metavar="FILE",
                help="3D model for used in alignment (Single SPIDER or MRC volume, or multi-volume HDF)")
        parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Parameter file with refinement info (free_param.par)")
        parser.add_option("-c",dest="ctf",type="string", metavar="FILE",
                help="Per-particle ctf-information file")
        parser.add_option("-o",dest="align",type="string", metavar="FILE",
                help="File with alignment info")
        parser.add_option("--prog",dest="prog",type="string", metavar="FILE",
                help="Program for 3D alignment: eman2, frealign")
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
	if not params['tilted']:
		print "\nWarning: no stack specified\n"
	elif not os.path.exists(params['tilted']):
		print "\nError: stack file '%s' does not exist\n" % params['tilted']
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
        if not params['ctf']:
                print "\nError: no ctf file specified"
                sys.exit()
        if not os.path.isfile(params['ctf']):
                print "\nError: ctf file %s does not exist\n" %params['ctf']
                sys.exit()
	#Check if output file exists
	f = params['align']
	test = '%s_format.par' %(f[:-4])
	if os.path.isfile(test):
                print "\nError: Output file %s already exists\n" %test
                sys.exit()
	

#========================
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#========================
def getEMANPath():        
        ### get the imagicroot directory        
        emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        

        if emanpath:                
                emanpath = emanpath.replace("EMAN2DIR=","")                
        if os.path.exists(emanpath):                        
                return emanpath        
        print "EMAN2 was not found, make sure eman2/2.05 is in your path"        
        sys.exit()

#========================
def Eman2Freali(az,alt,phi):

    t1 = Transform({"type":"eman","az":az,"alt":alt,"phi":phi,"mirror":False})
    d = t1.get_params("eman")   
    psi = d["phi"]+90   
    if psi >360:        
        psi = psi-360
    theta= d["alt"]
    phi = d["az"]-90
    return psi,theta,phi

#========================
def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

#========================
def eman2_sort(paramout,tilt,ctf,num_mod,debug):

	if debug is True:
		print 'eman2_sort():'
		print '		paramout = %s; tilt=%s; ctf=%s; num_mod=%s; debug=%s' %(paramout,tilt,ctf,num_mod,debug)

	#Sort particles by model(s)	
	if int(num_mod) == 1:
	
		if debug is True:
			print 'num_mod == 1'
       		param=open(paramout,'r')
       		count=1
       		text='%s_%02d.txt' %(tilt[:-4],0)
		c_o = '%s_model00.par' %(ctf[:-4])
		o1 = open(c_o,'w')
		y_o = '%s_model00' %(paramout)
		y1 = open(y_o,'w')

       		text=open(text,'w')

       		for line in param:
                	l=line.split()
                 	member=float(l[5])
			if debug is True:
				print l
                 	if member == 999:

                        	text.write("%s\n" %(count-1))
				
				c = linecache.getline(ctf,count)				
				y1.write('%s %s' %(str(count),line))				
				o1.write('%s' %(c))

 		  	count=count+1

      		text.close()
       		param.close()
	        cmd="proc2d %s %s_%02d.img list=%s_%02d.txt" %(tilt,tilt[:-4],0,tilt[:-4],0)
       		subprocess.Popen(cmd,shell=True).wait()

	else:

		for n in range(0,int(num_mod)):
			param=open(paramout,'r')
			c_o = '%s_model%02d.par' %(ctf[:-4],n)
                	o1 = open(c_o,'w')
			count=1
                	y_o = '%s_model%02d' %(paramout,n)
                	y1 = open(y_o,'w')
			text='%s_%02d.txt' %(tilt[:-4],n)	

			text=open(text,'w')

			for line in param:

				l=line.split()
				member=float(l[5])
	
				if member == n:

					text.write("%s\n" %(count-1))
					c = linecache.getline(ctf,count)
					y1.write('%s' %(line))
                         	        o1.write('%s' %(c))

				count=count+1
			text.close()
			param.close()
			cmd="proc2d %s %s_%02d.img list=%s_%02d.txt " %(tilt,tilt[:-4],n,tilt[:-4],n)
			subprocess.Popen(cmd,shell=True).wait()

def eman2_angConv(paramout,num_mod,ctf,mag,model,tilt,debug):
	mod_count = 0
	
	while mod_count < int(num_mod):

		print 'Working on model %s' %(mod_count)
			
		print '\n'                
		print 'Converting files into free-hand format'                
		print '\n'                

		parm='%s_model%02d' %(paramout,mod_count)
		if debug is True:
			print 'parm = %s' %(parm)

		f=open(parm,'r')
		out = open("%s_freeHand"%(parm),'w')
		count=1
		count2=1
		count=1
		for line in f:
					
			l = line.split()
			if debug is True:
				print l		
			parmPSI = float(l[0])
			parmTHETA = float(l[1])
			parmPHI = float(l[2])
			sx =(float(l[3]))
			sy =(float(l[4]))
			model1 = float(l[5])
			#Convert euler angles from EMAN2 to FREALIGN/SPIDER convention
			psi,theta,phi = Eman2Freali(parmPSI,parmTHETA,parmPHI)	
			out.write("%s 	%s	%s	%s	%s	%s\n"%(psi,theta,phi,sx,sy,model1))
		
		f.close()
		out.close()

		makeFH('%s_freeHand' %(parm),'%s_model%02d.par' %(ctf[:-4],int(mod_count)),mag,1,debug)

		eman2_mods(num_mod,model,mod_count,debug)

		im_to_mrc('%s_%02d.img' %(tilt[:-4],mod_count),debug)

		mod_count = mod_count + 1

def im_to_mrc(stack,debug):

        #Convert tilted particles to 3D-MRC format                

        # get box size
        im=EMData.read_images(stack,[0])
        nx = im[0].get_xsize()
        del im
	nimg = EMUtil.get_image_count(stack)

        img = EMData(nx,nx,nimg)
        img.write_image(stack[:-4]+'.mrc')

	i = 0

        while i < nimg:
                d = EMData()
                d.read_image(stack, i)
                region = Region(0, 0, i, nx, nx, 1)
                d.write_image(stack[:-4]+".mrc",0,EMUtil.get_image_ext_type("mrc"), False, region, EMUtil.EMDataType.EM_FLOAT, True)
        	i = i + 1

def eman2_mods(num_mod,model,mod_count,debug):

        #Convert model from HDF to MRC                
	if debug is True:
		print num_mod
		print model
		print mod_count
	
        if int(num_mod) > 1:

                cmd = 'e2proc3d.py --first=%s --last=%s %s %s_%03d.mrc' %(model,model[:-4],mod_count)
        	if debug is True:
			print cmd
	        subprocess.Popen(cmd,shell=True).wait()

	else:

                cmd = 'proc3d %s %s_%03d.mrc' %(model,model[:-4],int(mod_count))
		if debug is True:
			print cmd
                subprocess.Popen(cmd,shell=True).wait()

#==================
def makeFH(f,c,mag,div,debug):

	#Convert parameter file format with CTF info
	f1 = open(f,'r')
	fout = '%s_format.par' %(f[:-4])
	o1 = open(fout,'a')
	if debug is True:
		print 'c = %s' %(c)
	o1.write("C Frealign format parameter file created from Search_fspace parameter file\n")
	o1.write("C\n")
	o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

	count = 1

	for line in f1:

		l = line.split()
	
		if l[0] is 'C':
			continue
		if debug is True:
			print line
		psi = float(l[1])
		theta = float(l[2])
		phi = float(l[3])

		shiftx = float(l[4])/float(div)
		shifty = float(l[5])/float(div)

		ctf2 = linecache.getline(c,count)
		ctf = ctf2.split()
		df1 = float(ctf[0])
		df2 = float(ctf[1])
		astig = float(ctf[2])

		o1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f\n" %(count,psi,theta,phi,shiftx,shifty,float(mag),1,df1,df2,astig,50))

		count = count + 1

	o1.write("C\n")

#========================
def eman2(params):

	param = params['param']

	#Get current working directory
	script = sys.argv[0]
	cwd = '%s/lib' %(script[:-22])

	#Get parameter info: number of models
	p = open(param,'r')
	a = 'num_mod' 
	angl = grep(a,p)
	aL = angl.split()
	num_mod = aL[2]
	
	#Get parameter info: mag
        p = open(param,'r')
        a = 'mag'
        angl = grep(a,p)
        aL = angl.split()
        mag = aL[2]

	tilt = params['tilted']
	ctf = params['ctf']
	paramout = params['align']
	debug = params['debug']
	model = params['model']
	
	#Sort particles based up membership to model(s)
	eman2_sort(paramout,tilt,ctf,num_mod,debug)
	
	#Convert euler angles, model, and particles from EMAN2 to FREALIGN for each model
	eman2_angConv(paramout,num_mod,ctf,mag,model,tilt,debug)

	#Clean up
	mod = 0
	while mod < int(num_mod):
		
		cmd = 'rm %s_model%02d %s_model%02d_freeHand %s_model%02d.par %s_%02d.img %s_%02d.hed %s_%02d.txt' %(paramout,mod,paramout,mod,ctf[:-4],mod,tilt[:-4],mod,tilt[:-4],mod,tilt[:-4],mod)
		if debug is True:
			print cmd
		subprocess.Popen(cmd,shell=True).wait()
		
		mod = mod + 1
               
def frealign(params):

        param = params['param']

        #Get current working directory
        script = sys.argv[0]
        cwd = '%s/lib' %(script[:-22])

        #Get parameter info: mag
        p = open(param,'r')
        a = 'mag'
        angl = grep(a,p)
        aL = angl.split()
        mag = aL[2]

       #Get parameter info: pixel size
        p = open(param,'r')
        a = 'pix'
        angl = grep(a,p)
        aL = angl.split()
        pix = aL[2]

        tilt = params['tilted']
        ctf = params['ctf']
        paramout = params['align']
        debug = params['debug']
        model = params['model']
	f = params['align'] 
	
	makeFH(f,ctf,mag,pix,debug)
	im_to_mrc(tilt,debug)

if __name__ == "__main__":     
	getEMANPath()             
	from EMAN2 import *     
	from sparx  import *     
	params=setupParserOptions()     
	checkConflicts(params)

	if params['prog'] == 'eman2':
		eman2(params)

	if params['prog'] == 'frealign':
		frealign(params)

	if params['prog'] != 'eman2' and params['prog'] != 'frealign':
		print 'prog=%s unknown option specified' %(params['prog'])
