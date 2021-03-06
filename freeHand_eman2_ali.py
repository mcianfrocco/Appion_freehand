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
	parser.add_option("-t",dest="tilted",type="string",metavar="FILE",
                help="tilted stack (black particles in IMAGIC format)")
        parser.add_option("-m",dest="model",type="string",metavar="FILE",
                help="3D model(s) for alignment (Single SPIDER volume or multi-volume HDF file)")
        parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Parameter file with refinement info (free_param.par)")
        parser.add_option("-c",dest="ctf",type="string",metavar="FILE",
		help="CTF-information file for tilted particles; DEFAULT = -2.00 um")
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
		print "\nWarning: no untilted stack specified\n"
	elif not os.path.exists(params['untilted']):
		print "\nError: stack file '%s' does not exist\n" % params['untilted']
		sys.exit()
        if not params['tilted']:
                print "\nWarning: no tilted stack specified\n"
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
	if not os.path.isfile(params['ctf']):
                print "\nNo CTF-information specified for tilted stack; using 2 um as default\n"
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

#========================
def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

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
        print 'Converting stack into EMAN2 format'
        print '\n'

	#Filter particles to specified resolution limits
                
	cmd = 'proc2d %s %s_prep.img apix=%s hp=%s lp=%s' %(untilt,untilt[:-4],pix,min_res,max_res)
        subprocess.Popen(cmd,shell=True).wait()

       	cmd = '%s/EMAN2/up_head.py %s_prep.img %s' %(untilt[:-4],pix)
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

		cmd = 'mpirun -np 8 %s/refine.py start.hdf %s refine_eman2 --ou=%s --rs=1 --xr=%s --ts=%s --delta=%s --snr=%s --center=0 --maxit=1 --ref_a=S --sym=c1 --cutoff=%s --MPI --full_output' %(cwd,model,rad,sx,ts,ang,snr,cutoff)
		
		if debug is True:
			print cmd
               	subprocess.Popen(cmd,shell=True).wait()
	else:
	
		cmd = 'mpirun -np 8 %s/refine.py start.hdf %s refine_eman2 --ou=%s --rs=1 --xr=%s --ts=%s --delta=%s --snr=%s --center=0 --maxit=1 --ref_a=S --sym=c1 --cutoff=%s --MPI --full_output --sort' %(cwd,model,rad,sx,ts,ang,snr,cutoff)

                if debug is True:
                        print cmd
		subprocess.Popen(cmd,shell=True).wait()		
	
	#Clean up:
	cmd = 'rm logfile* start.hdf %s_prep.*' %(untilt[:-4])
	subprocess.Popen(cmd,shell=True).wait()

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
			parmPSI = float(l[1])
			parmTHETA = float(l[2])
			parmPHI = float(l[3])
			sx =(float(l[4]))
			sy =(float(l[5]))
			model1 = float(l[6])
			#Convert euler angles from EMAN2 to FREALIGN/SPIDER convention
			if debug is True:
				print 'parmPSI = %s	parmTHETA = %s	parmPHI = %s	' %(parmPSI,parmTHETA,parmPHI)
				
			psi,theta,phi = Eman2Freali(parmPSI,parmTHETA,parmPHI)	
			out.write("%s 	%s	%s	%s	%s	%s\n"%(psi,theta,phi,sx,sy,model1))
		
		f.close()
		out.close()

		makeFH_eman2('%s_freeHand' %(parm),'%s_model%02d.par' %(ctf[:-4],int(mod_count)),mag,1,debug)

		eman2_mods(num_mod,model,mod_count,debug)

		im_to_mrc('%s_%02d.img' %(tilt[:-4],mod_count),debug)

		mod_count = mod_count + 1

#=================
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

#============
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
def makeFH_eman2(f,c,mag,div,debug):

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

                if debug is True:
                        print line
                psi = float(l[0])
                theta = float(l[1])
                phi = float(l[2])

                shiftx = float(l[3])/float(div)
                shifty = float(l[4])/float(div)

                ctf2 = linecache.getline(c,count)
                ctf = ctf2.split()
                df1 = float(ctf[0])
                df2 = float(ctf[1])
                astig = float(ctf[2])

                o1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f\n" %(count,psi,theta,phi,shiftx,shifty,float(mag),1,df1,df2,astig,50))

                count = count + 1

        o1.write("C\n")

#=========
def eman2_conv(params,paramout):

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

if __name__ == "__main__":     
	getEMANPath()             
	getOPENMPIPath()
	from EMAN2 import *     
	from sparx  import *     
	params=setupParserOptions()     
	checkConflicts(params)
	params.add_option('--dir',
	align(params)
	eman2_conv(params)
