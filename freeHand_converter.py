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
                help="3D model for alignment (Single MRC volume)")
        parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Parameter file with refinement info (free_param.par)")
        parser.add_option("-c",dest="ctf",type="string", metavar="FILE",
                help="Per-particle ctf-information file")
        parser.add_option("-o",dest="align",type="string", metavar="FILE",
                help="File with alignment info")
        parser.add_option("--prog",dest="prog",type="string", metavar="FILE",
                help="Program for 3D alignment: eman2, [other 3D programs], etc.")
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
def eman2_sort(paramout,tilt,ctf,num_mod):

	#Sort particles by model(s)	
	if num_mod == 1:

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
				y1.write('%s' %(line))				
				o1.write('%s' %(c))

 		  	count=count+1

      		text.close()
       		param.close()
	        cmd="proc2d %s %s_%02d.img list=%s_%02d.txt" %(stack,new,0,new,0)
       		subprocess.Popen(cmd,shell=True).wait()

	else:

		for n in range(0,numMods):
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
			cmd="proc2d %s %s_%02d.img list=%s_%02d.txt " %(stack,new,n,new,n)
			subprocess.Popen(cmd,shell=True).wait()

def eman2_angConv(paramout,num_mod,ctf,mag):
	mod_count = 0
	
	while mod_count < num_mod:

		print 'Working on model %s' %(mod_count)
			
		print '\n'                
		print 'Converting files into free-hand format'                
		print '\n'                

		parm='%s_model%02d' %(paramout)

		f=open(parm,'r')
		out = open("%s_freeHand"%(parm),'w')
		count=1
		count2=1
		count=1
		for line in f:
					
			l = line.split()
		
			parmPSI = float(l[0])
			parmTHETA = float(l[1])
			parmPHI = float(l[2])
			sx =(float(l[3]))
			sy =(float(l[4]))
			model = float(l[5])
			#Convert euler angles from EMAN2 to FREALIGN/SPIDER convention
			psi,theta,phi = Eman2Freali(parmPSI,parmTHETA,parmPHI)	
			out.write("%s 	%s	%s	%s	%s	%s\n"%(psi,theta,phi,sx,sy,model))
		
		f.close()
		out.close()

		#Convert parameter file format with CTF and angular info                
		#cmd = '%s/make_freeHand_Param.py refine_eman2/%s_freeHand %s_model%02d.par %s' %(cwd,paramout,ctf[:-4],mod_count,mag)                
		eman2_makeFH('%s_freeHand' %(paramout),'%s_model_%02d.par' %(ctf,mod_count),mag)

		eman2_mods(num_mod,model,mod_count,debug)

		eman2_parts(tilt,mod_count,debug)

def eman2_parts(tilt,mod_count,debug):
        #Convert tilted particles to 3D-MRC format                

        cmd = '%s/im_to_mrc.py -f %s_%02d.img ' %(cwd,tilt[:-4],mod_count)
	if debug is True:
		print cmd
        subprocess.Popen(cmd,shell=True).wait()

def eman2_mods(num_mod,model,mod_count,debug):

        #Convert model from HDF to MRC                

        if num_mod > 1:

                cmd = '%s/e2proc3d_all.py %s %s' %(cwd,model,mod_count)
                subprocess.Popen(cmd,shell=True).wait()

                cmd = 'rm %s_%03d.spi' %(model[:-4],mod_count)
                subprocess.Popen(cmd,shell=True).wait()

                cmd = 'mv %s_%03d_mr.spi %s_%03d.spi' %(model[:-4],mod_count,model[:-4],mod_count)
                subprocess.Popen(cmd,shell=True).wait()

	else:
                if debug is True:
        	        print 'cp %s %s_%03d.spi' %(model,model[:-4],mod_count)

                cmd = 'cp %s %s_%03d.spi' %(model,model[:-4],mod_count)
                subprocess.Popen(cmd,shell=True).wait()


        cmd = 'proc3d %s_%03d.spi %s_%03d.mrc' %(model[:-4],mod_count,model[:-4],mod_count)
        subprocess.Popen(cmd,shell=True).wait()

        tot = file_len('refine_eman2/%s' %(paramout))

        #Convert tilted particles to 3D-MRC format                

        if debug is True:
        	print '%s/imagic_to_freeHand2.py -f %s_%02d.img --total=%s --box=%s' %(cwd,tilt[:-4],mod_count,tot,box)

        cmd = '%s/im_to_mrc.py -f %s_%02d.img ' %(cwd,tilt[:-4],mod_count)
        subprocess.Popen(cmd,shell=True).wait()

#==================
def eman2_makeFH(f,c,mag):

	#Convert parameter file format with CTF info
	f1 = open(f,'r')
	fout = '%s_format' %(f)
	o1 = open(fout,'a')

	o1.write("C Frealign format parameter file created from Search_fspace parameter file\n")
	o1.write("C\n")
	o1.write("C           PSI   THETA     PHI     SHX     SHY    MAG   FILM      DF1      DF2  ANGAST  CCMax\n")

	mag = float(sys.argv[3])

	count = 1

	for line in f1:

		l = line.split()
	
		psi = float(l[0])
		theta = float(l[1])
		phi = float(l[2])

		shiftx = float(l[3])
		shifty = float(l[4])

		ctf2 = linecache.getline(c,count)
		ctf = ctf2.split()
		df1 = float(ctf[0])
		df2 = float(ctf[1])
		astig = float(ctf[2])

		o1.write("%7d%8.3f%8.3f%8.3f%8.3f%8.3f%8.0f%6d%9.1f%9.1f%8.2f%7.2f\n" %(count,psi,theta,phi,shiftx,shifty,mag,1,df1,df2,astig,50))

		count = count + 1

	o1.write("C\n")

#========================
def eman2(params):

	param = params['params']

	#Get current working directory
	script = sys.argv[0]
	cwd = '%s/lib' %(script[:-22])

	#Get parameter info: number of models
	p = open(param,'r')
	a = 'num_mods' 
	angl = grep(a,p)
	aL = angl.split()
	num_mod = aL[2]
	
	#Get parameter info: mag
        p = open(param,'r')
        a = 'mag'
        angl = grep(a,p)
        aL = angl.split()
        mag = aL[2]

	tilt = params['tilt']
	ctf = params['ctf']
	paramout = params['align']
	debug = params['debug']
	
	#Sort particles based up membership to model(s)
	eman2sort(paramout,tilt,ctf,num_mod)
	
	#Convert euler angles, model, and particles from EMAN2 to FREALIGN for each model
	eman2_angConv(paramout,num_mod,ctf,mag,model,tilt)

               
if __name__ == "__main__":     
	getEMANPath()             
	from EMAN2 import *     
	from sparx  import *     
	params=setupParserOptions()     
	checkConflicts(params)
	
	if params['prog'] is 'eman2':
		eman2(params)
