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
        parser.set_usage("%prog -s <stack> -m <model> -p <parameter file> -c <ctf/alignment info>")
        parser.add_option("-s",dest="stack",type="string",metavar="FILE",
                help="Particle stack in Free-hand format (black,raw particles)")
        parser.add_option("-m",dest="model",type="string",metavar="FILE",
                help="3D model (Single MRC volume)")
        parser.add_option("-p",dest="param",type="string", metavar="FILE",
                help="Parameter file with refinement info (free_param.par)")
        parser.add_option("-c",dest="ctf",type="string", metavar="FILE",
                help="Per-particle CTF and alignment info (FREALIGN format)")
        parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 4:
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

def grep(string,list):
    expr = re.compile(string)
    for text in list:
        match = expr.search(text)
        if match != None:
            return match.string

def fastFree(params):
	debug = params['debug']
	param = params['param']
	stack = params['stack']
	model = params['model']
	ctf = params['ctf']

	#Get current working directory
	script = sys.argv[0]
	cwd = '%s/lib' %(script[:-17])

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
        procs = inc3[2]

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
        
	#Run Free-Hand test                
	info = linecache.getline(ctf,4)                

	i = info.split()                

	df1 = i[8]
        df2 = i[9]
        astig = i[10]
    
       	i = 1
	        
	iteration = 1 

	#Number of particles

	y = open(ctf,'r')
	tot = len(y.readlines())
	tot = tot - 4

	incr = str(round(tot/int(procs))+1)

	#Setup inputs for free-hand executable
	cmd = 'cp %s/fastfreehand_v1_01.exe .' %(cwd)					
	subprocess.Popen(cmd,shell=True).wait()
	exe = 'fastfreehand_v1_01.exe\n'
	p1 = '%s,%s,%s,%s\n' %(pix,snr,cs,volt)					
	p2 = '1,0\n'
	p3 = '%s\n' %(stack)
	p4 = '%s\n' %(model)
	p5 = '%s\n' %(ctf)
	#p6 = plots_.mrc (see below)
	p7 = '-%s,%s,-%s,%s\n' %(angSearch,angSearch,angSearch,angSearch)
	p8 = '%s,%s,%s\n' %(min_res,max_res,str(float(pix)*float(rad)))
	#p9 = first, last (see below)
	p10 = '%s\n' %(calc)
	p11 = '%s,%s,1,%s,%s,%s,0\n' %(incr[:-2],mag,df1,df2,astig)

        while i < int(tot): 
               	last = str(i + float(incr)-1)  
                last = last[:-2] 
                if i == 1:
        	        first = str(i) 
                else:
                        first = str(i)
	                first = first[:-2]
                if float(last) > int(tot): 
                        incr = int(incr[:-2]) - (int(last)- int(tot))
                        last = str(tot)
	
		p6 = 'model00_plots_CC_v101_%02d.mrc\n' %(iteration)
		p9 = '%s,%s\n' %(first,last)

		ff_cmd ='#!/bin/csh\n'
		ff_cmd +='fastfreehand_v1_01.exe << eot\n'
		ff_cmd +=p1
		ff_cmd +=p2
		ff_cmd +=p3
		ff_cmd +=p4
                ff_cmd +=p5
                ff_cmd +=p6
                ff_cmd +=p7
                ff_cmd +=p8
                ff_cmd +=p9
                ff_cmd +=p10
                ff_cmd +=p11
		ff_cmd +='eot\n'
		ff_cmd +='touch iteration%01d_finished\n' %(iteration)

		tmp = open('tmp%01d.csh'%(iteration),'w')
		tmp.write(ff_cmd)
		tmp.close()

		cmd = 'chmod +x tmp%01d.csh' %(iteration)
		subprocess.Popen(cmd,shell=True)
	
		cmd = './tmp%01d.csh' %(iteration)
		subprocess.Popen(cmd,shell=True)

		i = i + float(incr)
               	iteration = iteration + 1


def wait(params):

	param = params['param']

        #Free hand increment  
        p13 = open(param,'r')
        inc1 = 'incr'
        inc2 = grep(inc1,p13)
        inc3 = inc2.split()
        procs = inc3[2]

	i = 1
	
	while i<= int(procs):

		test = os.path.isfile('iteration%01d_finished'%(i))

		if test is False:
			time.sleep(5)
		if test is True:
			i = i + 1
	
	print 'Free-hand test completed for all particles'

	#Clean up:
	cmd = 'rm iteration?_finished tmp*'
	subprocess.Popen(cmd,shell=True)

def plot(params):
	
	param = params['param']
	out = params['out']
	debug = params['debug']
	model = params['model']
        #Free hand angular search
        p8 = open(param,'r')
        fs1 = 'freeHand_ang_search'
        fs2 = grep(fs1,p8)
        fs3 = fs2.split()
        angSearch = fs3[2]

        p18 = open(param,'r')
        fs1 = 'calc'
        fs2 = grep(fs1,p8)
        fs3 = fs2.split()
        calc = fs3[2]

        #Free hand angular search
        p9 = open(param,'r')
        fs1 = 'num_mod'
        fs2 = grep(fs1,p9)
        fs3 = fs2.split()
        mods = float(fs3[2])

	#Merge stacks:
	m = 0

	script = sys.argv[0]
	cwd = '%s/lib' %(script[:-22])
	
	while m < mods:

		num = len(glob.glob('model%02d_plots*.mrc'%(m)))

		i = 1
		while i <= int(num):

			cmd = '%s/mrc_to_im.b model%02d_plots_CC_v101_%s.mrc' %(cwd,m,i)
			if debug is True:
				print cmd
			subprocess.Popen(cmd,shell=True).wait()

			cmd = 'proc2d model%02d_plots_CC_v101_%s.img model%02d_plots_CC_v101_merge.img' %(m,i,m)
			if debug is True:
				print cmd
			subprocess.Popen(cmd,shell=True).wait()

			i = i + 1

		cmd = '%s/im_to_mrc.b model%02d_plots_CC_v101_merge.img' %(cwd,m)
		subprocess.Popen(cmd,shell=True).wait()

	        cmd = 'cp %s/totsumstack.exe .' %(cwd)
	        subprocess.Popen(cmd,shell=True).wait()

		cmd = 'cp %s/totsumstack_mult.csh .' %(cwd)
                subprocess.Popen(cmd,shell=True).wait()
	
	        cmd = './totsumstack_mult.csh model%02d' %(m) 
	        subprocess.Popen(cmd,shell=True).wait()

	        cmd = 'rm totsumstack.exe totsumstack_mult.csh'
	        subprocess.Popen(cmd,shell=True).wait()
	
		if calc is 'C':

			line1 = (float(angSearch)*2)/5
			line = line1/2
			if debug is True:
				print '%s = float(%s*2)/5' %(line1,angSearch)
				print '%s = %s / 2' %(line,line1)
				print '%s/npo_CC_wrap_mult.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)

		        cmd = '%s/npo_CC_wrap_mult.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)
		        subprocess.Popen(cmd,shell=True).wait()

			cmd = '%s/mrc_to_spi.b model%02d_plots_CC_v101_merge.mrc' %(cwd,m)
			subprocess.Popen(cmd,shell=True).wait()
	
			cmd = 'mkdir %s' %(out)
		       	subprocess.Popen(cmd,shell=True).wait()
	
			cmd = 'mv model%02d_plots_CC_v101* %s/' %(m,out)
		        subprocess.Popen(cmd,shell=True).wait()

			cmd = 'mv model%02d_frehand_CC.ps %s/' %(m,out)
		        subprocess.Popen(cmd,shell=True).wait()

	       		cmd = 'mv model%02d_averageplot_CC_v101.mrc %s/' %(m,out)
		        subprocess.Popen(cmd,shell=True).wait()

			cmd = 'cp %s %s/' %(model,out)
			subprocess.Popen(cmd,shell=True).wait()
	
		        cmd = 'cp %s %s/' %(param,out)
		        subprocess.Popen(cmd,shell=True).wait()
	
			tot = EMUtil.get_image_count('%s/model%02d_plots_CC_v101_merge.img' %(out,m)) 	
			n = int(angSearch)+1
			stack = '%s/model%02d_plots_CC_v101_merge.spi' %(out,m)
			peak(stack,tot,n)
			m = m + 1

                if calc is 'P':

                        line1 = (float(angSearch)*2)/5
                        line = line1/2
                        if debug is True:
                                print '%s = float(%s*2)/5' %(line1,angSearch)
                                print '%s = %s / 2' %(line,line1)
                                print '%s/npo_CC_wrap_mult_phase.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)

                        cmd = '%s/npo_CC_wrap_mult_phase.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = '%s/mrc_to_spi.b model%02d_plots_CC_v101_merge.mrc' %(cwd,m)
                        subprocess.Popen(cmd,shell=True).wait()

			cmd = 'mkdir %s' %(out)
                       	subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'mv model%02d_plots_CC_v101* %s/' %(m,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'mv model%02d_frehand_CC.ps %s/' %(m,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'mv model%02d_averageplot_CC_v101.mrc %s/' %(m,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'cp %s %s/' %(model,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'cp %s %s/' %(param,out)
                        subprocess.Popen(cmd,shell=True).wait()

                        tot = EMUtil.get_image_count('%s/model%02d_plots_CC_v101_merge.img' %(out,m))
                        n = int(angSearch)+1
                        stack = '%s/model%02d_plots_CC_v101_merge.spi' %(out,m)
                        peak(stack,tot,n)
                        m = m + 1

		cmd = 'mv *00.* %s' %(out)	
		subprocess.Popen(cmd,shell=True).wait()

	cmd = 'rm -r logfile* test.img test.hed model??*.mrc refine_eman2 z.plot start.hdf *_prep.img *_prep.hed '
 	subprocess.Popen(cmd,shell=True).wait()

	cmd = "cp %s/find_peaks_freeHand.spi %s" %(cwd,out)
	subprocess.Popen(cmd,shell=True).wait()

if __name__ == "__main__":     
	getEMANPath()             
	from EMAN2 import *     
	from sparx  import *     
	params=setupParserOptions()     
	fastFree(params)
