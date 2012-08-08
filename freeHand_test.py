#!/usr/bin/env python

import optparse
from sys import *
import os,sys,re
from optparse import OptionParser
import glob
import subprocess
from os import system
import linecache
import time

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
	parser.add_option("-x",action="store_true",dest="xmipp",default=False,
                help="Flag if data came from xmipp alignment")
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

def getCCP4Path():        
        ### get the openmpi directory        
        ccp4path = subprocess.Popen("env | grep CCP4_PATH", shell=True, stdout=subprocess.PIPE).stdout.read().strip()        

        if ccp4path:                
                ccp4path = ccp4path.replace("CCP4_PATH=","")      
        if os.path.exists(ccp4path):                        
                return ccp4path        
        print "ccp4 is not loaded, make sure it is in your path"        
        sys.exit()

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

def peak(stack,tot,cent):
	spifile = "currentSpiderScript.spi"
       	if os.path.isfile(spifile):
        	os.remove(spifile)
       	spi=open(spifile,'w')
	spicmd="SD IC NEW\n"
	spicmd+="incore\n"
	spicmd+="2,%s\n" %(tot)
	spicmd+="do lb1 [part] = 1, %s\n" %(tot)
	spicmd+="\n"
	spicmd+="PK [x] [y]\n"
	spicmd+="%s@{******[part]}\n"%(stack[:-4])
	spicmd+="(1,0)\n"
	spicmd+="[newX] = %s +[x]\n" %(cent)
	spicmd+="[newY] = %s +[y]\n" %(cent)
	spicmd+="SD IC [part] [newX] [newY]\n"
	spicmd+="incore\n"
	spicmd+="lb1\n"
	spicmd+="SD IC COPY\n"
	spicmd+="incore\n"
	spicmd+="%s_peak\n" %(stack[:-4])
	spicmd+="SD ICE\n"
	spicmd+="incore\n"
	spicmd+="en d\n"
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
	info = linecache.getline(ctf,5)                

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
	debug = params['debug']

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
	
	if debug is True:
		print 'Free-hand test completed for all particles'

	#Clean up:
	cmd = 'rm iteration?_finished tmp* fastfreehand_v1_01.exe'
	subprocess.Popen(cmd,shell=True)

def plotFH(params,ccp4_path):
	
	param = params['param']
	debug = params['debug']
	model = params['model']
	xmipp = params['xmipp']

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

	#Merge stacks:
	m = 0

	script = sys.argv[0]
	cwd = '%s/lib' %(script[:-17])
	
	num = len(glob.glob('model00_plots*.mrc'))

	i = 1
	while i <= int(num):

		cmd = 'e2proc2d.py model00_plots_CC_v101_%02d.mrc model00_plots_CC_v101_%02d.img --threed2twod' %(i,i)
		if debug is True:
			print cmd
		subprocess.Popen(cmd,shell=True).wait()

		cmd = 'rm model00_plots_CC_v101_%02d.mrc' %(i)
		subprocess.Popen(cmd,shell=True).wait()
	
		cmd = 'proc2d model00_plots_CC_v101_%02d.img model00_plots_CC_v101_merge.img' %(i)
		if debug is True:
			print cmd
		subprocess.Popen(cmd,shell=True).wait()

		i = i + 1

	if xmipp is True:
		

	cmd = 'e2proc2d.py model00_plots_CC_v101_merge.img model00_plots_CC_v101_merge.mrc --twod2threed'
	if debug is True:
		print cmd
	subprocess.Popen(cmd,shell=True).wait()
	
        cmd = 'cp %s/totsumstack.exe .' %(cwd)
        if debug is True:
                print cmd
	subprocess.Popen(cmd,shell=True).wait()

	totsum = '#!/bin/csh\n'
	totsum += 'totsumstack.exe << eot\n'
	totsum += 'model00_plots_CC_v101_merge.mrc\n' 
	totsum += 'model00_averageplot_CC_v101.mrc\n'
	totsum += 'eot\n'

        tmp = open('tmp.csh','w')
        tmp.write(totsum)
        tmp.close()

	if debug is True:
                print totsum

        cmd = 'chmod +x tmp.csh' 
        if debug is True:
                print cmd
	subprocess.Popen(cmd,shell=True)

        cmd = './tmp.csh' 
        subprocess.Popen(cmd,shell=True).wait()

        cmd = 'rm totsumstack.exe tmp.csh '
        subprocess.Popen(cmd,shell=True).wait()
	
	if calc is 'C':

		line1 = (float(angSearch)*2)/5
		line = line1/2
		if debug is True:
			print '%s = float(%s*2)/5' %(line1,angSearch)
			print '%s = %s / 2' %(line,line1)
			print '%s/npo_CC_wrap_mult.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)

		npo = '#!/bin/csh\n'
		npo += 'rm -f z.plot\n'
		npo += 'rm -f plot84.ps\n'
		npo += '%s/bin/npo mapin model00_averageplot_CC_v101.mrc plot z.plot << eof\n' %(ccp4_path)
		npo += 'NOTITLE\n'
		npo += 'MAP SCALE 1 INVERT\n'
		npo += '# For CCC\n'
		npo += 'CONTRS 0.0 to 1 by 0.002\n'
		npo += 'LIMITS 0 %s 0 %s 0 0\n' %(str(float(angSearch)*2),str(float(angSearch)*2))
		npo += 'SECTNS 0 0 1\n'
		npo += 'GRID  5 5\n'
		npo += 'GRID U DASHED 1.0 0.2 0 EVERY %s FULL\n' %(line)
		npo += 'GRID V DASHED 1.0 0.2 0 EVERY %s FULL\n' %(line)
		npo += 'PLOT Y\n'
		npo += 'eof\n'
		npo += '\n'
		npo += '%s/bin/pltdev -log -dev ps -abs -pen c -xp 3.1 -yp 3.1 -lan -i z.plot -o model00_average_frehand_CC.ps\n' %(ccp4_path)

	        tmp = open('tmp.csh','w')
        	tmp.write(npo)
	        tmp.close()

        	cmd = 'chmod +x tmp.csh'
	        subprocess.Popen(cmd,shell=True)

	        cmd = './tmp.csh'
	        subprocess.Popen(cmd,shell=True).wait()

        	cmd = 'rm tmp.csh z.plot'
	        subprocess.Popen(cmd,shell=True).wait()

		cmd = 'e2proc2d.py model00_plots_CC_v101_merge.img model00_plots_CC_v101_merge.spi' 
		subprocess.Popen(cmd,shell=True).wait()
	
		tot = EMUtil.get_image_count('model00_plots_CC_v101_merge.img') 	
		n = int(angSearch)+1
		stack = 'model00_plots_CC_v101_merge.spi' 
		peak(stack,tot,n)

		cmd = 'rm model00_plots_CC_v101_merge.img model00_plots_CC_v101_merge.hed model00_plots_CC_v101_merge.spi'
		subprocess.Popen(cmd,shell=True).wait()
	
	if calc is 'P':

                line1 = (float(angSearch)*2)/5
                line = line1/2
                if debug is True:
                        print '%s = float(%s*2)/5' %(line1,angSearch)
                        print '%s = %s / 2' %(line,line1)
                        print '%s/npo_CC_wrap_mult.csh %s model%02d %s' %(cwd,str(float(angSearch)*2),m,line)

                npo = '#!/bin/csh\n'
                npo += 'rm -f z.plot\n'
                npo += 'rm -f plot84.ps\n'
                npo += '%s/bin/npo mapin model00_averageplot_CC_v101.mrc plot z.plot << eof\n' %(ccp4_path)
                npo += 'NOTITLE\n'
                npo += 'MAP SCALE 1 INVERT\n'
                npo += '# For Pres\n'
                npo += 'CONTRS 77. to 86. by .3\n'
                npo += 'LIMITS 0 %s 0 %s 0 0\n' %(str(float(angSearch)*2),str(float(angSearch)*2))
                npo += 'SECTNS 0 0 1\n'
                npo += 'GRID  5 5\n'
                npo += 'GRID U DASHED 1.0 0.2 0 EVERY %s FULL\n' %(line)
                npo += 'GRID V DASHED 1.0 0.2 0 EVERY %s FULL\n' %(line)
                npo += 'PLOT Y\n'
                npo += 'eof\n'
		npo += '\n'
		npo += '%s/bin/pltdev -log -dev ps -abs -pen c -xp 3.1 -yp 3.1 -lan -i z.plot -o model00_average_frehand_CC.ps\nn' %(ccp4_path)

                tmp = open('tmp.csh','w')
                tmp.write(npo)
                tmp.close()

                cmd = 'chmod +x tmp.csh'
                subprocess.Popen(cmd,shell=True)

                cmd = './tmp.csh'
                subprocess.Popen(cmd,shell=True).wait()

                cmd = 'rm tmp.csh z.plot'
                subprocess.Popen(cmd,shell=True).wait()

                cmd = 'e2proc2d.py model00_plots_CC_v101_merge.img model00_plots_CC_v101_merge.spi'
                subprocess.Popen(cmd,shell=True).wait()

                tot = EMUtil.get_image_count('model00_plots_CC_v101_merge.img')  
                n = int(angSearch)+1
                stack = 'model00_plots_CC_v101_merge.spi' 
                peak(stack,tot,n)

		cmd = 'rm model00_plots_CC_v101_merge.img model00_plots_CC_v101_merge.hed model00_plots_CC_v101_merge.spi model00_plots_CC_v101_merge.mrc model00_plots*.img model00_plots*.hed'
                subprocess.Popen(cmd,shell=True).wait()

if __name__ == "__main__":     
	getEMANPath()             
	ccp4 = getCCP4Path()
	from EMAN2 import *     
	from sparx  import *     
	params=setupParserOptions()     
	fastFree(params)
	wait(params)
	plotFH(params,ccp4)
