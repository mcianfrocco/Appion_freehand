#!/usr/bin/env python



	paramout = 'paramout_%03d_00' %(int(ang))
		
	if debug is True:
		print'%s/sort_parts_by_mod.py -f %s -p refine_eman2/%s -c %s --num=%s' %(cwd,tilt,paramout,ctf,num_mod)

	cmd = '%s/sort_parts_by_mod.py -f %s -p refine_eman2/%s -c %s --num=%s' %(cwd,tilt,paramout,ctf,num_mod)
	subprocess.Popen(cmd,shell=True).wait()

	mod_count = 0
	
	while mod_count < num_mod:

		print 'Working on model %s' %(mod_count)
			
		print '\n'                
		print 'Converting files into free-hand format'                
		print '\n'                
		paramout = 'paramout_%03d_00_model%02d' %(int(ang),mod_count)                
	
		if debug is True:
			print '%s/EMAN_to_FREALIGN_fromParam_freeHand.py -p refine_eman2/%s' %(cwd,paramout)

		cmd = '%s/EMAN_to_FREALIGN_fromParam_freeHand.py -p refine_eman2/%s' %(cwd,paramout)                
		subprocess.Popen(cmd,shell=True).wait()                
		
		#Convert parameter file format with CTF and angular info                
		cmd = '%s/make_freeHand_Param.py refine_eman2/%s_freeHand %s_model%02d.par %s' %(cwd,paramout,ctf[:-4],mod_count,mag)                
		subprocess.Popen(cmd,shell=True).wait()                

		#Convert model from HDF to MRC                
				
		if num_mod > 1:

			if debug is True:
				print '%s/e2proc3d_all.py %s %s' %(cwd,model,mod_count)

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

               
		#Run Free-Hand test                
		info = linecache.getline('refine_eman2/%s_freeHand_format' %(paramout),4)                

		i = info.split()                

		#mag = i[6]                
		df1 = i[8]
	        df2 = i[9]
                astig = i[10]
    
	        i = 1
	        
		iteration = 1 

               	while i < int(tot): 
                	last = str(i + float(incr)-1)  
                        last = last[:-2] 
                        if i == 1:
                                 first = str(i) 
                        else:
                                 first = str(i)
	                         first = first[:-2]
                        if float(last) > int(tot): 
                                 incr = int(incr) - (int(last)- int(tot))
                                 last = str(tot)
			if last == int(tot):

				cmd= '%s/fastFreeHand_wrapper.csh %s %s %s %s %s_%02d.mrc %s_%03d.mrc refine_eman2/%s_freeHand_format %s %s %s %s %s %s %s %s %s %s %s %s model%02d %s' %(cwd,pix,snr,cs,volt,tilt[:-4],mod_count,model[:-4],mod_count,paramout,angSearch,min_res,max_res,str(float(pix)*float(rad)),first,last,incr,mag,df1,df2,astig,str(iteration),mod_count,calc)
                        	if debug is True:
					print cmd

				if debug is False:
      				        subprocess.Popen(cmd,shell=True).wait()

			else:

				if debug is True:
					print '%s/fastFreeHand_wrapper.csh %s %s %s %s %s_%02d.mrc %s_%03d.mrc refine_eman2/%s_freeHand_format %s %s %s %s %s %s %s %s %s %s %s %s model%02d %s' %(cwd,pix,snr,cs,volt,tilt[:-4],mod_count,model[:-4],mod_count,paramout,angSearch,min_res,max_res,str(float(pix)*float(rad)),first,last,incr,mag,df1,df2,astig,str(iteration),mod_count,calc)
				

                	        if debug is False:
					cmd= '%s/fastFreeHand_wrapper.csh %s %s %s %s %s_%02d.mrc %s_%03d.mrc refine_eman2/%s_freeHand_format %s %s %s %s %s %s %s %s %s %s %s %s model%02d %s' %(cwd,pix,snr,cs,volt,tilt[:-4],mod_count,model[:-4],mod_count,paramout,angSearch,min_res,max_res,str(float(pix)*float(rad)),first,last,incr,mag,df1,df2,astig,str(iteration),mod_count,calc)
                        		subprocess.Popen(cmd,shell=True)

       	                i = i + float(incr)
               	        iteration = iteration + 1


		mod_count = mod_count + 1

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


