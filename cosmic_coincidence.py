import ROOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

f=ROOT.TFile("out.root")
t=f.Get("raw_wf_tree")
p=f.Get("pulse_tree")

sum_rms=0
nadded=0

num_ch=5

y=0

kernel=np.load("spe_kernel.npy")
newkernel=np.zeros(2000)
finalkernel=np.zeros(2000)
for x in range(0,len(kernel)-48):
  newkernel[x]=kernel[x+48]
#finalkernel[0:14]=newkernel[0:14]
finalkernel=newkernel
plt.plot(finalkernel)
plt.show()

fftker=np.fft.rfft(finalkernel)
print len(fftker)
epsilon = 0.05*np.amax(np.abs(fftker))
epvec=np.ones(1001)
epvec=np.multiply(epsilon,np.ones(len(fftker)))
kernel_fourier=np.add(fftker,epvec)


print kernel_fourier

#plt.plot(fftker,"r-")
#plt.plot(kernel_fourier,"b-")
plt.show()
print "SHEEP"
plt.plot(np.fft.irfft(kernel_fourier))
plt.show()

pulses=[]
num_summed=0
sum_wf=np.zeros(2000)
sum_deco=np.zeros(2000)

print "WFM:",t.GetEntries()
print "PULSE:",p.GetEntries()

print p.GetEntries()
for entry in xrange(0,p.GetEntries()/num_ch):
	print "Entry: ",entry,"  of: ",p.GetEntries()
	p.GetEvent(entry*num_ch)

	if p.npulses==0:
		print "PULSE0:",entry
		continue

	amps=np.asarray(p.amp)
	if np.argmax(amps)!=0:
		print "NO AMP0",entry
		continue
	if p.reject_edge==True:
		print "BADEDGE",entry
		continue
	
	ch4s=np.asarray(p.tstart)
	ch4e=np.asarray(p.tend)
	ped4=p.ped
	amp4=p.event_amp
	t4=p.tstart[0]
		
	p.GetEvent(entry*num_ch+1)

	if p.npulses==0:
		continue

	ch13s=np.asarray(p.tstart)
	ch13e=np.asarray(p.tend)
	ped13=p.ped
	t13=p.tstart[0]

	fflag=True
	
	for s1 in range(0,len(ch4s)):
		for s2 in range(0,len(ch13s)):
			if ((ch13s[s2]>=ch4s[s1])&(ch13s[s2]<=ch4e[s1])&fflag):
				print "GOOD",entry
				t.GetEntry(entry*5)
				vch4=np.asarray(t.wf)
				vch4=np.add(vch4,-1.0*ped4)

				p.GetEntry(entry*5)
				mystart=np.asarray(p.tstart)
				myend=np.asarray(p.tend)

				matplotlib.rcParams.update({'font.size': 18})
				'''
				plt.plot(vch4,'b-',label="Waveform")
				plt.plot(mystart,vch4[mystart],'go',label="Start")
				plt.plot(myend,vch4[myend],'ro',label="End")
				plt.ylabel("Amplitude (ADC)")
				plt.xlabel("Time (16ns ticks)")
				legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
				'''
				#plt.show()
				#plt.plot(vch4,'b-')
				#plt.show()
				#plt.plot(ch4s,vch4[ch4s],'go')
				#plt.plot(ch4e,vch4[ch4e],'ro')
				#plt.show()
				#plt.cla()
				#t.GetEntry(entry*5+1)
				#vch13=np.asarray(t.wf)
				#vch13=np.add(vch13,-1.0*ped13)
				#vch13[0:t13]=0
				#vch13=np.roll(vch13,-1*int(t13))
				#time=np.arange(0,len(vch4))
				#plt.plot(time,vch4,'r-')
				#plt.plot(time,vch13,'b-')
				#plt.show()
				#plt.cla()		
				fflag=False
		

				wf_fft=np.fft.rfft(vch4)
				#deco=np.fft.irfft(np.divide(wf_fft,kernel_fourier))
				deco=np.fft.irfft(np.divide(wf_fft,kernel_fourier))
				print "LAMB"	
				#plt.plot(np.divide(vch4,np.amax(vch4)),'b-')
				#plt.plot(np.divide(deco,np.amax(deco)),'r-')
				#plt.show()

				'''
				plt.plot(fax,np.divide(np.abs(wf_fft),np.amax(np.abs(wf_fft))),'bo',label="Signal")
				plt.plot(fax,np.divide(np.abs(kernel_fourier),np.amax(np.abs(kernel_fourier))),'ro',label="Kernel")
				plt.xlabel("Frequency (MHz)")
				plt.ylabel("Normalized Amplitude")

				legend = plt.legend(loc='best', shadow=False, fontsize='x-large')

				plt.show()
				'''
				'''
				plt.plot(np.divide(vch4[830:910],np.amax(vch4[830:900])),'b-',label="Raw")
				plt.plot(np.divide(deco[830:910],np.amax(deco[830:900])),'r-',label="Deconvolved")
				plt.ylabel("Normalized Amplitude")
				plt.xlabel("Time (16ns ticks)")
				legend = plt.legend(loc='best', shadow=False, fontsize='x-large')
				plt.show()
				'''
				rolltime=int(t4)-10
				vch4[0:rolltime]=0
				vch4=np.roll(vch4,-1*int(rolltime))
				deco[0:rolltime]=0
				deco=np.roll(deco,-1*int(rolltime))
				for x in range(0,len(deco)):
					if deco[x]<0:
						deco[x]=0
				#plt.plot(vch4)
				#plt.plot(deco)
				#plt.show()
				if ((amp4>80.0)&(amp4<1000)):
					sum_wf=np.add(sum_wf,vch4)
					sum_deco=np.add(sum_deco,deco)
					num_summed+=1
print "NUM:",num_summed
plt.plot(np.arange(0,2000),sum_wf)
plt.show()
plt.plot(sum_deco)
plt.show()

outfile=ROOT.TFile("time_density.root","RECREATE")
density=ROOT.TH1D("density","density",2000,0,2000)
density_deco=ROOT.TH1D("density_deco","density_deco",2000,0,2000)

for x in range(0,2000):
	density.SetBinContent(x,sum_wf[x])
	density_deco.SetBinContent(x,sum_deco[x])

density.Write()
density_deco.Write()

outfile.Close()

np.savez("numpy_density.npy",sum_wf,sum_deco)


#sum_wf=np.zeros(2000)
#
#
#for event in range(0,len(co_events)):
#	entry=5*int(co_events[event])
#	t.GetEntry(entry)
#
#	vch4[0:t4[event]]=0
#	vch4=np.roll(vch4,-1*int(t4[event]))
#	print t4[event]
#
#	tch4=np.arange(0,len(vch4))	
#	t.GetEntry(entry+1)
#	vch13=np.asarray(t.wf)
#	vch13=np.add(vch13,-1.0*ped13[event])
#	vch13[0:t13[event]]=0
#	vch13=np.roll(vch13,-1*int(t13[event]))
#
#
#	sum_wf=np.add(sum_wf,vch4)
#	tch13=np.arange(0,len(vch13))	
#	#plt.cla()
#	#plt.plot(tch4,vch4,'r-')
#	#plt.plot(tch13,vch13,'b-')
#
#	#plt.show()
#
#
#plt.plot(tch4,sum_wf)
#plt.show()
#
##sum_wf=np.divide(sum_wf,np.amax(sum_wf))
#
#outfile= ROOT.TFile("channel_4_tconst.root","RECREATE")
#
#h4=ROOT.TH1D("h4","h4",2000,0,2000)
#
#for x in range(0,2000):
#	h4.SetBinContent(x+1,sum_wf[x])
#
#h4.Write()
#outfile.Write()
#outfile.Close()
#
#print "NEV: ",len(co_events)
#
#'''
#
#	while event
#	print j
#	tp.GetEntry(j)
#	print tp.ch
#	if (tp.ch==13):
#		print "LAMB"
#
#
#for x in xrange(0,t.GetEntries()):
#	t.GetEntry(x)
#	print x
#	print t.ch
#	print t.event
#	print
#	starts=[]
#	ends=[]
#
#	voltage=np.asarray(t.wf)
#	samples=np.arange(len(t.wf))
#	plt.plot(samples,voltage)
#	plt.show()
#
#'''
#
#
#'''
#	while tp.event==t.event:
#		if (t.ch==13):
#			starts.append(tp.tstart)
#			ends.append(tp.tend)
#			voltage=np.asarray(t.wf)
#			samples=np.arange(len(t.wf))
#			plt.plot(samples,voltage)
#			threshold=np.ones(len(t.wf))
#			pedmean=np.ones(len(t.wf))
#			threshold=np.multiply(threshold,tp.thresh)
#			#print "thresh: ",tp.thresh
#			pedmean=np.multiply(pedmean,tp.ped)
#			#print "mean: ",tp.ped
#			plt.plot(samples,threshold,'g-')
#			plt.plot(samples,pedmean,'r-')
#			#print starts
#			#print ends 
#			#print "LENGTH: ", len(voltage)
#			#print "STARTLEN: ", len(starts)
#			#print "ENDLEN: ", len(ends)
#			#print voltage
#			if len(voltage)<2000:
#				continue
#			#plt.plot(starts,voltage[starts],'ro')
#			#plt.plot(ends,voltage[ends],'go')
#				#if(tp.reject_edge==0):
#				#print "AREA: ",tp.area
#			plt.show()
#			plt.cla()
#		y+=1
#		tp.GetEvent(y)
#'''
#'''
#	if raw_input('')=='0':
#		#print "ADDED"
#		sum_rms+=np.std(voltage)
#		nadded+=1
#		#print sum_rms
#		#print nadded
#'''
#
##print sum_rms
##print nadded
#'''
#for entry in t:
#	voltage=np.asarray(entry.wf)
#	samples=np.arange(len(entry.wf))
#	plt.plot(samples,voltage)
#	plt.show()
#	#print t.GetEntry(0)
#	raw_input('')
#'''
