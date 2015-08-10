import ROOT
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy as sp
import scipy.interpolate as spint

f=ROOT.TFile("out.root")
t=f.Get("raw_wf_tree")
p=f.Get("pulse_tree")

sum_wf=np.zeros(100)
shortwf=np.zeros(100)
nsummed=0

ch4_amps=np.empty(0)
ch4_areas=np.empty(0)

amp_mean=3.82e01
amp_sigma=1.03e01
area_mean=2.23e02
area_sigma=9.30e01

cflag=True

for x in xrange(0,t.GetEntries()/5):
	if cflag==False:
		break
	t.GetEvent(x*5)
	p.GetEvent(x*5)
	if(p.npulses>0)&(p.reject_edge==False):
		myamps=np.asarray(p.amp)
		myareas=np.asarray(p.area)
	
		for y in xrange(0,len(myamps)):
			if cflag==False:
				break
			if ((myamps[y]>amp_mean-amp_sigma)&(myamps[y]<amp_mean+amp_sigma)):
				if ((myareas[y]>area_mean-area_sigma)&(myareas[y]<area_mean+area_sigma)):
					mywf=np.asarray(t.wf)
				



					mytstart=np.asarray(p.tstart)					
					mytend=np.asarray(p.tend)					
					print mytstart
					print mytstart[-1]
					overlap=False
					mytend=np.asarray(p.tend)					
					if (2000-mytend[-1]<50):
						overlap=True
					if (mytstart[0]<50):
						overlap=True

					if (y==len(mytstart)-1):
						pass
					elif(mytstart[y+1]-mytend[y])<=50:
						overlap=True	
					if (y==0):
						pass
					elif (mytstart[y]-mytend[y-1])<=50:
						overlap=True
					Start=np.zeros(1)
					End=np.zeros(1)
					Startwf=np.zeros(1)
					Endwf=np.zeros(1)
					Start[0]=mytstart[y]-50
					End[0]=mytend[y]+50

					if overlap==False:
						mywf=np.add(mywf,-1.0*p.ped)
						shortwf=mywf[Start[0]:Start[0]+100]
						shortwf=np.add(shortwf,-1.0*np.mean(shortwf[0:40]))
						sum_wf=np.add(sum_wf,shortwf)
						nsummed+=1
						if nsummed==100:
							cflag=False

print nsummed





sum_wf=np.divide(sum_wf,nsummed)
#plt.plot(xf,sum_wf)

matplotlib.rcParams.update({'font.size': 18})
plt.plot(np.arange(0,len(sum_wf)),sum_wf,'b-',label="Waveform")
plt.ylabel("Amplitude (ADC)")
plt.xlabel("Time (16ns ticks)")
plt.show()


np.save("spe_kernel.npy",sum_wf)



