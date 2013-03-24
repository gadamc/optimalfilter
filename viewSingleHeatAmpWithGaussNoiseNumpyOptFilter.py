import KDataPy.signals
import ROOT
import numpy as np
import plotpulse
import rootpy
rootpy.log.basic_config_colorized()
from random import gauss
from rootpy.io import open as ropen
from rootpy.tree import Tree, TreeChain
from rootpy.plotting import Hist
import optfilter
import matplotlib.pyplot as plt
plt.ion()

simnoise = KDataPy.signals.GaussianNoise(length=512, width=30)
simsignaltype = KDataPy.signals.HeatSignal
channel_name = 'chal'

simsignalAmp = -300

myOptimalFilter = optfilter.optimalFilter()
myOptimalFilter.window = optfilter.tukeyWindow(512, 0.0).window

#create the noise power data array and give it to the kampsite
noise_power_numpy = simnoise.width**2 * np.ones( len(simnoise)/2 + 1, dtype=np.float64 )
#noise_power_numpy = np.load('randomWindowedNoisePower.npy')  

myOptimalFilter.noisepower = noise_power_numpy


#give the kampsite a template
_templatePulse = simsignaltype()
_templatePulse /= -1.0*_templatePulse().min() #force the peak to be == -1

myOptimalFilter.setTemplate(_templatePulse)

myOptimalFilter.calcKernel()

print len(_templatePulse), len(myOptimalFilter.templatefft)
print len(myOptimalFilter.kernel), len(myOptimalFilter.noisepower), len(myOptimalFilter.templatefft)

simpulse = simsignaltype()

_neededForScale = -1.0
if simpulse().min() != 0:
  _neededForScale = simpulse().min()
simpulse.setpar(1, simsignalAmp/(-1.*_neededForScale))  #this scale factor is here to force the amplitude to be equal to simsignalAmp

print simsignalAmp, 'scaled to ', simsignalAmp/(-1.*_neededForScale)

while True:

  #if this is ArbNoise - generate a new noise instance
  if isinstance(simnoise, KDataPy.signals.ArbitraryNoise):
    simnoise.generate()

  #generate a trace
  simsignal = simpulse + simnoise

  myOptimalFilter.calcAmp(simsignal)

  fig1 = plt.figure(1)
  plt.subplot(6,1,1)
  plt.cla()
  plt.plot(simsignal)
  plt.plot(myOptimalFilter.input_window)
  plt.plot(_templatePulse * myOptimalFilter.amp_estimator.max())
  plt.plot(simpulse)

  plt.subplot(6,1,2)
  plt.cla()
  plt.loglog(myOptimalFilter.templatePower)
  plt.loglog(myOptimalFilter.noisepower)
  plt.loglog( np.abs(myOptimalFilter.kernel)**2/len(simsignal) )

  plt.subplot(6,1,3)
  plt.cla()
  plt.loglog( np.abs(myOptimalFilter.inputfft)**2/len(simsignal))
  plt.loglog( np.abs(myOptimalFilter.amp_estimator_integrand)**2/len(simsignal))

  plt.subplot(6,1,4)
  plt.cla()
  plt.plot( myOptimalFilter.amp_estimator)


  chi2Function = myOptimalFilter._chi2b

  chi2_time = np.zeros(myOptimalFilter.N)
  for time in range(myOptimalFilter.N):  #loop over time bins
    chi2_time[time] = chi2Function(time, 300)

  chi2_amp = np.zeros(600)
  for a_i in range(600):  #loop over amp bins
    chi2_amp[a_i] = chi2Function(myOptimalFilter.N/2, a_i)  #hack -- fixed time position... works for simulation 


  plt.subplot(6,1,5)
  plt.cla()
  #chi2 = myOptimalFilter.chi2(myOptimalFilter._chi2b)
  #chi2 = myOptimalFilter.chi2_amp()
  #plt.plot(chi2)
  #chi2b = myOptimalFilter.chi2_amp(chi2func = myOptimalFilter._chi2b)
  plt.plot(chi2_amp)

  plt.subplot(6,1,6)
  plt.cla()
  #chi2 = myOptimalFilter.chi2(myOptimalFilter._chi2b)
  #chi2 = myOptimalFilter.chi2_time(amp = np.abs(simsignalAmp), chi2func = myOptimalFilter._chi2b)
  #plt.plot(chi2/len(simsignal))
  #chi2b = myOptimalFilter.chi2_time(chi2func = myOptimalFilter._chi2b)
  plt.plot(chi2_time)
  

  print len(myOptimalFilter.amp_estimator), len(myOptimalFilter.amp_estimator_integrand)
  print myOptimalFilter.N
  print 'amplitude estimation max', myOptimalFilter.amp_estimator.max(), 'bin', myOptimalFilter.amp_estimator.argmax()
  print 'expected min chi2', myOptimalFilter._chi2( np.abs(simsignalAmp), len(simsignal)/2)
  print 'expected min chi2b', myOptimalFilter._chi2b( np.abs(simsignalAmp), len(simsignal)/2)
  print 'theoretical variance', myOptimalFilter.variance()
  raw_input()



