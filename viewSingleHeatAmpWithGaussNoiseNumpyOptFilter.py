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

# print len(_templatePulse), len(myOptimalFilter.templatefft)
# print len(myOptimalFilter.kernel), len(myOptimalFilter.noisepower), len(myOptimalFilter.templatefft)

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


  



