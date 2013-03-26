import signals
import optfilter
import rootpy
rootpy.log.basic_config_colorized()
from random import gauss
from rootpy.io import open as ropen
from rootpy.tree import Tree, TreeChain
from rootpy.plotting import Hist
import numpy as np
import random

outfilename = 'ampInAmpOutHeatArbitraryNumpyFilter.root'
rootFileOut = ropen(outfilename, 'recreate')
tree = Tree("t")
branches = {
  'amp_in':'F', 'amp_out': 'F', 'amp_out_maxtime': 'F'
}
tree.create_branches(branches)


numEvents = 1000
_pulselength = 512
_pulseStartTime = _pulselength/2
_tukey_alpha = 0.3

amps = [-10, -20, -50, -100, -200, -400, -500, -800, -1000, -1400, -1800, -2000]
#amps = [-100]

#signal
simpulse = signals.HeatSignal(length = _pulselength)
scaleFactor = -1.0*simpulse().min()


#noise
noise_power_numpy = 5000*np.genfromtxt('chalnoise.txt', unpack=True)  
simnoise = signals.ArbitraryNoise(noise_power_numpy)


#set up the optimal filter
myOptimalFilter = optfilter.optimalFilter()
myOptimalFilter.window = optfilter.tukeyWindow(_pulselength, _tukey_alpha).window
myOptimalFilter.noisepower = noise_power_numpy
tempPulse = simpulse / scaleFactor
myOptimalFilter.setTemplate(tempPulse)
myOptimalFilter.calcKernel()


for anAmp in amps:
  print 'running', anAmp
  counter  = 0
  simpulse.setpar(1, anAmp/scaleFactor)  #this scale factor is here to force the amplitude to be equal to simsignalAmp

  while True:

    #generate a trace - but 
    randomStart = _pulseStartTime + random.uniform(-10, 80) #smear the start time
    simpulse.setpar(0, randomStart )

    simnoise.generate()

    simsignal = simnoise + simpulse

    myOptimalFilter.calcAmp(simsignal)
    tree.amp_in = anAmp
    tree.amp_out = myOptimalFilter.amp_estimator.max()
    tree.amp_out_maxtime = myOptimalFilter.amp_estimator.argmax()
    tree.Fill()

    counter += 1
    if counter == numEvents: 
      break

tree.write()
rootFileOut.close()

