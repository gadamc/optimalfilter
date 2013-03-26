import signals
import rootpy
rootpy.log.basic_config_colorized()
from random import gauss
from rootpy.io import open as ropen
from rootpy.tree import Tree, TreeChain
from rootpy.plotting import Hist
import numpy as np
import random
import ROOT


#output ROOT file
outfilename = 'ampInAmpOutHeatGaussKDataFilter.root'
rootFileOut = ropen(outfilename, 'recreate')
tree = Tree("t")
branches = {
  'amp_in':'F', 'amp_out': 'F', 'amp_out_maxtime': 'F'
}
tree.create_branches(branches)

#signal/analyis characteristics
numEvents = 1000
_pulselength = 512
_pulseStartTime = _pulselength/2
_tukey_alpha = 0.3

amps = [-10, -20, -50, -100, -200, -400, -500, -800, -1000, -1400, -1800, -2000]
#amps = [-100]


#signal
#test different pulse shape
parameters = [514.08/2.016,-1, 21.04/2.016/2.4213 , 19.4/2.016 ,0.12, 0.88*129.54/2.016]
simpulse = signals.HeatSignal(_pulselength, parameters)
scaleFactor = -1.0*simpulse().min()


#noise
#gaussian white noise with windowing by tukey window
simnoise = signals.GaussianNoise(length=_pulselength, width=10)
noise_power_numpy = simnoise.width**2 * (723.55578301951266/900.0) * np.ones( len(simnoise)/2 + 1, dtype=np.float64 )  



#set up the optimal filter
channel_name = 'chal'
ROOT.gSystem.Load('libkamping')
ROOT.gSystem.Load('libkpta')
ROOT.gSystem.Load('libkds')
cham = ROOT.KChamonixKAmpSite()
cham.CreateHeatWindow(_pulselength, _tukey_alpha)

#  noise
noise_power_vp = ROOT.std.vector('double')(len(noise_power_numpy))
for i in range(len(noise_power_numpy)):
  noise_power_vp[i] =  noise_power_numpy[i]
cham.SetNoisePower('chal', noise_power_vp)

#  template
_template_Pulse_vp = ROOT.std.vector('double')(len(simpulse))
for i in range(len(simpulse)):
  _template_Pulse_vp[i] = simpulse[i]/scaleFactor
cham.SetTemplate('chal', _template_Pulse_vp, 0 ,0)


#make a KRawBoloPulseRecord
pRaw = ROOT.KRawBoloPulseRecord()
pRaw.SetPulseLength(simpulse.length)
pRaw.SetPulseTimeWidth( 1e6 ) #in nanoseconds
pRaw.SetIsHeatPulse(True)
pRaw.SetPretriggerSize(simpulse.length/2)
pRaw.SetChannelName('chal')

trace = ROOT.std.vector('short')(len(simpulse))

for anAmp in amps:
  print 'running', anAmp
  counter  = 0
  simpulse.setpar(1, anAmp/scaleFactor)  #this scale factor is here to force the amplitude to be equal to simsignalAmp

  while True:

    #generate a trace - but 
    randomStart = _pulseStartTime + random.uniform(-10, 80) #smear the start time
    simpulse.setpar(0, randomStart )

    simsignal = simnoise + simpulse
    
    for i in range(len(simsignal)):
      trace[i] =  int(simsignal[i])

    pRaw.SetTrace(trace)
    optkamper = cham.GetOptimalKamper(pRaw)
    optkamper.SetIonPulseStartTime(0.0);  
    optkamper.SetUseMinimizer(False)  #don't look for the minimum of chi**2

    optkamper.MakeKamp(pRaw)
    resultsMap = optkamper.GetResults()
      

    tree.amp_in = anAmp
    tree.amp_out = resultsMap['amp'].fValue
    tree.amp_out_maxtime = resultsMap['peakPosition'].fValue
    tree.Fill()

    counter += 1
    if counter == numEvents: 
      break

tree.write()
rootFileOut.close()

