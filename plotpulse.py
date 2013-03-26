import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gSystem.Load("libkpta")
import KDataPy.util as kutil
import scipy.optimize

plt.ion()

def runAndPlotOptFilter(event, pulse, pta=None, **kwargs):  

  cham = kwargs['cham']
  ionpt = kwargs.get('ionpulsestarttime', 0)

  if pulse.GetChannelName() not in kwargs['chanlist']: return

  print pulse.GetChannelName()

  optkamper = cham.GetOptimalKamper(pulse)
  thewindow = optkamper.GetWindow()
  #print thewindow.GetWindowSize()
  #print thewindow.GetOutputPulseSize()
  preproc = optkamper.GetPreProcessor()
  #print preproc.GetInputPulseSize(), preproc.GetOutputPulseSize()
  #print pulse.GetPulseLength()

  if kwargs.get('runCalculation', True):      
    if pulse.GetChannelName().startswith('chal'):
      optkamper.SetIonPulseStartTime(ionpt);  
    optkamper.MakeKamp(pulse)

  #print thewindow.GetWindowSize()
  #print thewindow.GetOutputPulseSize()
  #print preproc.GetInputPulseSize(), preproc.GetOutputPulseSize()
  #print pulse.GetPulseLength()

  resultsMap = optkamper.GetResults()
  for key, val in resultsMap:
    print key, val.fValue


  rawpulse = np.array(pulse.GetTrace())
  thewindow = optkamper.GetWindow()
  windowpulse = kutil.get_as_nparray(thewindow.GetWindow(), thewindow.GetWindowSize() )
  windowoutput = kutil.get_out( thewindow )

  optfilter = optkamper.GetOptimalFilter()

  noisepower = kutil.get_as_nparray(optfilter.GetNoiseSpectrum(), optfilter.GetNoiseSpectrumSize())
  templatepower = kutil.get_as_nparray(optfilter.GetTemplatePower(), optfilter.GetTemplatePowerSize())

  hc2p = ROOT.KHalfComplexPower()
  hc2p.SetInputPulse( optfilter.GetOptimalFilter(), optfilter.GetOptimalFilterSize())
  hc2p.RunProcess()
  optfilpower = kutil.get_out(hc2p)

  hc2p.SetInputPulse(optfilter.GetInputPulse(), optfilter.GetInputPulseSize())
  hc2p.RunProcess()
  thispulsepower = kutil.get_out( hc2p)

  hc2p.SetInputPulse( optfilter.GetOptFilterAndSignal(), optfilter.GetOptFilterAndSignalSize())
  hc2p.RunProcess()
  optfilter_andpulsepower = kutil.get_out(hc2p)

  ampestimator = kutil.get_out( optfilter)


  res = scipy.optimize.minimize(optkamper.Chi2Functor, [np.abs(windowpulse[ampestimator.argmax()]), ampestimator.argmax()], method='Nelder-Mead')
  chi2_min_amp = res.x[0]
  chi2_min_time = res.x[1]
  print res

  print 'chi2 start time after minimization', chi2_min_time
  print 'chi2 amp after minimization', chi2_min_amp
  print 'chi2 at mins', optfilter.GetChiSquared(chi2_min_amp, chi2_min_time)
  print 'amplitude theoretical variance', 1./optfilter.GetAmpEstimatorDenominator()

  print 'kdata: calcullate chi2 versus time at amplitude', chi2_min_amp
  stepsize = 0.01
  timearray = np.array([atime*stepsize for atime in range( int(pulse.GetPulseLength()/stepsize)) ])
  chi2_time = np.array([optfilter.GetChiSquared(chi2_min_amp, time) for time in timearray])
  
  print 'kdata: calcullate chi2 versus amplitude at time', chi2_min_time
  stepsize = 0.01
  amps = np.array([a*stepsize for a in range( 0, int(int(2*ampestimator.max())/stepsize)) ])
  chi2_amp = np.array([ optfilter.GetChiSquared(a_i, chi2_min_time) for a_i in amps ])

  print 'chi2 amp from scan', chi2_amp.argmin()*stepsize
  print 'chi2 start time from scan', chi2_time.argmin()

  


  theFig = plt.figure( kwargs.get('figure', 1))
  theWidth = 1

  plt.axis('off')

  axes  = plt.subplot(6,1,1)
  plt.cla()
  
  plt.plot(rawpulse, linewidth=theWidth)
  plt.plot(windowoutput, linewidth=theWidth)
  if kwargs.has_key('template'):
    scalefactor = -1*ampestimator.min()
    if ampestimator.max() > scalefactor:
      scalefactor = ampestimator.max()
    plt.plot( np.array(kwargs['template'])* scalefactor)

  axes = plt.subplot(6,1,2)

  plt.cla()
  plt.loglog(noisepower, linewidth=theWidth)
  plt.loglog(templatepower, linewidth=theWidth)
  plt.loglog(optfilpower, linewidth=theWidth)


  axes = plt.subplot(6,1,3)

  plt.cla()
  plt.loglog(thispulsepower, linewidth=theWidth)
  plt.loglog(optfilter_andpulsepower, linewidth=theWidth)


  axes = plt.subplot(6,1,4)

  plt.cla()
  plt.plot(ampestimator, linewidth=theWidth)


  axes = plt.subplot(6,1,5)

  plt.cla()
  plt.plot( amps, chi2_amp, linewidth=theWidth)


  axes = plt.subplot(6,1,6)

  plt.cla()
  plt.plot(timearray, chi2_time, linewidth=theWidth)

  plt.show()
  if kwargs.get('wait', True):
    raw_input()

