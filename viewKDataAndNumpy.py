import signals
import random
import ROOT
import numpy as np
import plotpulse
#import rootpy
#rootpy.log.basic_config_colorized()
#from random import gauss
# from rootpy.io import open as ropen
# from rootpy.tree import Tree, TreeChain
# from rootpy.plotting import Hist
import optfilter
import matplotlib.pyplot as plt
plt.ion()
import KDataPy.util as kut
import scipy.optimize

_pulselength = 512
_tukey_alpha = 0.3
simsignalAmp = -30
_pulseStartTime = _pulselength/2.0

#create the noise power data array and give it to the kampsite/filter

# simnoise = signals.GaussianNoise(length=_pulselength, width=10)
# #noise_power_numpy = simnoise.width**2 * np.ones( len(simnoise)/2 + 1, dtype=np.float64 )
# noise_power_numpy = simnoise.width**2 * (723.55578301951266/900.0) * np.ones( len(simnoise)/2 + 1, dtype=np.float64 )  #gaussian white noise with windowing by tukey window

noise_power_numpy = 1000*np.genfromtxt('chalnoise.txt', unpack=True)  
simnoise = signals.ArbitraryNoise(noise_power_numpy)


simsignaltype = signals.HeatSignal
#simsignaltype = signals.BBv2IonSignal

channel_name = 'chal'

myOptimalFilter = optfilter.optimalFilter()
myOptimalFilter.window = optfilter.tukeyWindow(_pulselength, _tukey_alpha).window

ROOT.gSystem.Load('libkamping')
ROOT.gSystem.Load('libkpta')
ROOT.gSystem.Load('libkds')
cham = ROOT.KChamonixKAmpSite()
cham.CreateHeatWindow(_pulselength, _tukey_alpha)





myOptimalFilter.noisepower = noise_power_numpy
noise_power_vp = ROOT.std.vector('double')(len(noise_power_numpy))
for i in range(len(noise_power_numpy)):
  noise_power_vp[i] =  noise_power_numpy[i]
cham.SetNoisePower('chal', noise_power_vp)



#give the kampsite a template
_templatePulse = simsignaltype(length = _pulselength)
_templatePulse /= -1.0*_templatePulse().min() #force the peak to be == -1
myOptimalFilter.setTemplate(_templatePulse)

_template_Pulse_vp = ROOT.std.vector('double')(len(_templatePulse))
for i in range(len(_templatePulse)):
  _template_Pulse_vp[i] = _templatePulse[i]
cham.SetTemplate('chal', _template_Pulse_vp, 0 ,0)


myOptimalFilter.calcKernel()

#print len(_templatePulse), len(myOptimalFilter.templatefft)
#print len(myOptimalFilter.kernel), len(myOptimalFilter.noisepower), len(myOptimalFilter.templatefft)

simpulse = simsignaltype(length = _pulselength)

_neededForScale = simpulse().min()
#print simsignalAmp, 'scaled to ', simsignalAmp/(-1.*_neededForScale)
simpulse.setpar(1, simsignalAmp/(-1.*_neededForScale))  #this scale factor is here to force the amplitude to be equal to simsignalAmp
#print 'simulated pulse amplitude: ', -1* simpulse().min()
# simpulse.setpar(0, 514.08/2.016 + 150)

#make a KRawBoloPulseRecord
pRaw = ROOT.KRawBoloPulseRecord()
pRaw.SetPulseLength(simpulse.length)
pRaw.SetPulseTimeWidth( 1e6 ) #in nanoseconds
pRaw.SetIsHeatPulse(True)
pRaw.SetPretriggerSize(simpulse.length/2)
pRaw.SetChannelName('chal')

r2hc = ROOT.KRealToHalfComplexDFT()
hc2r = ROOT.KHalfComplexToRealDFT()

while True:

  #if this is ArbNoise - generate a new noise instance
  if isinstance(simnoise, signals.ArbitraryNoise):
    simnoise.generate()

  #generate a trace
  randomAmp = simsignalAmp + random.gauss(0, 10)  #smear the amplitude a bit
  randomStart = _pulseStartTime + random.gauss(0, 40) 
  print 'pulse amp / start time: ', randomAmp, randomStart

  simpulse.setpar(0, randomStart )
  simpulse.setpar(1, randomAmp/(-1.*_neededForScale) )

  #randomAmp = simsignalAmp
  simsignal = simnoise + simpulse

  trace = ROOT.std.vector('short')(len(simsignal))
  for i in range(len(simsignal)):
    trace[i] =  int(simsignal[i])

  pRaw.SetTrace(trace)
  #print pRaw.GetPulseLength(), len(simsignal), trace.size()

  plotpulse.runAndPlotOptFilter('', pRaw, None, cham=cham, ionpulsestarttime = 0.003, chanlist = ['chal'], template = _template_Pulse_vp, wait=False,  tempAmp = np.abs(randomAmp) )


  myOptimalFilter.calcAmp(simsignal)


  fig1 = plt.figure(2)
  

  chi2Function = myOptimalFilter.chi2

  print 'numpy chi2 minimization'
  simpleAmp  = np.abs(myOptimalFilter.input_window[myOptimalFilter.amp_estimator.argmax()])
  res = scipy.optimize.minimize(myOptimalFilter.chi2functor, [simpleAmp, myOptimalFilter.amp_estimator.argmax()], method='Nelder-Mead')
  chi2_min_amp = res.x[0]
  chi2_min_time = res.x[1]
  print res

  #chi2_time = np.zeros(myOptimalFilter.N)
  print 'numpy: calcullate chi2 versus time at amplitude', chi2_min_amp
  chi2_time = np.array([ chi2Function(chi2_min_amp, time) for time in range(myOptimalFilter.N)])
  # for time in range(myOptimalFilter.N):  #loop over time bins
  #   print (str(time) + ' ')
  #   chi2_time[time] = chi2Function(np.abs(simsignalAmp), time)

  #chi2_amp = np.zeros(2*np.abs(simsignalAmp))
  print 'numpy: calcullate chi2 versus amplitude at time', chi2_min_time
  chi2_amp = np.array([ chi2Function(a_i, chi2_min_time) for a_i in range( int(-2*np.abs(randomAmp)), int(2*np.abs(randomAmp)) ) ])
  # for a_i in range(2*np.abs(simsignalAmp)):  #loop over amp bins
  #   print(str(a_i) + ' ')
  #   chi2_amp[a_i] = chi2Function(a_i, myOptimalFilter.N/2)  #hack -- fixed time position... works for simulation 

  bestFitPulse = _templatePulse
  bestFitPulse = simsignaltype(length = _pulselength)
  bestFitPulse.setpar(0, chi2_min_time)
  bestFitPulse.setpar(1, chi2_min_amp/_neededForScale)

  plt.subplot(6,1,1)
  plt.cla()
  plt.plot(simsignal)
  plt.plot(myOptimalFilter.input_window)
  plt.plot(simpulse)
  plt.plot(bestFitPulse)
  

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
  #plt.plot( myOptimalFilter.amp_estimator2)


  plt.subplot(6,1,5)
  plt.cla()
  #chi2 = myOptimalFilter.chi2(myOptimalFilter._chi2b)
  #chi2 = myOptimalFilter.chi2_amp()
  #plt.plot(chi2)
  #chi2b = myOptimalFilter.chi2_amp(chi2func = myOptimalFilter._chi2b)
  plt.plot( range( int(-2*np.abs(randomAmp)), int(2*np.abs(randomAmp)) ) , chi2_amp)

  plt.subplot(6,1,6)
  plt.cla()
  #chi2 = myOptimalFilter.chi2(myOptimalFilter._chi2b)
  #chi2 = myOptimalFilter.chi2_time(amp = np.abs(simsignalAmp), chi2func = myOptimalFilter._chi2b)
  #plt.plot(chi2/len(simsignal))
  #chi2b = myOptimalFilter.chi2_time(chi2func = myOptimalFilter._chi2b)
  plt.plot(chi2_time)
  

  #print len(myOptimalFilter.amp_estimator), len(myOptimalFilter.amp_estimator_integrand)
  #print myOptimalFilter.N
  print 'amplitude estimation max', chi2_min_amp, 'at time', chi2_min_time
  #print 'amplitude estimation2 max', myOptimalFilter.amp_estimator2.max(), 'bin', myOptimalFilter.amp_estimator2.argmax()
  print 'chi2', myOptimalFilter.chi2( chi2_min_amp, chi2_min_time)/(myOptimalFilter.N-2)
  print 'amplitude theoretical variance', myOptimalFilter.variance()


  chamOptKamper = cham.GetOptimalKamper(pRaw)

  #cross-correlation between the signal and the template -- for gaussian white noise this should
  #be exactly the same as the optimal filter... simply a matched filter
  #   inverse FourierTransform (template(f).conjugate() * signal(f))  ==  template(t) (X) signal(t)
  # ccfig = plt.figure(3)
  # plt.cla()
  # plt.plot(np.fft.irfft( myOptimalFilter.templatefft.conjugate() * np.fft.rfft(simsignal) / myOptimalFilter.templatePower.sum()) )


  # #just to be sure - test also with kdata tools
  # r2hc.SetInputPulse(trace)
  # r2hc.RunProcess()
  # trace_fft = kut.get_out(r2hc)  #in half-complex notation

  # trace_fft_vp = ROOT.std.vector('double')(len(trace_fft))
  # for i in range(len(trace_fft)):
  #   trace_fft_vp[i] = trace_fft[i]

  # template_kdatafft =  kut.get_as_nparray(chamOptKamper.GetOptimalFilter().GetTemplateDFT(), chamOptKamper.GetOptimalFilter().GetTemplateDFTSize()) 
  # template_kdatafft_conjugate = ROOT.std.vector('double')(len(template_kdatafft))
  # for i in range(len(template_kdatafft)/2 + 1):
  #   template_kdatafft_conjugate[i] = template_kdatafft[i]/myOptimalFilter.templatePower.sum()
  # j = len(template_kdatafft) - 1
  # for i in range(len(template_kdatafft)/2 + 1, len(template_kdatafft)):
  #   template_kdatafft_conjugate[i] = -1.0*template_kdatafft[i]/myOptimalFilter.templatePower.sum()
  #   j -= 1

  # product = ROOT.std.vector('double')(len(template_kdatafft))
  # complexAlgebra = ROOT.KHalfComplexArray();
  # complexAlgebra.Multiply(product, template_kdatafft_conjugate, trace_fft_vp)
  


  # hc2r.SetInputPulse(product)
  # hc2r.RunProcess()
  # plt.plot(kut.get_out(hc2r))



  #test for differences in the KChamoinixKAmpSite/OptimalKamper and the Numpy-based optimal camper
  #difffig = plt.figure(4)

  chamOptFilter = chamOptKamper.GetOptimalFilter()


  raw_input()


# print myOptimalFilter._phase(100, 256)

# print chamOptFilter.GetPhase(100, 256)
# print myOptimalFilter.N
# print chamOptFilter.GetPhase(1, 1)


# def getchi2diff(amp, time):
#   return np.array([ chamOptFilter.GetChiSquareElement(amp, time, k) - myOptimalFilter._chi2Elementb(amp, time, k) for k in range(chamOptFilter.GetNoiseSpectrumSize()) ])

# def printchi2diff(amp, time):
#   for k in range(chamOptFilter.GetNoiseSpectrumSize()):
#     print 'kdata'
#     chamOptFilter.GetChiSquareElement(amp, time, k, True)
#     print 'numpy'
#     myOptimalFilter._chi2Elementb(amp, time, k, True)
#     raw_input()


