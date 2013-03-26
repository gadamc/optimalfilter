import numpy
import math
import cmath

class tukeyWindow(object):

  def __init__(self, length, alpha):


    self.window = numpy.ones(length, dtype=numpy.float64)
    temp = alpha*(length-1.0)/2.0
    for i in range(int(temp)):
      self.window[i] = 0.5 + 0.5*math.cos( math.pi *( i/temp - 1.0) )
    
    for i in range( int(length - temp), length):
      self.window[i] = 0.5 + 0.5*math.cos( math.pi *( i/temp - 2.0/alpha - 1.0) )



class optimalFilter(object):
  '''
    to use this class you must give it a window function, a noise power spectrum,
    a template pulse and instruct your instance to calculate the optimal 
    filter kernel. 

    for example:
    
    myOptFil = optimalFilter()
    myOptFil.window = tukeyWindow(512, 0.3).window
    myOptFil.noisepower = anoisepower  #anoisepower is a numpy array
    myOptFil.setTemplate( aTemplate )
    myOptFil.calcKernel() #you should call calcKernel if you ever change the template or the noise_power

    #for each pulse,

    myOptFil.calcAmp(aPulse)

    ampEstimator = myOptFil.amp_estimator

    all arrays in this class should be numpy arrays.

  '''

  def setTemplate(self, templatePulse, makeShift = True):
    
    tempPulseWindow = templatePulse * self.window
    
    if makeShift:
      shift = len(tempPulseWindow)/2
      for i in range(len(tempPulseWindow)):
        if tempPulseWindow[i] != 0:
          shift = i
          break

      tempPulseWindow = numpy.concatenate( (tempPulseWindow[shift:] , numpy.zeros(shift, dtype=numpy.float64)) )

      #make sure that its an even number -- calculations may depend on this...
      if len(tempPulseWindow) % 2:
        tempPulseWindow = numpy.concatenate( ( tempPulseWindow, numpy.zeros(1, dtype=numpy.float64)) )

    
    self.templatefft = numpy.fft.rfft(tempPulseWindow)

    self.N = len(templatePulse)
    self.templatePower = numpy.abs(self.templatefft)**2 / float(self.N)


  def calcKernel(self):
    
    denom_intergrand = self.templatePower / self.noisepower

    self.denom = denom_intergrand.sum()

    self.kernel = (self.templatefft.conjugate()/self.noisepower) / self.denom



  def calcAmp(self, inputPulse):

    self.input_window = inputPulse * self.window

    self.inputfft = numpy.fft.rfft(self.input_window)

    self.amp_estimator_integrand = self.kernel*self.inputfft

    #inverse fourier transform
    self.amp_estimator = numpy.fft.irfft(self.amp_estimator_integrand)
    
    #...or explicit ingegral calculation.

    # self.amp_estimator = numpy.empty(self.N)
    # for time in range(self.N):
    #   self.amp_estimator[time] = -self.amp_estimator_integrand[0].real
    #   for freq_k in range(len(self.inputfft)):
    #     g = self.amp_estimator_integrand[freq_k] * cmath.exp(1j*self._phase(freq_k, time))
    #     self.amp_estimator[time] += (g.conjugate() + g).real
    # self.amp_estimator /= self.N

    # self.amp_estimator = numpy.empty(self.N)
    # print -1*len(self.inputfft)+1
    # for time in range(self.N):
    #   g = 0
    #   for freq_k in range(-1*len(self.inputfft)+1, len(self.inputfft)):
    #     if freq_k > 0:
    #       g += self.amp_estimator_integrand[freq_k] * cmath.exp(1j*self._phase(freq_k, time))
    #     else: 
    #       g += self.amp_estimator_integrand[numpy.abs(freq_k)].conjugate() * cmath.exp(1j*self._phase(freq_k, time)) 
    #   self.amp_estimator[time] = g.real
    # self.amp_estimator /= self.N


    #both, inverse fourier transform and the explicit integral calculation seem to give the same answer
    #I also observe the same answer in the C++/KData optimal filter.

    #HOWEVER - the amplitude estimation seems to be off by about ~sqrt(pi)????!!! and i can't figure out why. 

  def _phase(self, k, time):
    '''
    returns 2 * pi * f_k * time, where f_k = k / N
    '''
    return 2.0*cmath.pi*float(k)*float(time)/float(self.N)


  def _chi2Element(self, amp, index, freq_k):

    #this is the numerator in the chi^2 integral
    diff = self.inputfft[freq_k] - float(amp)*cmath.exp(-1j*self._phase(freq_k, index))*self.templatefft[freq_k]
    
    return 2.0*numpy.abs(diff)**2/self.noisepower[freq_k]

  # # for testing, this was the explicitly written chi2 function
  # def _chi2Element(self, amp, index, freq_k, debug = False):

  #   sigre = float(self.inputfft[freq_k].real)
  #   sigim = float(self.inputfft[freq_k].imag)

  #   tempre = float(self.templatefft[freq_k].real)
  #   tempim = float(self.templatefft[freq_k].imag)
    
  #   sig2 = sigre*sigre + sigim*sigim
  #   temp2 = tempre*tempre + tempim*tempim
  #   phase = self._phase(freq_k, index)

  #   chi2 = (sig2 + float(amp) * float(amp) * temp2) / self.noisepower[freq_k] 
  #   chi2 -= 2.0 * amp * (sigre*tempre + sigim*tempim)*math.cos( phase ) / self.noisepower[freq_k] 
  #   chi2 -= 2.0 * amp * (tempim*sigre - sigim*tempre)*math.sin( phase ) / self.noisepower[freq_k] 

  #   if debug:
  #     print 'amp, index, freq_k, sig_re, sig_im, temp_re, temp_im, sig**2, temp**2, phase, cos(phase), sin(phase), noisepower[k]'
  #     print amp, index, freq_k, sigre, sigim, tempre, tempim, sig2, temp2, phase, math.cos(phase), math.sin(phase), self.noisepower[freq_k]
  #     print 'chi2 line 1', (sig2 + float(amp) * float(amp) * temp2) / self.noisepower[freq_k] 
  #     print 'chi2 line 2', -2.0 * amp * (sigre*tempre + sigim*tempim)*math.cos( phase ) / self.noisepower[freq_k] 
  #     print 'chi2 line 3', -2.0 * amp * (tempim*sigre - sigim*tempre)*math.sin( phase ) / self.noisepower[freq_k] 

  #   return 2.0*chi2

  def chi2(self, amp, index):
    '''
      calculate chi**2 for a particular pulse amplitude, amp, and at a particular sample point, index.
    '''
    
    chi2 = numpy.array([self._chi2Element(amp, index, k) for k in range(len(self.noisepower))])

    return chi2.sum()/self.N

  def chi2functor(self, x):
    '''
      calculate chi**2 for a particular pulse amplitude, x[0], and at a particular sample point, x[1].
      This is just a wrapper around the chi2 method which can be used with SciPy's optimize/minimize tools.
    '''
    
    return self.chi2(x[0], x[1])

  

  def variance(self):
    return 1./(self.denom)
