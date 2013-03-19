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
    to use this class you must give it a window function. for example
    myOptFil = optimalFilter()
    myOptFil.window = tukeyWindow(512, 0.3).window

    also, you must give it the noise power
    myOptFil.noisepower = anoisepower  #anoisepower is a numpy array

    after you call setTemplate, you should call calcKernel just once before you 
    can use the class to estimate amplitudes of singals with the calcAmp method.

    setTemplate and calcAmp take numpy arrays as inputs

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


  def calcKernel(self):
    self.N = (2* (len(self.templatefft) - 1))
    
    self.templatePower = numpy.abs(self.templatefft)**2 / float(self.N)
    denom_intergrand = self.templatePower / self.noisepower
    self.denom = denom_intergrand.sum()

    self.kernel = (self.templatefft.conjugate()/self.noisepower) / self.denom



  def calcAmp(self, inputPulse):

    self.input_window = inputPulse * self.window

    self.inputfft = numpy.fft.rfft(self.input_window)

    self.amp_estimator_integrand = self.kernel*self.inputfft

    #!!! why doesn't the inverse fourier transform produce the correct results!
    self.amp_estimator = numpy.fft.irfft(self.amp_estimator_integrand)/numpy.sqrt(math.pi)

    # instead, the integration over the real parts does work. 
    self.amp_estimator2 = numpy.zeros(self.N)
    for time in range(self.N):
      for freq_k in range(len(self.inputfft)):
        self.amp_estimator2[time] += self.amp_estimator_integrand[freq_k].real * cmath.exp(1j*self._phase(freq_k, time)).real / self.N

  def _phase(self, k, time):
    '''
    returns 2 * pi * f_k * time, where f_k = k / N
    '''
    return 2.0*cmath.pi*float(k)*float(time)/float(self.N)


  def _chi2Element(self, amp, index, freq_k):

    #this is the numerator in the chi^2 integral
    diff = self.inputfft[freq_k] - float(amp)*cmath.exp(-1j*self._phase(freq_k, index))*self.templatefft[freq_k]
    
    return numpy.abs(diff)**2/self.noisepower[freq_k]/self.N

  # for testing, this was the explicitly written chi2 function
  # def _chi2Elementb(self, amp, index, freq_k, debug = False):

  #   sigre = float(self.inputfft[freq_k].real)
  #   sigim = float(self.inputfft[freq_k].imag)

  #   tempre = float(self.templatefft[freq_k].real)
  #   tempim = float(self.templatefft[freq_k].imag)
    
  #   sig2 = sigre*sigre + sigim*sigim
  #   temp2 = tempre*tempre + tempim*tempim
  #   phase = self._phase(freq_k, index)

  #   chi2 = (sig2 + float(amp) * float(amp) * temp2) / self.noisepower[freq_k] 
  #   chi2 -= 2.0 * amp * (sigre*tempre + sigim*tempim)*math.cos( phase ) / self.noisepower[freq_k] 
  #   chi2 += 2.0 * amp * (tempim*sigre - sigim*tempre)*math.sin( phase ) / self.noisepower[freq_k] 

  #   if debug:
  #     print 'amp, index, freq_k, sig_re, sig_im, temp_re, temp_im, sig**2, temp**2, phase, cos(phase), sin(phase), noisepower[k]'
  #     print amp, index, freq_k, sigre, sigim, tempre, tempim, sig2, temp2, phase, math.cos(phase), math.sin(phase), self.noisepower[freq_k]
  #     print 'chi2 line 1', (sig2 + float(amp) * float(amp) * temp2) / self.noisepower[freq_k] 
  #     print 'chi2 line 2', -2.0 * amp * (sigre*tempre + sigim*tempim)*math.cos( phase ) / self.noisepower[freq_k] 
  #     print 'chi2 line 3', 2.0 * amp * (tempim*sigre - sigim*tempre)*math.sin( phase ) / self.noisepower[freq_k] 

  #   return chi2

  def _chi2(self, amp, index):
    '''
      calculate chi**2 for a particular pulse amplitude, amp, and at a particular sample point, index.
    '''
    
    chi2 = numpy.array([self._chi2Element(amp, index, k) for k in range(len(self.noisepower))])

    return chi2.sum()/self.N

  def chi2_time(self, amp = None, chi2func = None):  
    if chi2func is None:
      chi2func = self._chi2
    if amp is None:
      amp = self.amp_estimator.max()
      if amp < numpy.abs(self.amp_estimator.min()):
        amp = self.amp_estimator.min()

    return numpy.array( [ chi2func(i, amp) for i in range(self.N)] )

  def chi2_amp(self, time = None, chi2func = None, ampMin = -1000, ampMax = 1000 numSteps = 2000):
    
    if chi2func is None:
      chi2func = self._chi2
    if time is None:
      time = numpy.abs(self.amp_estimator).argmax()
    
    stepSize = ampMax - ampMin/numSteps
   
    return numpy.array([ chi2func(time, ampMin + i * stepSize) for i in range(numSteps) ])


  def variance(self):
    return self.denom * self.N
