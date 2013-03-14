import numpy
import math


class tukeyWindow(object):

  def __init__(self, length, alpha):


    self.window = numpy.ones(length, dtype=numpy.float64)
    temp = int(alpha*(length-1.0)/2.0);
    for i in range(temp):
      self.window[i] = 0.5 + 0.5*math.cos( math.pi *( float(i/float(temp)) - 1.0))
    
    for i in range(length - temp, length):
      self.window[i] = 0.5 + 0.5*math.cos( math.pi *( float(i/float(temp)) - 2.0/alpha - 1))



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

      tempPulseWindowShift = numpy.concatenate( (tempPulseWindow[shift:] , numpy.zeros(shift, dtype=numpy.float64)) )

      #make sure that its an even number
      if len(tempPulseWindowShift) % 2:
        tempPulseWindowShift = numpy.concatenate( ( tempPulseWindowShift, numpy.zeros(1, dtype=numpy.float64)) )

      self.templatefft = numpy.fft.rfft(tempPulseWindowShift)
    else:
      self.templatefft = numpy.fft.rfft(tempPulseWindow)


  def calcKernel(self):
    N = (2* (len(self.templatefft) - 1))
    
    self.templatePower = numpy.abs(self.templatefft)**2 / N
    denom_intergrand = self.templatePower / self.noisepower
    self.denom = denom_intergrand.sum()

    self.kernel = (self.templatefft/self.noisepower) / self.denom


  def calcAmp(self, inputPulse):

    self.input_window = inputPulse * self.window

    self.inputfft = numpy.fft.rfft(self.input_window)

    self.amp_estimator_integrand = self.kernel*self.inputfft 

    self.amp_estimator = numpy.fft.irfft(self.amp_estimator_integrand)

  def chi2(self, time):
    pass


