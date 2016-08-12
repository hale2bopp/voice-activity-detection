import pyaudio
import numpy as np
from scipy.fftpack import dct
from scipy.signal import argrelmax
import heapq
import struct

# enter the required frequency resolution that you would like
reqdRes = 30 #Hz
fs = 44100
CHANNELS = 2
WIDTH = 2
# required resolution = fs/N, where N  is the size of the fft
# find the closest power of 2 greater than fs/Ns
fftLen = 2**(np.ceil(np.log2(fs/float(reqdRes))))
CHUNK = 4096

p = pyaudio.PyAudio()
stream = p.open(format=p.get_format_from_width(WIDTH),
	        channels=CHANNELS,
	        rate=fs,
	        input=True,
	        output=True,
	        frames_per_buffer=CHUNK,
		)

# Find peaks in the spectrogram after smoothing each time frame data
def find_peaks(cepstr):
	# Find the 3 greatest relative maxima and choose the index of the first peak
	arr1 = heapq.nlargest(3,range(len(cepstr)),cepstr.take) # D			
	# return the index of the second maximum		
	return int(arr1[1])

# The arrays to plot frequency and time
full_freq = np.linspace(0,fs,fftLen)
R = []
threshold = 10
freqs = np.fft.fftfreq(np.int16(fftLen), 1/float(fs))

# check past 10 frames for fundamental frequency to check the variation in the measured fundamental
# if the fundamental varies a lot, another way to tell that it's noise
background = 10000000
while(1):	
	data = stream.read(CHUNK)
	input_wave = np.fromstring(data, dtype=np.int16)
# take FFT of newest frame 
	R = np.abs(np.fft.fft(np.multiply(input_wave[0:fftLen], np.hamming(fftLen))))[0:(fftLen/2)]/np.sqrt(fftLen)
# find average power in the frame	
	pwr = (1/float(fftLen))*np.sum(np.square(input_wave[0:fftLen]))
# find the power of the frame with minimum power while recording (assumption - minimum power would be that of the background noise)	
	if pwr<background:
		background = pwr
# find the ratio of power of the current frame (pwr) to the minimum power recorded till now (background) 	
# in order to have a relative idea of how loud the input is. rel_pwr will always be >1 as if pwr is less than background (i.e. the power 
# of the current frame is less than the minimum), then the minimum power (background) will be updated.
	rel_pwr = pwr/float(background)
# find the peak frequency
	max_val = find_peaks(R)
	funda_freq = freqs[max_val]
# decide whether noise, music or voice based on relative power and peak frequency
# threshold is decided by checking the rel_pwr variable as some noise is brought near the microphone
	if rel_pwr>threshold and funda_freq > 370:
		print 'music'
		print funda_freq
	elif rel_pwr<=threshold:
		print 'noise'
		print funda_freq
	elif rel_pwr>threshold and funda_freq <= 370:
		print 'voice'
		print funda_freq
	print ' '

stream.stop_stream()
stream.close()
p.terminate() 

