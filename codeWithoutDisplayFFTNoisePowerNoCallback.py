import pyaudio
import numpy as np
from scipy.fftpack import dct
from scipy import signal
from scipy.signal import argrelmax
from scipy.stats.mstats import gmean
import heapq
import struct

reqdRes = 60 #Hz
fs = 44100
RATE = 44100
CHANNELS = 2
WIDTH = 2
fftLen = 2**(np.ceil(np.log2(fs/float(reqdRes))))
CHUNK = 4096#np.int32(fftLen)+1000
# check for 
p = pyaudio.PyAudio()
stream = p.open(format=p.get_format_from_width(WIDTH),
	        channels=CHANNELS,
	        rate=RATE,
	        input=True,
	        output=True,
	        frames_per_buffer=CHUNK,
		)


# Find peaks in the spectrogram after smoothing each time frame data
def find_peaks(cepstr):
	value = []
	N  = 3    # Filter order
	Wn = 0.1 # Cutoff frequency
	B, A = signal.butter(N, Wn, output='ba')
	arr1 = []
	maxn_rows = []
	relmax_rows = []
	# Smooth data
	smooth_data = np.array(signal.filtfilt(B,A, cepstr))	# A
	# find indices of the relative maxima
	relmax_rows = argrelmax(smooth_data)	        # B 
	values = np.array(smooth_data[relmax_rows])# C
	# Find the 3 greatest relative maxima and choose the index of the second peak after zero
	arr1 = heapq.nlargest(3,range(len(values)),values.take) # D			
	try:
		value.append(int(relmax_rows[0][arr1[1]]))
	except IndexError:
		value.append(0)
	return value, relmax_rows

print fftLen

# The arrays to plot frequency and time
full_freq = np.linspace(0,fs/6,fftLen/2)
R = []
threshold = 10

# check past 10 frames for fundamental frequency to check the variation in the measured fundamental
# if the fundamental varies a lot, another way to tell that it's noise
background = 10000000
while(1):	
	data = stream.read(CHUNK)
	shorts = struct.unpack('h'*(CHUNK*2),data)
	input_wave = np.array(list(shorts),dtype = float)
# take FFT of newest frame 
	R = np.abs(dct(np.multiply(input_wave[0:fftLen], np.hamming(fftLen))))/np.sqrt(fftLen)
	R = R[0:fftLen/2]	
	pwr = np.sum(np.square(input_wave))
# find the power of the frame with minimum power to use it as a baseline	
	if pwr<background:
		background = pwr
	rel_pwr = pwr/float(background)
# find the peak frequency
	max_val,relmax_val= find_peaks(R)
	funda_freq = full_freq[max_val]
# decide whether noise, music or voice based on relative power and peak frequency
	if rel_pwr>threshold and funda_freq > 370:
		print 'music'
		print funda_freq
	elif rel_pwr<=threshold:
		print 'noise'
	elif rel_pwr>threshold and funda_freq <= 370:
		print 'voice'
		print funda_freq
	print ' '

stream.stop_stream()
stream.close()
p.terminate() 

