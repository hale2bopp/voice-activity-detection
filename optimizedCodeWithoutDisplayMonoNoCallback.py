import pyaudio
import numpy as np
from scipy.fftpack import fft, dct,idct
from scipy import signal
from scipy.signal import argrelmax
from scipy.stats.mstats import gmean
import heapq
import struct



reqdRes = 30 #Hz

fs = 44100
RATE = 44100
CHANNELS = 1
WIDTH = 2
fftLen = 2**(np.ceil(np.log2(fs/float(reqdRes))))
CHUNK = 4096#np.int32(fftLen)+1000

# check for 
p = pyaudio.PyAudio()

a = p.get_device_count()

for i in range (0,a):
	if p.get_device_info_by_host_api_device_index(0,i).get('maxInputChannels')==CHANNELS:
	        print "Input Device id ", i, " - ", p.get_device_info_by_host_api_device_index(0,i).get('name')
		deviceNum = i		
		 
devinfo = p.get_device_info_by_index(deviceNum)
print "Selected device is ",devinfo.get('name')
# List of available audio devices


stream = p.open(format=p.get_format_from_width(WIDTH),
	        channels=CHANNELS,
	        rate=RATE,
	        input=True,
	        output=True,
	        input_device_index=deviceNum,
	        frames_per_buffer=CHUNK
		)


# Find peaks in the spectrogram/cepstrogram after smoothing each time frame data
def find_peaks(cepstr):
	value = []
	N  = 3    # Filter order
	Wn = 0.1 # Cutoff frequency
	B, A = signal.butter(N, Wn, output='ba')
	arr1 = []
	maxn_rows = []
	relmax_rows = []
	smooth_data = np.array(signal.filtfilt(B,A, cepstr))	# A
	relmax_rows = argrelmax(smooth_data)	        # B 
	values = np.array(smooth_data[relmax_rows])# C
	arr1 = heapq.nlargest(3,range(len(values)),values.take) # D			
	try:
		value.append(int(relmax_rows[0][arr1[1]]))
	except IndexError:
		value.append(0)
	return value, relmax_rows

print fftLen

# The arrays to plot frequency and time
full_freq = np.linspace(0,fs,fftLen)
R = []
LSFM_bufSize = 10	
LSFM_buf=[0]*LSFM_bufSize
threshold = -30

# check past 10 frames for fundamental frequency to check the variation in the measured fundamental
# if the fundamental varies a lot, another way to tell that it's noise

while(1):	
	data = stream.read(CHUNK)
	shorts = struct.unpack('h'*(CHUNK),data)
	input_wave = np.array(list(shorts),dtype = float)
# update the FFT buffer and give space for a new entry	
	LSFM_buf.pop(0)

# take FFT of newest frame and compute cepstrum	
	R = np.abs(dct(np.multiply(input_wave[0:fftLen], np.hamming(fftLen))))/np.sqrt(fftLen)	
	AM = np.mean(R)
	GM = gmean(R)
	LSFM_buf.append(100*np.log10(GM/float(AM)))
# Compute spectral flatness coefficient
	max_LSFM = max(LSFM_buf)
	min_LSFM = min(LSFM_buf)
	
	max_val,relmax_val= find_peaks(R)
	funda_freq = full_freq[max_val]
# print fundamental frequency
	if funda_freq > fs/2.0:
		funda_freq = fs-funda_freq
		print funda_freq ,' | ',  LSFM_buf[0]
	else:	
		print funda_freq,' | ',  LSFM_buf[0]
# decide whether noise, music or voice
	if max_LSFM<=threshold and funda_freq > 370:
		print 'music'
	elif max_LSFM>threshold:
		print 'noise'
	elif max_LSFM<=threshold and funda_freq <= 370:
		print 'voice'
	print ' '

stream.stop_stream()
stream.close()
p.terminate() 

