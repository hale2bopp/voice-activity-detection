# X = feature vector with x0 =1 and x1 = power
# y : 1 (noise) or y=0 (voice/music) characterised by the power of the signal
# Theta is a vector with theta0 = 1 and thetha1 =?
# have an array with the noise
#
#|++++++++++++++++++++++++++++++++++++++++++++++++|
#| power 	|     voice/music(0) or noise (1) |     
#|++++++++++++++++++++++++++++++++++++++++++++++++|
#| 	5	|	1			  |
#| 	11	|	0			  |		
#|	12	|	0			  |	
#|	4	|	1			  |	
#|++++++++++++++++++++++++++++++++++++++++++++++++|
#
# Minimize J(0) h(0) as g(0TX)
# m : number of training examples
# n: dimensionality of the problem.
# y is a vector with zeros or ones.  m X 1 dimensional vector.
# X is the power vector, also n+1 x 1 dimensional.
# must find theta, in order to find the hypothesis h(0Tx)
# theta found by running gradient descent and minimizing cost function
# theta is a n+1x1 vector, theta0  and theta1 to be determined. give it some initial values and just see what happens	
# cost = (1/2*m)*np.sum( (g(theta,x)[i]-y[i])**2 )

import pyaudio
import numpy as np
from scipy.fftpack import dct
from scipy import signal
from scipy.signal import argrelmax
from scipy.stats.mstats import gmean
import heapq
import struct
import matplotlib.pyplot as plt
reqdRes = 60 #Hz
fs = 44100
RATE = 44100
CHANNELS = 2
WIDTH = 2
fftLen = np.int32(2**(np.ceil(np.log2(fs/float(reqdRes)))))
CHUNK = 4096#np.int32(fftLen)+1000
record_time = 2 #s
num_samples = record_time*fs
num_frames = np.int32(num_samples/fftLen)

# check for 
p = pyaudio.PyAudio()
stream = p.open(format=p.get_format_from_width(WIDTH),
	        channels=CHANNELS,
	        rate=RATE,
	        input=True,
	        output=True,
	        frames_per_buffer=CHUNK,
		)
def g(theta,x):
	return 1/float(1+np.exp(-1*np.dot(theta,x)))
	
def cost(x,y,theta):
	m=len(x[0])
	 
	cost_1 = 0
	print ' '
	print m
	for i in range(m):	
		print 'yalla',( y[i]*np.log(abs(g(theta,x[:,i])))+  (1-y[i])*np.log(abs(1-g(theta,x[:,i]))) )	
		cost_1+= (-1/m)*( y[i]*np.log(abs(g(theta,x[:,i])))+  (1-y[i])*np.log(abs(1-g(theta,x[:,i]))) )
		print 'cost_1',cost_1
	return cost_1

def grad_desc(theta,x,y,alpha):
	m=len(x)
	J = cost(x,y,theta)
	print J
	plt.plot(J)
	temp0=0
	temp1=0
	for i in range(m):
		temp0 += g(theta,x[:,i]) *x[0,i] 	
		temp1 += g(theta,x[:,i]) *x[1,i] 
	temp0 = theta[0] - temp0*(alpha/float(m))
	temp1 = theta[1] - temp1*(alpha/float(m))
	theta[0]=temp0
	theta[1]=temp1
	return theta,J

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

# The arrays to plot frequency and time
full_freq = np.linspace(0,fs/6,fftLen/2)
R = []
threshold = 10
frame_marker = np.array([0 for i in range(num_frames)])
rel_pwr=np.array([0 for i in range(num_frames)])

# check past 10 frames for fundamental frequency to check the variation in the measured fundamental
# if the fundamental varies a lot, another way to tell that it's noise
background = 10000000
maxIter=100
counter =0

total_wave = []
for j in range(num_frames):
	data = stream.read(CHUNK)
	shorts = struct.unpack('h'*(CHUNK*2),data)
	input_wave = np.array(list(shorts),dtype = float)
	pwr=np.sum(np.square(input_wave))
# find the power of the frame with minimum power to use it as a baseline	
	if pwr<background:
		background = pwr
	rel_pwr[counter] = pwr/float(background)
# decide whether noise, music or voice based on relative power and peak frequency
	if rel_pwr[counter]>threshold: 
		frame_marker[counter] = 0
		print 'music'
	elif rel_pwr[counter]<=threshold:
		frame_marker[counter] = 100
		print 'noise'
	print ' '
	counter+=1


x = np.array([np.ones((num_frames,1)),np.array(rel_pwr)])
theta = np.array([1.2,3])
J = np.zeros(maxIter)

for i in range(maxIter):	
	theta,J[i] = grad_desc(theta,x,frame_marker,0.005)

iterplot = np.arange(0,maxIter)	
plt.plot(iterplot,J)
plt.show()
while(1):
	data = stream.read(CHUNK)
	shorts = struct.unpack('h'*(CHUNK*2),data)
	input_wave = np.array(list(shorts),dtype = float)
	total_wave.append(input_wave[0:fftLen])
# take FFT of newest frame 
	R = np.abs(dct(np.multiply(input_wave[0:fftLen], np.hamming(fftLen))))/np.sqrt(fftLen)
	R = R[0:fftLen/2]	
	pwr=np.sum(np.square(input_wave))
# find the power of the frame with minimum power to use it as a baseline	
	if pwr<background:
		background = pwr
	rel_pwr = pwr/float(background)
# find the peak frequency
	max_val,relmax_val= find_peaks(R)
	funda_freq = full_freq[max_val]
# decide whether noise, music or voice based on relative power and peak frequency
	val = g(theta,np.array([1,rel_pwr]))
	if val<=0.5 and funda_freq > 370:		
		print 'music'
		print funda_freq
	elif val>0.5:
		print 'noise'
	elif val<=0.5 and funda_freq <= 370:
		print 'voice'
		print funda_freq
	print ' '

stream.stop_stream()
stream.close()
p.terminate() 

