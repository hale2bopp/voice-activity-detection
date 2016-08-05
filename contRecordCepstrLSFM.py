import pyaudio
import cv2
import numpy as np
from scipy.io import wavfile
from scipy.fftpack import fft, dct,idct
from scipy import signal
import matplotlib.pyplot as plt
import math
#from scipy.signal import argrelextrema
from scipy.signal import argrelmax
from scipy.stats.mstats import gmean
import heapq
import struct
from matplotlib import animation

t = 0 # time iterator
fs = 44100
CHUNK = 8192
RATE = 44100
CHANNELS = 2
WIDTH = 2

p = pyaudio.PyAudio()

a = p.get_device_count()

for i in range (0,a):
        if p.get_device_info_by_host_api_device_index(0,i).get('maxInputChannels')==2:
                print "Input Device id ", i, " - ", p.get_device_info_by_host_api_device_index(0,i).get('name')
		deviceNum = i		
		 
devinfo = p.get_device_info_by_index(deviceNum)
print "Selected device is ",devinfo.get('name')

stream = p.open(format=p.get_format_from_width(WIDTH),
                channels=CHANNELS,
                rate=RATE,
                input=True,
                output=True,
                input_device_index=deviceNum,
                frames_per_buffer=CHUNK)



# Just takes a double channel and makes it into a single channel
def make_single_channel(reqdLen,samp_rate, snd):
	if np.size(snd[0])==2:
		if len(snd)>(reqdLen*samp_rate)-1:
			input_wave = snd[0:(reqdLen*samp_rate)-1,0]
		else:
			input_wave = snd[:,0]
	elif np.size(snd[0])==1: 
		if len(snd)>(reqdLen*samp_rate)-1:
			input_wave = snd[0:(reqdLen*samp_rate)-1]
		else:
			input_wave = snd[:]
	return input_wave
# had an array with smoothed data : A
# formed an array with the indices of all local maxima in the original array : B
# formed an array with all the values of the local maxima from B : C
# formed an array with the indices of n largest local maxima from B. Indices correspond to C :D
# if i take the second value of D as 3, it means that the third local maxima was the second largest
# thus, i must take the fourth value of B to find the index of the second largest maximum in the original array A
# i.e., B[D[1]]
# return value is finally the index of the maximum value in that time frame


# Find peaks in the spectrogram/cepstrogram after smoothing each time frame data
def find_peaks(cepstr,start,end,c_or_s):
	rows = len(cepstr)
	cols = len(cepstr[0])
	value = []
	N  = 3    # Filter order
	Wn = 0.1 # Cutoff frequency
	B, A = signal.butter(N, Wn, output='ba')
	arr1 = [0 for i in range(rows)]
	maxn_rows = [[0 for j in range(1)] for i in range(rows)]
	relmax_rows = [[0 for j in range(1)] for i in range(rows)]
	for i in range(rows):
		smooth_data = signal.filtfilt(B,A, cepstr[i])	# A
		relmax_rows[i] = argrelmax(smooth_data)	        # B 
		values = [smooth_data[j] for j in relmax_rows[i]]# C
		values = np.array(list(values[0]))
#		print values
		arr1 = heapq.nlargest(3,range(len(values)),values.take) # D			
#		print arr1
		if c_or_s == 0: # if you are trying to find peak of cepstogram choose second greatest peak (ignore peak at 0)
			value.append(int(relmax_rows[i][0][arr1[1]]))
		elif c_or_s == 1: # if you want peak of spectogram then find the greatest peak
			value.append(int(relmax_rows[i][0][arr1[0]]))
#	print value
	return value[start:end], relmax_rows[start:end]

def divLists(a,b):
	if np.size(a) >1 and np.size(b) >1:
		return [i/float(j) for i,j in zip(a,b)]
	elif np.size(a) ==1 or np.size(b) ==1:
		return a/float(b)

def init():
    line.set_data([], [])
    return line

def geomet(arr):
	sum1 = 0
	inv = 1.0/len(arr)
	for j in arr:
		sum1+=np.log(j)	
#	print sum1
	mean = np.exp(sum1*inv)
	return mean






reqdRes = 20 #Hz
fftLen = 2**(np.ceil(np.log2(fs/float(reqdRes))))
print fftLen

# 1 millisecond = 1 second (fs samples) / 1000 
ms1 = 0*fs/1000
# 20 ms = 1 second (fs samples) / 50 
ms20 = 10*fs/1000

factor = 4.0
# The arrays to plot frequency and time
f = np.linspace(0,fs/(2*factor),fftLen/factor)
q = (1/float(fs))*np.arange(ms1,ms20)
full_freq = np.linspace(0,fs,fftLen)
rows = 500
cols = fftLen

bufSize = 1

R = [[0 for j in range(0,int(fftLen))] for i in range(0,bufSize)]
C = [[0 for j in range(0,int(fftLen))] for i in range(0,bufSize)]

fig = plt.figure()
ax = plt.axes(xlim=(0, 5), ylim=(-100, 100))
line, = ax.plot([], [], lw=2)

frame=0.0*np.ones((rows,cols,3));



x=[]
y=[]
initial_LSFM_buf = []
LSFM_bufSize = 10
LSFM_buf=[0 for i in range(LSFM_bufSize)]
ii=0
while(1):	
	data = stream.read(CHUNK)
	shorts = struct.unpack('h'*(CHUNK*2),data)
	samples = np.array(list(shorts),dtype = float)
#	print max(samples),min(samples)
#add noise
	noise = np.random.normal(-1,1,len(samples))
	input_wave = samples#+noise
#	SNR = 10*np.log(np.sum(np.abs(samples))/np.sum(np.abs(noise)))/np.log(10.0) #dB
#	print SNR
# update the FFT buffer and give space for a new entry	
	for i in range(1,bufSize):
		R[i][:] = R[i-1][:]

	for i in range(1,LSFM_bufSize):
		LSFM_buf[i] = LSFM_buf[i-1]

# take FFT of newest frame and compute cepstrum	
	R[0][:] = 5*np.log10(np.abs(dct(input_wave[0:fftLen])))/np.sqrt(fftLen)	
	print type(R),len(R[0])
	C[0][:] = np.log10(np.abs(dct(R[0][:])))/np.sqrt(882) 
	R = np.array(R)
	C = np.array(C)
	
	V_c = C[:,ms1:ms20]
	
#	print 'rows',len(V_c),'cols',len(V_c[0])
	# average spectral power over 10 frames
	S = 1/float(bufSize)*np.sum(np.square(np.abs(R)),0)	
#	print type(S),len(S)
# AM over frequency bins	
	AM = (1/float(fftLen))*np.sum(S)
#	print AM	
# GM over frequency bins	
	GM = gmean(S)
	print GM, geomet(S)
# Compute spectral flatness coefficient
	LSFM_buf[0] = 1000*np.sum(np.log(AM/GM))	
	max_LSFM = max(LSFM_buf)
	min_LSFM = min(LSFM_buf)
	threshold = (max_LSFM+min_LSFM)/2.0
	print 'max_LSFM',max_LSFM,'min_LSFM', min_LSFM,'threshold',threshold
	if ii<10:
		ii+=1
		initial_LSFM_buf.append(1000*np.sum(np.log(AM/GM)))	
	background = np.mean(initial_LSFM_buf)
#plot LSFM vs time for 8 seconds 
# 

#	x.append(t)
#	y.append(LSFM)
#	plt.axes(xlim=(0, t), ylim=(-100, 100))
#	plt.plot(x, y)
#	plt.pause(0.05)
#	line.set_data(x, y)
#	t+=0.1
#	ii+=1
#	if ii> RATE*8:
#		break

	max_val,relmax_val= find_peaks(R,0,1,1)
#	funda_freq = 1/float(q[max_val])
	funda_freq = full_freq[max_val]
	if funda_freq > fs/2.0:
		funda_freq = fs-funda_freq
		print funda_freq ,' | ',  LSFM_buf[0]
	else:	
		print funda_freq,' | ',  LSFM_buf[0]
	if LSFM_buf[0] > threshold:
		print 'noise'
	elif LSFM_buf[0] <= threshold and funda_freq > 200:
		print 'music'
	elif LSFM_buf[0] <= threshold and funda_freq <= 200:
		print 'voice'
	    #shift "frame" 1 up:

	frame[0:(rows-1),:]=frame[1:rows,:]; 

#Color mapping:
	#Red:
	frame[rows-1,:,2]=R[0]
	#Green:
	frame[rows-1,:,1]=np.abs(1-2*R[0])
	#Blue:
	frame[rows-1,:,0]=(1.0-R[0])
	#frame[rows-1,:,0]=frame[rows-1,:,1]**3
	# Display the resulting frame
	cv2.imshow('frame',frame)
	#Keep window open until key 'q' is pressed:
	if cv2.waitKey(1) & 0xFF == ord('q'):
		break
cv2.destroyAllWindows()


stream.stop_stream()
stream.close()
p.terminate() 

