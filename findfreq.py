#!/usr/bin/python
import scipy.io.wavfile as wav
from numpy import empty, append, conj, sin, arange, pi, real, array, zeros
from scipy import fft, ifft
from scipy.stats import binom
from pylab import plot, subplot, show, xlim, ylim, axis, autoscale
import sys
import operator

def zero_padded(arr, start, end):
	for i in xrange(start, end):
		if i < 0 or i >= len(arr):
			yield 0
		else:
			yield arr[i]

def smooth(arr, n):
	avgs = empty(len(arr))
	pm = binom.pmf(xrange(n*2),n*2,0.5)
	for i in xrange(len(avgs)):
		subarr = zero_padded(arr,i-n,i+n)
		avgs[i] = sum(map(operator.mul, pm, subarr))
	return avgs

def correlation(left, right):
	# pad both left and right to the same length
	cor_len = max(len(left), len(right))
	return ifft(fft(left, cor_len) * fft(conj(right), cor_len))

def show_sig(rtransform, start, end, rate, subsample):
	interval = end - start
	subsample = subsample[:interval]
	intpersec = rate/interval
	sig = filter(lambda x:x>1, rtransform)
	if not len(sig):
		return
	print start/float(rate)
	print sig
	subplot(3,1,1)
	axis([start, end, -0.5, 0.5])
	plot(xrange(start, end), subsample)
	subplot(3,1,2)
	xlim(0,rate)
	autoscale(axis='y')
	plot(xrange(0, rate), rtransform)
	subplot(3,1,3)
	xlim(0,rate)
	autoscale(axis='y')
	plot(xrange(0, rate), smooth(rtransform, 10))
	show()

def get_subsample(sample, start, end, window):
	start = max(0,start)
	end = min(len(sample),end)

	if not window:
		return array(sample[start:end,0])

	subsample = empty(end-start)
	mid = (start + end) / 2
	for i in xrange(start, end):
		window_factor = 1 - abs(i-mid)/float(mid-start)
		subsample[i-start] = (sample[i][0] * window_factor)
	return subsample

def main():
	(rate, sample) = wav.read(sys.argv[1])
	sample = sample/(2.0**15)
	intpersec = 50
	interval = rate/intpersec
	for i in xrange(len(sample)/interval):
		subsample = get_subsample(sample, i*interval, (i+1)*interval)
		transform = fft(subsample, rate)
		rtransform = abs(transform)

		# step 2: get autocorrelation
		# step 3: fix autocorrelation accounting for varying areas of correlation (given zero padding)
		# step 4: find peaks in autocorrelation
		# step 5: given peaks, perform FFTs
		# step 6: find peaks in FFTs
		# step 7: given FFT peaks, calculate actual peak frequencies
		show_sig(rtransform, i*interval, (i+1)*interval, rate, subsample)


# vector-like buffer that doubles itself when appending subsample data
class SampleBuffer:
	def __init__(self):
		self.flush()

	def append(self, subsample):
		if self.buf is None:
			self.buf = array(subsample)
		else:
			while self.size + len(subsample) > len(self.buf):
				self.buf = append(self.buf, zeros(len(self.buf)))
			self.buf[self.size:self.size+len(subsample)] = subsample
			self.size += len(subsample)

	def get_size(self):
		return self.size

	def get_arr(self):
		return self.buf[:self.size]

	def flush(self):
		self.size = 0
		self.buf = None


def get_waveform(n, nz, period):
	w = empty(n)
	w[:nz] = sin(arange(nz)*2*pi/float(period))
	return w

def corr_main():
	(rate, sample) = wav.read(sys.argv[1])
	sample = sample/(2.0**15)
	intpersec = 100
	interval = rate/intpersec

	subsample = None
	samplebuf = SampleBuffer()
	for i in xrange(len(sample)/interval):
		# save previous subsample
		if subsample is not None:
			samplebuf.append(subsample)
		# get new subsample, do not apply any windows
		subsample = get_subsample(sample, i*interval, (i+1)*interval, window=False)
		# if subsample is too quiet, move on
		if max(subsample) < 0.05:
			continue

		# is there enough buffer data to perform an adequate comparison?
		if samplebuf.get_size() > len(subsample) and samplebuf.get_arr() is not None:
			samplebufarr = samplebuf.get_arr()

			# get the correlation between current subsample and previous subsample data
			norm_correlation = array(real(correlation(samplebufarr, subsample)))
			# chop off data from the end when the subsample is wrapping around the previous data
			# TODO: verify that this is not wrapping the previous subsample data around current subsample
			norm_correlation.resize(len(norm_correlation) - len(subsample))

			# get the index and magnitude of maximum correlation
			(max_cor_i, max_cor) = max(enumerate(norm_correlation), key=lambda x: x[1])

			# at the index of max correlation, get the maximum possible value for the correlation
			# given each of the waveforms involved
			norm_factor = sum(b**2 for b in samplebufarr[max_cor_i:max_cor_i+len(subsample)])
			norm_factor = max(norm_factor, sum(c**2 for c in subsample))
			# normalize the correlation as a fraction of the maximum possible correlation
			norm_correlation = norm_correlation/norm_factor
			max_cor = max_cor/norm_factor

			subplot(3,1,1)
			plot(subsample)
			subplot(3,1,2)
			plot(samplebufarr)
			subplot(3,1,3)
			plot(norm_correlation)
			show()

			if max_cor < 0.3:
				# not enough correlation, flush previous data and start building new buffer
				samplebuf.flush()
				# TODO: this scheme allows for two unrelated waveforms to be sent to buffer,
				# if a flush just happened

if __name__ == '__main__':
	#main()
	corr_main()
