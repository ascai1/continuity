#!/usr/bin/python
import scipy.io.wavfile as wav
from numpy import copy, real
from scipy import fft
from pylab import plot, subplot, show, xlim, ylim, axis, autoscale
import sys
import operator
from samplebuffer import SampleBuffer
from analysisutils import smooth, correlation, get_subsample


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


def corr_main():
	(rate, sample) = wav.read(sys.argv[1])
	sample = sample/(2.0**15)
	intpersec = 100
	interval = rate/intpersec

	subsample = None
	samplebuf = SampleBuffer(rate)
	for i in xrange(len(sample)/interval):
		# save previous subsample
		if subsample is not None:
			samplebuf.append(subsample)

		# get new subsample, do not apply any windows
		subsample = get_subsample(sample, i*interval, (i+1)*interval, window=False)
		buffer_waveform = samplebuf.get_arr()
		max_cor = None

		if max(subsample) < 0.05:
			subsample = None
		else:
			current_waveform = subsample
			if samplebuf.get_size() <= len(subsample):
				current_waveform = subsample[len(subsample)/2:]

			if buffer_waveform is not None:
				# get the correlation between current subsample and previous subsample data
				norm_correlation = copy(real(correlation(buffer_waveform, current_waveform)))
				# chop off data from the end when the subsample is wrapping around the previous data
				# TODO: verify that this is not wrapping the previous subsample data around current subsample
				norm_correlation.resize(len(norm_correlation) - len(current_waveform))

				# get the index and magnitude of maximum correlation
				(max_cor_i, max_cor) = max(enumerate(norm_correlation), key=lambda x: x[1])

				# at the index of max correlation, get the maximum possible value for the correlation
				# given each of the waveforms involved
				norm_factor = sum(b**2 for b in buffer_waveform[max_cor_i:max_cor_i+len(current_waveform)])
				norm_factor = max(norm_factor, sum(c**2 for c in current_waveform))
				# normalize the correlation as a fraction of the maximum possible correlation
				norm_correlation = norm_correlation/norm_factor
				max_cor = max_cor/norm_factor

		if subsample is None and buffer_waveform is not None:
			#samplebuf.print_info(i*interval, delta=5)
			samplebuf.flush()
		elif max_cor is not None and max_cor < 0.27:
			if len(buffer_waveform) > rate/12:
				samplebuf.print_info(i*interval, subsample, norm_correlation, delta=10)
				sys.stdout.flush()
				samplebuf.plot(i*interval, subsample, norm_correlation)
			samplebuf.flush()

if __name__ == '__main__':
	#main()
	corr_main()
