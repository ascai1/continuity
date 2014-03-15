#!/usr/bin/python
from numpy import copy, real
from scipy import fft
from pylab import plot, subplot, show, xlim, ylim, axis, autoscale
import sys
import operator
from samplebuffer import SampleBuffer
from analysisutils import smooth, correlation, get_subsample, midsection
from math import sqrt


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

def get_regions(rate,sample,num_max_regions):
	output = list()

	sample = sample/(2.0**15)
	intpersec = 100
	interval = rate/intpersec

	subsample = None
	samplebuf = SampleBuffer(rate)

	for i in xrange(len(sample)/interval):
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
				# truncate the buffer to avoid calculating with the entire stored waveform
				trunc_buffer_waveform = midsection(buffer_waveform, len(current_waveform)*4)
				# get the correlation between current subsample and previous subsample data
				norm_correlation = copy(real(correlation(trunc_buffer_waveform, current_waveform)))
				# normalize by dividing by geometric mean of waveform energies
				# www-prima.imag.fr/jlc/papers/IAS95-martin.pdf
				norm_factor = sum(b**2 for b in trunc_buffer_waveform)
				norm_factor = sqrt(norm_factor*sum(c**2 for c in current_waveform))

				# chop off data from the end when the subsample is wrapping around the previous data
				# TODO: verify that this is not wrapping the previous subsample data around current subsample
				max_cor = max(norm_correlation[:-len(current_waveform)])/norm_factor

		end_period = subsample is None and buffer_waveform is not None
		end_period = end_period or max_cor is not None and max_cor < 0.35
		if end_period:
			if len(buffer_waveform) > len(current_waveform) * 4:
				output.append((i*interval-len(buffer_waveform),buffer_waveform))
				if num_max_regions is not None and num_max_regions == len(output):
					break
			samplebuf.flush()

		# save previous subsample
		if subsample is not None:
			samplebuf.append(subsample)

	return output
