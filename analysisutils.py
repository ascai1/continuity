from numpy import empty, append, conj, sin, arange, pi, real, copy, zeros
from scipy import fft, ifft
from scipy.stats import binom

def get_subsample(sample, start, end, window):
	start = max(0,start)
	end = min(len(sample),end)

	if not window:
		return copy(sample[start:end,0])
	return triangle_window(sample[start:end,0])

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

def find_peaks(arr, delta):
	maxima = list()

	if not len(arr):
		return maxima

	loc_max = arr[0]
	loc_max_pos = 0

	loc_min = arr[0]
	min_next = False

	for i in xrange(1, len(arr)):
		if arr[i] > loc_max:
			loc_max = arr[i]
			loc_max_pos = i
		elif arr[i] < loc_min:
			loc_min = arr[i]
			loc_min_pos = i

		if not min_next and arr[i] < loc_max - delta:
			maxima.append((loc_max_pos, loc_max))
			loc_min = arr[i]
			loc_min_pos = i
			min_next = True
		elif min_next and arr[i] > loc_min + delta:
			loc_max = arr[i]
			loc_max_pos = i
			min_next = False

	return maxima

def triangle_window(arr):
	mid = len(arr) / 2.0
	result = copy(arr)
	for i in xrange(len(arr)):
		result[i] *= 1 - abs(i-mid)/mid
	return result

def get_padded_waveform(n, nz, period):
	w = empty(n)
	w[:nz] = sin(arange(nz)*2*pi/float(period))
	return w

