from numpy import array, arange, append, zeros
from scipy import fft
from pylab import plot, subplot, show
from analysisutils import find_peaks, triangle_window

# vector-like buffer that doubles itself when appending subsample data
class SampleBuffer:
	def __init__(self, rate):
		self.rate = rate
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
		if self.buf is None:
			return None
		return self.buf[:self.size]

	def flush(self):
		self.size = 0
		self.buf = None

	def print_info(self, now, subsample=None, norm_correlation=None, delta=1):
		if self.buf is None:
			return
		now = float(now)
		buffer_waveform = self.get_arr()
		print "({0}, {1})".format((now - len(buffer_waveform)) / self.rate, now / self.rate)
		rtransform = abs(fft(buffer_waveform))[:len(buffer_waveform)/2]
		peaks = find_peaks(rtransform, delta)
		peaks.sort(key=lambda x: -x[1])
		for peak in peaks:
			peak = (peak[0] * self.rate / float(len(buffer_waveform)), peak[1])
			print peak

	def plot(self, now, subsample=None, norm_correlation=None):
		if self.buf is None:
			return
		buffer_waveform = triangle_window(self.get_arr())
		if subsample is None or norm_correlation is None:
			subplot(2,1,1)
			x_axis = arange(now - len(buffer_waveform), now) / float(self.rate)
			plot(x_axis, buffer_waveform)
			subplot(2,1,2)
			x_axis = arange(0, self.rate, float(self.rate) / len(buffer_waveform))
			plot(x_axis, abs(fft(buffer_waveform)))
			show()
		else:
			subplot(4,1,1)
			plot(subsample)
			subplot(4,1,2)
			x_axis = arange(now - len(buffer_waveform), now) / float(self.rate)
			plot(x_axis, buffer_waveform)
			subplot(4,1,3)
			plot(norm_correlation)
			subplot(4,1,4)
			x_axis = arange(0, self.rate, float(self.rate) / len(buffer_waveform))
			plot(x_axis, abs(fft(buffer_waveform)))
			show()
