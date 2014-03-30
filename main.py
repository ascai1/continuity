#!/usr/bin/python
from findfreq import get_regions
from analysisutils import find_peaks
from notes import get_note
import scipy.io.wavfile as wav
import sys
from scipy import fft

def corr_main():
	(rate, sample) = wav.read(sys.argv[1])
	regions = get_regions(rate, sample,400)
	for start, region in regions:
		print start/float(rate), (start + len(region))/float(rate)

	for _, region in regions:
		print get_note(region, rate)[0]

if __name__ == '__main__':
	corr_main()
