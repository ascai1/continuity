from scipy import fft
from numpy import log
from pylab import subplot, plot, log, copy, show
from analysisutils import find_peaks

# Get notes from regions passed in from get_regions.
# Currently this only supports detecting a single frequency.
# We need to experiment or look up how to match regions with notes.
def get_note(region, sample_rate):

	# Piano note.
	a4_pitch = 440

	# A list of frequencies can be found here:
	# http://en.wikipedia.org/wiki/Piano_key_frequencies

	# A sample of FFTs per piano key can be found here:
	# https://www.youtube.com/watch?v=5xjD6SRY8Pg

	# The perceived note that we hear and associate with a certain frequency
	# are actually a series of harmonic peaks. The fundamental frequency
	# might be missing from the FFT.

	# There are three ways that I know of to get a note through the FFT:
	# Look at peaks: I don't think this is happening without prior knowledge
	# to what the note will be (which we will have...)
	# Maybe we could look into some sort of machine learning for this.
	# Autocorrelation: Never looked into this.
	# Harmonic product spectrum: I'm using this for now.
	# Basic idea -- http://cnx.org/content/m11714/latest/

	# TBH I think this is the wrong approach, but I don't know what approach
	# to take given the notes we have from the MusicXML file.

	# Compress the FFT to 1/2, 1/3 & 1/4 to its size.
	# Could be made more efficient or something but whatever.
	max_harmonics = 4
	original_freqs = fft(region, sample_rate)
	hps_result = copy(original_freqs)

	for i in xrange(2, max_harmonics + 1):
		compressed_freqs = copy(original_freqs[::i])
		hps_result = hps_result[:len(compressed_freqs)]
		hps_result *= compressed_freqs

	# Find resulting peak here.
	return find_peaks(hps_result, 0.5)[0] # i dunno lol