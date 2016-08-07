# Convert RA and dec from HH:MM:SS and deg:min:sec to decimal degrees and vice versa

from __future__ import division, print_function
from astropy import units as u
from astropy.coordinates import SkyCoord

def hmsdms2decimal(hmsdms):
	'''
	Input should be a string in '14h12m00s +53d00m00s' format
	'''
	c = SkyCoord(hmsdms, frame='icrs')
	return(c.to_string('decimal'))

def decimal2hmsdms(ra, dec):
	'''
	Input should be two numbers -- RA and dec in decimal degrees
	'''
	c = SkyCoord(ra, dec, frame='icrs' ,unit='deg')
	return(c.to_string('hmsdms'))


# Examples:
if __name__ == "__main__":
	print(hmsdms2decimal('14h12m00s +53d00m00s'))
	print(decimal2hmsdms(213., 53.))
