#!/usr/bin/env python
import sys, os
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import ultracamutils
import matplotlib.pyplot
import argparse, subprocess
import numpy, math
import classes
import ultraspecClasses
from trm import ultracam
from trm.ultracam.UErrors import PowerOnOffError, UendError, UltracamError
import ultracam_shift
import time, datetime
import json
import Image
import ucamObjectClass
from photutils import datasets
from photutils import daofind
from photutils import aperture_photometry, CircularAperture, psf_photometry, GaussianPSF
import astropy.table, astropy.io
from astropy.stats import median_absolute_deviation as mad


if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Uses the Astrometry.net API to try to solve the field')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-p', '--preview', action='store_true', help='Show image previews with Matplotlib')
	parser.add_argument('-s', '--stack', action='store_true', help='Stack the images in the preview window')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('--startframe', default=1, type=int, help='Start frame. \'1\' is the default')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames. No parameter means all frames, or from startframe to the end of the run')
	parser.add_argument('-t', '--sleep', default=0, type=int, help='Sleep time (in seconds) between frames. \'0\' is the default')
	parser.add_argument('--xyls', action='store_true', help='Write an XYLS (FITS) file output catalog that can be used as input to Astronomy.net')
	
	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	
	runInfo = ultraspecClasses.runInfo(arg.runname)
	runInfo.loadFromJSON(config.RUNINFO)
	
	debug.write(runInfo, 2)
	
	coordinates = (runInfo.ra, runInfo.dec)
	debug.write(coordinates, 2)
	debug.write(ultracamutils.toSexagesimal(coordinates), 2)
		
	runDate, runID = ultracamutils.separateRunNameAndDate(arg.runname)
	
	""" Check that the stacked image is available
	"""
	(runDate, runNumber) = ultracamutils.separateRunNameAndDate(arg.runname)
	
	workingFolder = ultracamutils.addPaths(config.WORKINGDIR, runDate)
	
	stackedImageFilename = workingFolder + '/' + runNumber + ".png"
	if os.path.isfile(stackedImageFilename):
		print "Found - ", stackedImageFilename 
	else: 
		print "No stacked image found at:", stackedImageFilename
		
	solutionOutputFile = workingFolder + '/' + runNumber + "_wcs_solution.fits"
	FITSSolutionOutputFile = workingFolder + '/' + runNumber + "_wcs_solved_image.fits"
	# Run astrometryClient
	astrometryCommand = ['astrometryClient.py']
	astrometryCommand.append("-kpadlqljoevlogqik")
	astrometryCommand.append("-u" + stackedImageFilename)
	astrometryCommand.append("--wcs=" + solutionOutputFile)
	astrometryCommand.append("--wcsfits=" + FITSSolutionOutputFile)
	astrometryCommand.append("-w")
	astrometryCommand.append("--ra=" + str(runInfo.ra))
	astrometryCommand.append("--dec=" + str(runInfo.dec))
	astrometryCommand.append("--radius=1")
	debug.write(astrometryCommand, level =2)
	
	result = subprocess.call(astrometryCommand)
	
	print "Result:", result
	
