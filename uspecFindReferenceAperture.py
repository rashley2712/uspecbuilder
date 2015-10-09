#!/usr/bin/env python
import sys
"""try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
"""
import os
import ultracamutils
import matplotlib.pyplot
import argparse
import numpy, math
import classes
import ultraspecClasses
from trm import ultracam
from trm.ultracam.UErrors import PowerOnOffError, UendError, UltracamError
import ultracam_shift
import time, datetime
import json
from scipy import ndimage
from PIL import Image
import ucamObjectClass
from photutils import datasets
from photutils import daofind
from photutils import aperture_photometry, CircularAperture, psf_photometry, GaussianPSF
import astropy.table, astropy.io
from astropy.stats import median_absolute_deviation as mad
import astropy.stats.sigma_clipping
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel
from photutils.detection import detect_sources
from scipy.ndimage import binary_dilation
from photutils.background import Background
from photutils.morphology import (centroid_com, centroid_1dg, centroid_2dg)
from photutils import CircularAperture
from photutils import CircularAnnulus
import photutils
import ppgplot
import scipy.optimize

def shift_func(output_coords, xoffset, yoffset):
	return (output_coords[0] - yoffset, output_coords[1] - xoffset)

def gaussian(x, a, b, c, d):
	return a * numpy.exp(- (x - b)**2 / (2 * c**2)) + d
	
def shifting_gaussian(x, b):
	global a, c, d
	return a * numpy.exp(- (x - b)**2 / (2 * c**2)) + d
	

def determineFullFrameSize(windows):
	leftestPixel = 1057
	rightestPixel = 0
	topestPixel = 0 
	bottomestPixel = 1040
	for w in windows:
		if w.xll/w.xbin < leftestPixel: leftestPixel = w.xll/w.xbin
		if w.yll/w.ybin < bottomestPixel: bottomestPixel = w.yll/w.ybin
		if (w.xll/w.xbin + w.nx) > rightestPixel: rightestPixel = w.xll/w.xbin + w.nx
		if (w.yll/w.ybin + w.ny) > topestPixel: topestPixel = w.yll/w.ybin + w.ny
			
	return leftestPixel, bottomestPixel, rightestPixel, topestPixel
		

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Performs aperture photometry for Ultraspec runs. [dd-mm-yyyy/runxxx.dat]')
	parser.add_argument('runname', type=str, help='Ultraspec run name  [eg 2013-07-21/run010]')
	parser.add_argument('-p', '--preview', action='store_true', help='Show image previews with Matplotlib.')
	parser.add_argument('-s', '--stack', action='store_true', help='Stack the images in the preview window.')
	parser.add_argument('--noshift', action='store_true', help='Don''t apply a linear shift to the stacked images to correct for drift.')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf.')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only.')
	parser.add_argument('--startframe', default=1, type=int, help='Start frame. \'1\' is the default.')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames. No parameter means all frames, or from startframe to the end of the run.')
	parser.add_argument('-t', '--sleep', default=0, type=float, help='Sleep time (in seconds) between frames. \'0\' is the default.')
	parser.add_argument('-w', '--watch', type=int, help='Watch an aperture... specify the number.')
	parser.add_argument('--xyls', action='store_true', help='Write an XYLS (FITS) file output catalog that can be used as input to Astronomy.net.')
	parser.add_argument('--usefirstframe', action='store_true', help='Use the first frame of the run. Usually the first frame will be discarded.')
	parser.add_argument('-i', '--keyimages', action='store_true', help='Show some key images during the processing of this run.')
	parser.add_argument('--createconfig', action='store_true', help='Use this option to create a default configuration file.')
	parser.add_argument('--apertures', type=int, help='Number of reference apertures to use. The default is defined in the configuration file.')
	
	arg = parser.parse_args()
	
	if arg.createconfig:
		ultracamutils.createConfigFile()
		sys.exit()
	
	config = ultracamutils.readConfigFile(arg.configfile)
	
	innerSkyRadius = float(config.INNER_SKY)
	outerSkyRadius = float(config.OUTER_SKY)
	apertureRadius = float(config.APERTURE_RADIUS)
	polynomialDegree = int(config.POLY_DEGREE)
	numReferenceApertures = int(config.REF_APERTURES)
	if arg.apertures!=None: numReferenceApertures = arg.apertures
	applyShift = True
	if (arg.noshift):
		applyShift = False
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	debug.write(arg, level = 2)
	debug.write("Astropy version %s"%(astropy.__version__), level = 3)
	
	sourcesFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sources.csv"
	debug.write("Loading source list from: " + sourcesFilename, 2)	
	sourceList = ultraspecClasses.sourceList()
	success = sourceList.loadFromCSV(sourcesFilename)
	if (not success):
		debug.error("Unable to open the list of sources. Have you run uspecCreateSourceMap yet?")
		sys.exit()
	else:
		debug.write("Loaded %d sources from the CSV file."%sourceList.getNumSources(), 2)
	referenceApertures = ultraspecClasses.referenceApertures()
	referenceApertures.initFromSourceList(sourceList, max=numReferenceApertures)
	debug.write("Number of reference apertures we are going to use is %d."%len(referenceApertures.sources), 2)
	
	margins = 10
	
	runInfo = ultraspecClasses.runInfo(arg.runname)
	found = runInfo.loadFromJSON(config.RUNINFO)
	if not found:
		debug.write("Could not get info for this run from the ultra.json file.", level = 1)
		xmlRead = runInfo.loadFromXML(config.ULTRASPECRAW)
		
	debug.write(runInfo, 2)
		
	runFilename = ultracamutils.addPaths(config.ULTRASPECRAW, arg.runname)

	debug.write("Opening the Ultraspec raw file at: " + runFilename, level = 2)
	
	runDate, runID = ultracamutils.separateRunNameAndDate(arg.runname)
	
	""" Check that the working folders and the output folders are there
	"""
	(runDate, runNumber) = ultracamutils.separateRunNameAndDate(arg.runname)
	
	workingFolder = ultracamutils.addPaths(config.WORKINGDIR, runDate)
	ultracamutils.createFolder(workingFolder)
	
	startFrame = arg.startframe
	if startFrame<1:
		debug.error("startframe cannot be less than 1")
		sys.exit()

	rdat  = ultracam.Rdata(runFilename, startFrame, server=False)

	maximumFrames = rdat.ntotal()
	debug.write("Total number of frames in the run is %d"%maximumFrames, level = 2 )
	if startFrame>maximumFrames:
		debug.error("startframe " + str(startFrame) + ", is beyond the end of the run, which has only " + str(maximumFrames) + " frames in it.")
		sys.exit()
		
	frameRange = maximumFrames - startFrame + 1
	
	if arg.numframes!=None:
		requestedNumFrames = arg.numframes
		if requestedNumFrames<(frameRange):
			frameRange = requestedNumFrames
	
	startTime = datetime.datetime.now()
	timeLeftString = "??:??"
	""" Run through all the frames in the .dat file.
	"""
			
	fullFrame = numpy.zeros((1057, 1040))
			
	allWindows = []

	ccdFrame = rdat()
	frameWindows = ccdFrame[0]
		
	for windowIndex, w in enumerate(frameWindows):
		# Set up some info about the window sizes and extents
		window = ultraspecClasses.window()
		window.setExtents(w.llx, w.lly, w.nx, w.ny)
		window.setBinning(w.xbin, w.ybin)
		
		image = w._data
		if (arg.usefirstframe):
			window.setData(image)
			bkg_sigma = 1.48 * mad(image)
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			window.setSourcesAvoidBorders(sources)	
		else: 
			window.setBlankData(image)
		
		
		allWindows.append(window)
		
	(xmin, ymin, xmax, ymax) = determineFullFrameSize(allWindows)
	fullFramexsize = xmax - xmin
	fullFrameysize = ymax - ymin
	
	
	""" Set up the PGPLOT windows """
	xyPositionPlot = {}
	xyPositionPlot['pgplotHandle'] = ppgplot.pgopen('/xs')
	xyPositionPlot['yLimit'] = 1.0
	xyPositionPlot['numXYPanels'] = len(referenceApertures.sources)
	ppgplot.pgpap(6.18, 1.618)
	ppgplot.pgsubp(1, xyPositionPlot['numXYPanels'])
	ppgplot.pgsci(5)
	for panel in range(xyPositionPlot['numXYPanels']):
		currentSize = ppgplot.pgqch()
		ppgplot.pgsch(1)
		yLimit = xyPositionPlot['yLimit']
		ppgplot.pgenv(startFrame, startFrame + frameRange, -yLimit, yLimit, 0, -2)
		ppgplot.pgbox('A', 0.0, 0, 'BCG', 0.0, 0)
		ppgplot.pglab("", "%d"%panel, "")
		ppgplot.pgsch(currentSize)
	
	ppgplot.pgask(False)
	ppgplot.pgsci(1)
	
	
	if (arg.preview):		
		bitmapView = {}
		bitmapView['pgplotHandle'] = ppgplot.pgopen('/xs')
		ppgplot.pgpap(8, 1)
		
		ppgplot.pgenv(0.,fullFramexsize,0.,fullFrameysize, 1, 0)
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		ppgplot.pgsfs(2)
	
	if (arg.watch!=None) and (arg.watch<numReferenceApertures):
		watch = arg.watch
		watchView = {}
		watchView['pgplotHandle'] = ppgplot.pgopen('/xs')
		ppgplot.pgpap(10, 1)
		ppgplot.pgsvp(0.1, 0.7, 0.3, 0.9)
		ppgplot.pgswin(-margins, margins, -margins, margins)
		# ppgplot.pgenv(-margins, margins, -margins, margins, 1, 0)
		# ppgplot.pgenv(-margins, margins, 0, 10, 0, 0)
		watchView['pgPlotTransform'] = [-11, 1, 0, -11, 0, 1]
	else:
		watch=-1
	
	""" End of PGPLOT set up """
	
					
	xValues = []
	yValues = []	
	yAxisMax= 100	
	for frameIndex in range(2, frameRange + 1):
		framesToGo = frameRange - frameIndex
		currentTime = datetime.datetime.now()
		trueFrameNumber = startFrame + frameIndex - 1
		completionPercent = (float(frameIndex) / float(frameRange) * 100.)
		timePassed = ultracamutils.timedeltaTotalSeconds(currentTime - startTime)
		totalTime = timePassed * 100. / completionPercent
		etaTime = startTime + datetime.timedelta(seconds = totalTime)
		timeLeft = etaTime - currentTime
		(hours, mins, secs) = ultracamutils.timedeltaHoursMinsSeconds(timeLeft)
		timeLeftString = str(hours).zfill(2) + ":" + str(mins).zfill(2) + ":" + str(secs).zfill(2)
		
		ccdFrame = rdat()
		
		statusString = "\r%s Frame: [%d/%d]"%(timeLeftString, trueFrameNumber, frameRange)
		sys.stdout.write(statusString)
		sys.stdout.flush()
		
		windows = ccdFrame[0]
		for windowIndex, w in enumerate(windows):
			image = w._data		
			allWindows[windowIndex].setData(image)
		
		
		if arg.preview:
			ppgplot.pgslct(bitmapView['pgplotHandle']) 
			ppgplot.pgbbuf()
			fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
			for w in allWindows:
				if (arg.stack):
					boostedImage = ultracamutils.percentiles(w.stackedData, 20, 99)
				else:
					boostedImage = ultracamutils.percentiles(w.data, 20, 99)
				xll = w.xll/w.xbin - xmin
				xsize = w.nx
				yll = w.yll/w.ybin - ymin
				ysize = w.ny
				fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage		
			
			rows, cols  = numpy.shape(fullFrame)
			
			# Draw the grayscale bitmap
			ppgplot.pggray(fullFrame, 0, cols-1 , 0, rows-1 , 0, 255, pgPlotTransform)
	
			# Draw the full reference aperture list
			ppgplot.pgsci(3)
			for s in sourceList.getSources():
				(x, y) = s.abs_position
				ppgplot.pgcirc(x, y, 10)
		
			# ppgplot.pgslct(bitmapView)
			ppgplot.pgsci(2)

				
		for index, s in enumerate(referenceApertures.getSources()):
			window = allWindows[s.windowIndex]
			center = s.latestPosition
			xcenterInt = int(center[0])
			xcenterOffset = center[0] - margins
			#xcenterOffset = xcenterInt - margins
			
			ycenterInt = int(center[1])
			ycenterOffset = center[1] - margins
			#ycenterOffset = ycenterInt - margins
			
			if (ycenterOffset<0) or (xcenterOffset<0): continue
			zoomImageData = window.data[ycenterInt-margins:ycenterInt+margins, xcenterInt-margins:xcenterInt+margins]
			#(xcen, ycen) = photutils.morphology.centroid_2dg(zoomImageData, error=None, mask=None)
			
			xCollapsed = numpy.sum(zoomImageData, 0)
			yCollapsed = numpy.sum(zoomImageData, 1)
			
			xPeak = numpy.argmax(xCollapsed)
			yPeak = numpy.argmax(yCollapsed)
			if (yPeak==len(yCollapsed)-1) or (yPeak == 0) \
				or (xPeak==len(xCollapsed)-1) or (xPeak == 0):
				print "Peak too close to the edge [%d, %d] ... skipping this aperture [%d] for this frame [%d]."%(xPeak,yPeak, index, trueFrameNumber)
				continue

			# Fit quadratic polynomials to the collapsed profiles in the X-Y axis
			xMat = [ [ (xPeak-1)**2, (xPeak-1) , 1 ] , \
			         [ (xPeak)**2,   (xPeak)   , 1 ] , \
			         [ (xPeak+1)**2, (xPeak+1) , 1 ] ]
			yMat = [ xCollapsed[xPeak-1], xCollapsed[xPeak], xCollapsed[xPeak+1] ]
			(a, b, c) = numpy.linalg.solve(xMat, yMat)
			xPoly = (a, b, c)
			newxPeak = -1.0 * b / (2.0 * a)

			xMat = [ [ (yPeak-1)**2, (yPeak-1) , 1 ] , \
			         [ (yPeak)**2,   (yPeak)   , 1 ] , \
			         [ (yPeak+1)**2, (yPeak+1) , 1 ] ]
			yMat = [ yCollapsed[yPeak-1], yCollapsed[yPeak], yCollapsed[yPeak+1] ]
			(a, b, c) = numpy.linalg.solve(xMat, yMat)
			yPoly = (a, b, c)
			newyPeak = -1.0 * b / (2.0 * a)

			# Fit a Gaussian with pre-defined FWHM to the collapsed profiles.
			# X - direction
			fwhm = 2.5
			c = fwhm / 2.355
			b = 0.0
			baseLevel = numpy.median(xCollapsed)
			a = numpy.max(xCollapsed) - baseLevel
			d = baseLevel
			numGaussianPoints = len(xCollapsed)
			xGaussian = [float(x) * 2*margins/numGaussianPoints - margins for x in range(numGaussianPoints)]
			result, covariance = scipy.optimize.curve_fit(shifting_gaussian, xGaussian, xCollapsed, b)
			xBestOffset = result[0]
			xBestOffsetError = numpy.sqrt(numpy.diag(covariance))[0]
			debug.write("x-offset: %f [%f]"%(xBestOffset, xBestOffsetError), 3)
			xGaussianFit = [gaussian(x, a, xBestOffset, c, d) for x in xGaussian]
			xcen = xBestOffset + margins
			
			# Y - direction
			fwhm = 2.5
			c = fwhm / 2.355
			b = 0.0
			baseLevel = numpy.median(yCollapsed)
			a = numpy.max(yCollapsed) - baseLevel
			d = baseLevel
			numGaussianPoints = len(yCollapsed)
			yGaussian = [float(x) * 2*margins/numGaussianPoints - margins for x in range(numGaussianPoints)]
			result, covariance = scipy.optimize.curve_fit(shifting_gaussian, yGaussian, yCollapsed, b)
			yBestOffset = result[0]
			yBestOffsetError = numpy.sqrt(numpy.diag(covariance))[0]
			debug.write("y-offset: %f [%f]"%(yBestOffset, yBestOffsetError), 3)
			yGaussianFit = [gaussian(x, a, yBestOffset, c, d) for x in yGaussian]
			ycen = yBestOffset + margins
			
			#print "Centroid method: (%f, %f)   vs   Quadratic fit: (%f, %f)"%(xcen, ycen, newxPeak, newyPeak) 
			if index==watch:
				ppgplot.pgslct(watchView['pgplotHandle'])
				ppgplot.pgbbuf()
				ppgplot.pgeras()
				
				# Image of the watched source
				ppgplot.pgsvp(0.1, 0.7, 0.3, 0.9)
				ppgplot.pgswin(-margins, margins, -margins, margins)
				zRows, zCols = numpy.shape(zoomImageData)
				preview = ultracamutils.percentiles(zoomImageData, 20, 99)
				ppgplot.pggray(preview, 0, zCols-1, 0, zRows-1, 0, 255, watchView['pgPlotTransform'])
			
				# Graph of the x-direction
				ppgplot.pgsvp(0.1, 0.7, 0.1, 0.2)
				yMax = numpy.max(xCollapsed) * 1.1
				yMin = numpy.min(xCollapsed) * 0.9
				ppgplot.pgswin(-margins, margins, yMin, yMax)
				ppgplot.pgsci(1)
				ppgplot.pgbox('ABI', 1.0, 10, 'ABI', 0.0, 0)
				xPoints = [x - margins for x in range(len(xCollapsed))]
				ppgplot.pgsci(2)
				ppgplot.pgbin(xPoints, xCollapsed, True)
				numPolyPoints = 50
				xFit = [float(i) * len(xCollapsed)/numPolyPoints for i in range(numPolyPoints)]
				yFit = [xPoly[0]*x*x + xPoly[1]*x + xPoly[2] for x in xFit]
				ppgplot.pgsci(3)
				ppgplot.pgline([x - margins for x in xFit], yFit)
				ppgplot.pgsls(2)
				ppgplot.pgline([newxPeak-margins, newxPeak-margins], [yMin, yMax]) 
				ppgplot.pgsci(4)
				ppgplot.pgline([xBestOffset, xBestOffset], [yMin, yMax])
				ppgplot.pgsls(1)
				ppgplot.pgline(xGaussian, xGaussianFit)
				
				# Graph of the y-direction
				ppgplot.pgsvp(0.8, 0.9, 0.3, 0.9)
				yMax = numpy.max(yCollapsed) * 1.1
				yMin = numpy.min(yCollapsed) * 0.9
				ppgplot.pgswin(yMin, yMax, -margins, margins)
				ppgplot.pgsci(1)
				ppgplot.pgbox('ABI', 0.0, 0, 'ABI', 1.0, 10)
				xPoints = [x - margins for x in range(len(yCollapsed))]
				ppgplot.pgsci(2)
				ppgplot.pgbin(yCollapsed, xPoints, True)
				numPolyPoints = 50
				xFit = [float(i) * len(xCollapsed)/numPolyPoints for i in range(numPolyPoints)]
				yFit = [yPoly[0]*x*x + yPoly[1]*x + yPoly[2] for x in xFit]
				ppgplot.pgsci(3)
				ppgplot.pgline(yFit, [x - margins for x in xFit])
				ppgplot.pgsls(2)
				ppgplot.pgline([yMin, yMax], [newyPeak-margins, newyPeak-margins]) 
				ppgplot.pgsci(4)
				ppgplot.pgline([yMin, yMax], [yBestOffset, yBestOffset])
				ppgplot.pgsls(1)
				ppgplot.pgline(yGaussianFit, yGaussian)
				
				
				ppgplot.pgsci(1)
				ppgplot.pgebuf()
			
			xcen+= xcenterOffset
			ycen+= ycenterOffset
			
			xError = xBestOffsetError
			yError = yBestOffsetError
			s.setLatestPosition(trueFrameNumber, (xcen, ycen), errors = (xError, yError))
			apertures = CircularAperture((xcen, ycen), r=apertureRadius)
			annulus_apertures = CircularAnnulus((xcen, ycen), r_in=innerSkyRadius, r_out=outerSkyRadius)
			
			# Draw the re-positioned apertures
			xll = window.xll/window.xbin - xmin
			yll = window.yll/window.ybin - ymin
			if arg.preview:
				ppgplot.pgslct(bitmapView['pgplotHandle'])
				plotx= xcen + xll
				ploty= ycen + yll
				#print xll, yll, center, xcen, ycen
				ppgplot.pgcirc(plotx, ploty, apertureRadius)
				ppgplot.pgcirc(plotx, ploty, innerSkyRadius)
				ppgplot.pgcirc(plotx, ploty, outerSkyRadius)
				ppgplot.pgptxt(plotx-10, ploty-10, 0, 0, str(index)) 
				ppgplot.pgebuf()
	
		
		ppgplot.pgslct(xyPositionPlot['pgplotHandle'])
		for panel, aperture in enumerate(referenceApertures.getSources()):
			xPosition, yPosition = numpy.subtract(aperture.latestPosition, aperture.position)
			# Check if we need to re-scale the vertical axis
			yLimit = xyPositionPlot['yLimit']
			if abs(xPosition)>yLimit or abs(yPosition)>yLimit:
				yLimit*=1.2
				xyPositionPlot['yLimit'] = yLimit
				ppgplot.pgsubp(1, xyPositionPlot['numXYPanels'])
				ppgplot.pgsci(5)
				ppgplot.pgeras()
				for p in range(xyPositionPlot['numXYPanels']):
					ppgplot.pgenv(startFrame, startFrame + frameRange, -yLimit, yLimit, 0, 0)
					ppgplot.pgbox('A', 1.0, 10, 'BCG', 0.0, 0)
					ppgplot.pglab("", "%d"%p, "")
					ppgplot.pgsch(currentSize)
					
				for p, a in enumerate(referenceApertures.getSources()):
					xValues = [log['frameNumber'] for log in a.positionLog]
					yValues = [log['position'][0] - a.position[0] for log in a.positionLog]
					yErrors = [log['positionError'][0] for log in a.positionLog] 
					ppgplot.pgsci(2)
					ppgplot.pgpanl(1, p + 1)
					ppgplot.pgpt(xValues, yValues, 1)
					ppgplot.pgerry(xValues, yValues + yErrors, yValues + yErrors, 0)
					yValues = [log['position'][1] - a.position[1] for log in a.positionLog]
					ppgplot.pgsci(3)
					ppgplot.pgpanl(1, p + 1)
					ppgplot.pgpt(xValues, yValues, 1)
				
					
			shortXArray = [trueFrameNumber]
			shortYArray = [xPosition]
			ppgplot.pgsci(2)
			ppgplot.pgpanl(1, panel + 1)
			ppgplot.pgpt(shortXArray, shortYArray, 1)
			ppgplot.pgsci(3)
			shortYArray = [yPosition]
			ppgplot.pgpt(shortXArray, shortYArray, 1)
			
		if arg.sleep!=0:
			time.sleep(arg.sleep)
	sys.stdout.write("\rProcessed %d frames      \n"%frameRange)
	sys.stdout.flush()
	
	if arg.preview:
		ppgplot.pgslct(bitmapView['pgplotHandle'])
		ppgplot.pgclos()
		
	if watch!=-1:
		ppgplot.pgslct(watchView['pgplotHandle'])
		ppgplot.pgclos()
	
	# Sort the apertures by coverage
	referenceApertures.calculateFrameCoverage(frameRange-1)   # The '-1' is because we don't get photometry from the first frame
	referenceApertures.sortByCoverage()
	minCoverage = 80.0
	referenceApertures.limitCoverage(minCoverage)
	referenceApertures.sortByFlux()
	
	# Fit polynomials to the positions of the reference apertures
	for p, referenceAperture in enumerate(referenceApertures.getSources()):
		xValues = [log['frameNumber'] for log in referenceAperture.positionLog]
		yValues = [log['position'][0] - referenceAperture.position[0] for log in referenceAperture.positionLog]
		polynomial = numpy.polyfit(xValues, yValues, polynomialDegree)
		# Draw the polynomial
		ppgplot.pgslct(xyPositionPlot['pgplotHandle'])
		ppgplot.pgsci(2)
		ppgplot.pgpanl(1, p + 1)
		yValues = numpy.polyval(polynomial, xValues)
		ppgplot.pgpt(xValues, yValues, 1)
		referenceAperture.setPolynomial("x", polynomialDegree, polynomial)
	
		yValues = [log['position'][1] - referenceAperture.position[1] for log in referenceAperture.positionLog]
		polynomial = numpy.polyfit(xValues, yValues, polynomialDegree)
		# Draw the polynomial
		ppgplot.pgsci(3)
		ppgplot.pgslct(xyPositionPlot['pgplotHandle'])
		ppgplot.pgpanl(1, p + 1)
		yValues = numpy.polyval(polynomial, xValues)
		ppgplot.pgpt(xValues, yValues, 1)
		referenceAperture.setPolynomial("y", polynomialDegree, polynomial)
		
	
	ppgplot.pgslct(xyPositionPlot['pgplotHandle'])
	ppgplot.pgclos()
	
	
	# Write the reference aperture data as a CSV file
	referenceAperture = referenceApertures.getSources()[0]
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_reference_aperture.csv"
	outputFile = open(outputFilename, 'w')
	
	lineString = "Frame, Window, X, Y, X_ABS, Y_ABS\n"
	outputFile.write(lineString)
	apertureData = referenceAperture.getData()
	windowIndex = referenceAperture.windowIndex
	window = allWindows[windowIndex]
	xll = window.xll/window.xbin - xmin
	yll = window.yll/window.ybin - ymin
		
	for d in apertureData:
		xAbs = d[2] + xll
		yAbs = d[3] + yll
		lineString = "%d, %d, %f, %f, %f, %f\n"%(d[0], windowIndex, d[2], d[3], xAbs, yAbs)
		outputFile.write(lineString)
	outputFile.close()
	
	sys.exit()
	
	
