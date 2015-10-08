#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
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
import ppgplot

def shift_func(output_coords, xoffset, yoffset):
	return (output_coords[0] - yoffset, output_coords[1] - xoffset)

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
	
	parser = argparse.ArgumentParser(description='Reads the ULTRASPEC [dd-mm-yyyy/runxxx.dat] files and prepares the run for aperture identification. Also produces previews of the images.')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-p', '--preview', action='store_true', help='Show image previews with Matplotlib')
	parser.add_argument('-s', '--stack', action='store_true', help='Stack the images in the preview window')
	parser.add_argument('--noshift', action='store_true', help='Don''t apply a linear shift to the stacked images to correct for drift.')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('--startframe', default=1, type=int, help='Start frame. \'1\' is the default')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames. No parameter means all frames, or from startframe to the end of the run')
	parser.add_argument('-t', '--sleep', default=0, type=int, help='Sleep time (in seconds) between frames. \'0\' is the default')
	parser.add_argument('--xyls', action='store_true', help='Write an XYLS (FITS) file output catalog that can be used as input to Astronomy.net')
	parser.add_argument('--usefirstframe', action='store_true', help='Use the first frame of the run. Usually the first frame will be discarded.')
	parser.add_argument('-i', '--keyimages', action='store_true', help='Show some key images during the processing of this run.')
	
	arg = parser.parse_args()
	
	config = ultracamutils.readConfigFile(arg.configfile)
	
	applyShift = True
	if (arg.noshift):
		applyShift = False
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	debug.write(arg, level = 2)
	debug.write("Astropy version %s"%(astropy.__version__), level = 3)
	
	
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
		
	if maximumFrames<10:
		debug.error("The total number of frames in this run is less than 10. We need more to create a source map. Exiting.")
		sys.exit()
	
	""" We are going to use the first 10 frames to build the original source map that is merely used as a guide for creating the main stacked image.
	"""	
	
	debug.write("Examining the first 10 frames of the run...", level = 2)
	originalStackedFrame = numpy.zeros((1057, 1040))
	
	ccdFrame = rdat()
	window = ccdFrame[0]
	allWindows = []	
	for windowIndex, w in enumerate(window):
		# Set up some info about the window sizes and extents
		window = ultraspecClasses.window()
		window.setExtents(w.llx, w.lly, w.nx, w.ny)
		window.setBinning(w.xbin, w.ybin)
		
		image = w._data
		if (arg.usefirstframe):
			window.setData(image)
			bkg_sigma = 1.48 * mad(image)
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			window.setSources(sources)	
		else: 
			window.setBlankData(image)
		
		allWindows.append(window)
		
	(xmin, ymin, xmax, ymax) = determineFullFrameSize(allWindows)
	fullFramexsize = xmax - xmin
	fullFrameysize = ymax - ymin
	
	# Stack up the next 10 frames - (numbers 2 to 11)
	for frame in range(10):
		ccdFrame = rdat()
		frameWindows = ccdFrame[0]
		
		for windowIndex, w in enumerate(frameWindows):
			image = w._data
			allWindows[windowIndex].addData(image)
			
		
	# Reconstruct a full frame from the windows	
	boostedFullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		boostedImage = ultracamutils.percentiles(w.stackedData, 10, 99.8)
		image = w.stackedData
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		boostedFullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + image
		
		bkg_sigma = 1.48 * mad(image)
		sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma) 
		w.setSourcesAvoidBorders(sources)	
		
		
	# Get the source list from this image
	# Combine the sources from all of the windows
	allSources = []
	for index, w in enumerate(allWindows):
		xll = w.xll/w.xbin - xmin
		yll = w.yll/w.ybin - ymin
		sources = w.getSources()
		positions = zip(sources['xcentroid'], sources['ycentroid'], sources['flux'])
		new_positions = [(x + xll, y + yll, flux) for (x, y, flux) in positions]
		allSources+=new_positions
		
		
	# Sort these sources in order of brightness and take the top 40%
	allSources = sorted(allSources, key=lambda object: object[2], reverse = True)
	numSources = len(allSources)
	maxSources = int(round((numSources)*0.4))
	debug.write("Number of sources: %d, number of top sources: %d"%(numSources, maxSources), 2)
	if maxSources<1:
		debug.write("WARNING: Not enough sources for shift calculation, proceeding in '--noshift' mode.", 1)
		applyShift = False
	else:
		topSources = allSources[0:maxSources]
		masterApertureList = [ (x, y) for (x, y, flux) in topSources]
		
	
	# Plot the preview frame
	if (arg.keyimages):
		#stackedFigure = matplotlib.pyplot.figure(figsize=(10, 10))
		#matplotlib.pyplot.title("Initial 10 frame stacked image")
		stackedPreview = ppgplot.pgopen('/xs')
		ppgplot.pgenv(0.,fullFramexsize,0.,fullFrameysize, 1, 0)
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		ppgplot.pglab("x", "y", "Initial 10 frame stacked image.")
	
		# Display the image on the user's screen
		image = matplotlib.pyplot.imshow(boostedFullFrame, cmap='gray_r')
		for s in allSources:
			x, y = s[0], s[1]
			matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 10, color='green', fill=False, linewidth=1.0))
		if applyShift:
			for s in topSources:
				x, y = s[0], s[1]
				matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 10, color='blue', fill=False, linewidth=1.0))
		
		rows, cols = numpy.shape(boostedFullFrame)
		ppgplot.pggray(boostedFullFrame, 0, cols-1, 0, rows-1, 0, 255, pgPlotTransform)
		ppgplot.pgsfs(2)   # Set fill style to 'outline'
		ppgplot.pgsci(3)   # Set the colour to 'green'
		for s in allSources:
			x, y = s[0], s[1]
			ppgplot.pgcirc(x,y, 10)
		
		
	""" End of the prework """

	rdat.set(1)		# Reset back to the first frame
	
	frameRange = maximumFrames - startFrame + 1
	
	if arg.numframes!=None:
		requestedNumFrames = arg.numframes
		if requestedNumFrames<(frameRange):
			frameRange = requestedNumFrames
	
	startTime = datetime.datetime.now()
	timeLeftString = "??:??"
	""" Run through all the frames in the .dat file.
	"""
	if arg.preview:
		mainPreview = ppgplot.pgopen('/xs')
		ppgplot.pgenv(0.,fullFramexsize,0.,fullFrameysize, 1, 0)
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		if arg.stack:
			ppgplot.pglab("x", "y", "Stacked image.")
		else:
			ppgplot.pglab("x", "y", "Frame image.")
		if applyShift:
			zoomedAperture = ppgplot.pgopen('/xs')
			ppgplot.pgenv(0., 20 ,0., 20, 1, 0)
			pgPlotTransform = [0, 1, 0, 0, 0, 1]
			ppgplot.pglab("x", "y", "Zoom on reference source.")		
			
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
			bkg_sigma = 1.48 * mad(image)
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			allWindows[windowIndex].setSourcesAvoidBorders(sources)	
			
		# Combine the sources in all of the windows
		allSources = []
		for index, w in enumerate(allWindows):
			xll = w.xll/w.xbin - xmin
			yll = w.yll/w.ybin - ymin
			sources = w.getSources()
			positions = zip(sources['xcentroid'], sources['ycentroid'], sources['flux'])
			new_positions = [(x + xll, y + yll, flux) for (x, y, flux) in positions]
			allSources+=new_positions

		allSources = sorted(allSources, key=lambda object: object[2], reverse = True)
		# Remove the flux column from the source list. We don't need it anymore. 
		tempSources = [ (x, y) for (x, y, flux) in allSources]
		allSources = tempSources
		
		if (applyShift):
			oldCatalog = numpy.array(masterApertureList)
			newCatalog = numpy.array(allSources)

			psize  = 0.1
			fwhm   = 4.
			dmax   = 10.
			mmax   = 10.

			(gaussImage, xp, yp, xr, yr) = ultracam_shift.vimage(oldCatalog, newCatalog, dmax, psize, fwhm)
			debug.write("Applying offset: (%2.2f, %2.2f)"%(xr, yr), level = 3)

		for windowIndex, w in enumerate(windows):
			image = w._data
			allWindows[windowIndex].setData(image)
			if (applyShift): image = ndimage.interpolation.shift(image, (-1.0*yr, -1.0*xr), order = 4 )
			allWindows[windowIndex].addToStack(image)
				
		if arg.preview: 
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
					
			ppgplot.pgslct(mainPreview)
			(rows, cols) = numpy.shape(fullFrame)
			ppgplot.pggray(fullFrame, 0, cols-1, 0, rows-1, 0, 255, pgPlotTransform)
			ppgplot.pgsci(3)
			ppgplot.pgsfs(2)
			for s in allSources:
				(x, y) = s
				ppgplot.pgcirc(x, y, 10)
				
			# Now also draw the zoomed in region around the first aperture
			if (applyShift):
				(x, y) = masterApertureList[0]
				croppedFrame = fullFrame[y-10:y+10, x-9:x+10]
				ppgplot.pgslct(zoomedAperture)
				rows, cols = numpy.shape(croppedFrame)
				ppgplot.pggray(croppedFrame, 0, cols-1, 0, rows-1, 0, 255, pgPlotTransform)
				ppgplot.pgsci(3)
				ppgplot.pgsfs(2)
				ppgplot.pgcirc(11, 11, 1)
				xLine = [11, 11+xr]
				yLine = [11, 11+yr]
				ppgplot.pgline(xLine, yLine)
				
			
			
			
		if arg.sleep!=0:
			time.sleep(arg.sleep)
	sys.stdout.write("\rProcessed %d frames      \n"%frameRange)
	sys.stdout.flush()
	
	ppgplot.pgclos()
		
	################################################################################################
	""" We have run through all of the images now. """
	################################################################################################
	
	allSources = []
	sourceList = ultraspecClasses.sourceList()
	for index, w in enumerate(allWindows):
		xll = w.xll/w.xbin - xmin
		yll = w.yll/w.ybin - ymin
		image = w.stackedData
		mean, median, std = sigma_clipped_stats(image, sigma=3.0)
		maximum = numpy.max(image)
		minimum = numpy.min(image)
		debug.write("Mean: %f,  Median: %f,  Std (clipped 3sigma):%f"%(mean, median, std) , 2)
		debug.write("Minimum: %f,  Maximum: %f"%(minimum, maximum), 2)
		threshold = median + (std * 2.)
		segm_img = detect_sources(image, threshold, npixels=5)
		mask = segm_img.astype(numpy.bool)
		mean, median, std = sigma_clipped_stats(image, sigma=3.0, mask=mask) 
		debug.write("After source masking", 2)
		debug.write("Mean: %f,  Median: %f,  Std (clipped 3sigma): %f"%(mean, median, std), 2)
		selem = numpy.ones((5, 5))    # dilate using a 5x5 box
		mask2 = binary_dilation(mask, selem)
		mean, median, std = sigma_clipped_stats(image, sigma=3.0, mask=mask2)
		debug.write("After dilation", 2)
		debug.write("Mean: %f,  Median: %f,  Std (clipped 3sigma): %f"%(mean, median, std), 2)
		
		# Check the window image for any areas that should be masked...
		lowerLimitBkg = median - std*5.
		debug.write("5 sigma below the median is the lowerLimitBkg for the mask: %f"%(lowerLimitBkg), 2)
		mask = (image < lowerLimitBkg)
		maskBitmap = numpy.zeros(numpy.shape(mask))
		maskBitmap = 255 * (mask) 
		if (arg.keyimages):
			maskImage = matplotlib.pyplot.figure(figsize=(10, 10))
			matplotlib.pyplot.title("Mask for window:%d"%(index))
			matplotlib.pyplot.imshow(maskBitmap, origin='lower', cmap='Greys_r',  interpolation = 'nearest')
			matplotlib.pyplot.show(block=False)
		bkg = Background(image, (10, 10), filter_shape=(3, 3), method='median', mask=mask)
		background = bkg.background
		if (arg.keyimages):
			bgImage = matplotlib.pyplot.figure(figsize=(10, 10))
			matplotlib.pyplot.title("Fitted background, window:%d"%(index))
			matplotlib.pyplot.imshow(background, origin='lower', cmap='Greys_r',  interpolation = 'nearest')
			matplotlib.pyplot.show(block=False)
		image = image - background
		
		# Final stage source detection
		bkg_sigma = 1.48 * mad(image)
		sigmaThreshold = float(config.SIGMA_THRESHOLD)* float(bkg_sigma)
		debug.write("Threshold for source detection is %f sigma or %f counts."%(float(config.SIGMA_THRESHOLD), sigmaThreshold), 2)
		sources = daofind(image, fwhm=4.0, threshold=sigmaThreshold)   
		
		w.setSourcesAvoidBorders(sources)
		w.BGSubtractedImage = image	
		
		sources = w.getSources()
		for s in sources:
			position = (s['xcentroid'], s['ycentroid'])
			sourceObject = ultraspecClasses.source(0, position, index)
			sourceObject.setDAOPhotData(s['sharpness'], s['roundness1'], s['roundness2'], s['npix'], s['sky'], s['peak'], s['flux'], s['mag'])
			sourceObject.setOffsets((xll, yll))
			sourceList.addSource(sourceObject)
		positions = zip(sources['xcentroid'], sources['ycentroid'], sources['flux'])
		new_positions = [(x + xll, y + yll, flux) for (x, y, flux) in positions]
		allSources+=new_positions
		
	# Get the final stacked image
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		imageData = w.stackedData
		boostedImageData = ultracamutils.percentiles(w.BGSubtractedImage, 40, 99.8)
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + imageData
		boostedFullFrame[yll:yll+ysize, xll:xll+xsize] = boostedFullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImageData

	
	#tempSources = [ (x, y) for (x, y, flux) in allSources]
	#allSources = tempSources	
	allSources = sorted(allSources, key=lambda object: object[2], reverse = True)
	if (arg.keyimages):
		finalFigure = matplotlib.pyplot.figure(figsize=(10, 10))
		matplotlib.pyplot.title("Final stacked image")
		matplotlib.pyplot.imshow(boostedFullFrame, cmap='gray_r')
		for s in allSources:
			(x, y) = s[0], s[1]
			matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 5, color='blue', fill=False, linewidth=1.0))
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show()
	
	# Output the source list for debug purposes
	if (arg.debuglevel>1):
		sourceString = "Sources"
		for s in allSources:
			sourceString+= "\n(%3.2f, %3.2f) %.2f"%(s[0], s[1], s[2])
		debug.write(sourceString, 2)
		
	
	sourceList.sortByFlux()
	sourcesFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sources.csv"
	debug.write("Writing source list to CSV file: " + sourcesFilename, 2)
	sourceList.writeToCSV(sourcesFilename)
	
	# Write the XYLS FITS file
	if (arg.xyls):
		IDs = []
		x_values = []
		y_values = []
		fluxes = []
		
		for num, s in enumerate(allSources):
			IDs.append(num)
			x_values.append(s[0])
			y_values.append(s[1])
			fluxes.append(s[2])
			

		FITSFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sources.xyls"
		debug.write("Writing FITS file: " + FITSFilename, level=2)
		col1 = astropy.io.fits.Column(name='ID', format='I', array=IDs)
		col2 = astropy.io.fits.Column(name='X', format='E', array=x_values)
		col3 = astropy.io.fits.Column(name='Y', format='E', array=y_values)
		col4 = astropy.io.fits.Column(name='FLUX', format='E', array=fluxes)
		cols = astropy.io.fits.ColDefs([col1, col2, col3, col4])	
		#tbhdu =astropy.io.fits.new_table(cols)
		tbhdu =astropy.io.fits.TableHDU.from_columns(cols)
		
		prihdr = astropy.io.fits.Header()
		prihdr['TARGET'] = runInfo.target
		prihdr['RA'] = runInfo.ra
		prihdr['DEC'] = runInfo.dec
		prihdr['COMMENT'] = "This file created by uspecCreateSourceMap.py from the Ultracam pipeline."
		prihdr['RUNIDENT'] = arg.runname
		prihdr['WIDTH'] = fullFramexsize
		prihdr['HEIGHT'] = fullFrameysize
		
		prihdu = astropy.io.fits.PrimaryHDU(header=prihdr)
		thdulist = astropy.io.fits.HDUList([prihdu, tbhdu])
		thdulist.writeto(FITSFilename, clobber=True)
	
	# Generate the stacked image for writing to disc
	stackedFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Stacked image")
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		boostedImage = ultracamutils.percentiles(w.BGSubtractedImage, 10, 99.5)
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + ".png"
	# Write the image data with PIL library, rather than matplotlib
	imgData = numpy.rot90(fullFrame, 3)
	imgSize = numpy.shape(imgData)
	imgLength = imgSize[0] * imgSize[1]
	testData = numpy.reshape(imgData, imgLength, order="F")
	img = Image.new("L", imgSize)
	palette = []
	for i in range(256):
		palette.extend((i, i, i)) # grey scale
		img.putpalette(palette)
	img.putdata(testData)
	debug.write("Writing PNG file: " + outputFilename, level = 2) 
	img.save(outputFilename, clobber=True)
	
	palette = []
	for i in range(256):
		palette.extend((255-i, 255-i, 255-i)) # inverse grey scale
		img.putpalette(palette)
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_inverted.png"
	debug.write("Writing PNG file: " + outputFilename, level = 2) 
	img.save(outputFilename, "PNG", clobber=True)
	
	# Write out the stacked image as a non-normalised FITS image
	FITSFilename =  ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_stacked.fits"
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		imageData = w.BGSubtractedImage
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + imageData

	ra = runInfo.ra  # Convert RA to degrees
	dec = runInfo.dec
	fieldScaleX = -8.3E-05
	fieldScaleY = 8.3E-05
	
	prihdr = astropy.io.fits.Header()
	prihdr['COMMENT'] = "This file created by the Ultracam pipeline."
	prihdr['TARGET'] = runInfo.target
	prihdr['COMMENT'] = runInfo.comment
	prihdr['EQUINOX'] = 2000
	prihdr['RADECSYS'] = "FK5"
	prihdr['CTYPE1'] = "RA---TAN"
	prihdr['CTYPE2'] = "DEC--TAN"
	prihdr['CRPIX1'] = fullFramexsize/2
	prihdr['CRPIX2'] = fullFrameysize/2
	prihdr['CRVAL1'] = ra
	prihdr['CRVAL2'] = dec
	prihdr['CDELT1'] = fieldScaleX
	prihdr['CDELT2'] = fieldScaleY
	
	hdu = astropy.io.fits.PrimaryHDU(fullFrame, header=prihdr)
	
	hdulist = astropy.io.fits.HDUList([hdu])
	hdulist.writeto(FITSFilename, clobber=True)
			
	sys.exit()
	
