import numpy, json, classes
import astropy.io.fits
import os, subprocess, math, re
from PIL import Image,ImageDraw,ImageFont

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

def getRunsByDay(rawFilePath, date):
	
	try:
		fileList = os.listdir(rawFilePath + '/' + date) 
	except OSError:
		print "No data for the date: %s"%date
		return []
		
	runs_re = re.compile(r'run[0-9][0-9][0-9].xml')
	
	todaysRuns = []
	for f in fileList:
			r = runs_re.search(f)
			if (r):
				runNumber = r.group(0)[:6]
				runname = date + '/' + runNumber
				todaysRuns.append(runname)
	
	return todaysRuns
	

def readConfigFile(filename):
    """ Reads the config file name ucambuilder.conf and returns a configObject
    """
    configuration = classes.configObject()
    
    try:
        configfile = open(filename, "r")
    except IOError:
        print "Warning: Cannot find the config file %s ... will use default values"%filename
        return configuration

    for line in configfile:
        if(line[0]!="#"):      # Ignore lines that are comments
            if len(line.split())>=2:
                option = line.split()[0]
                value = line.split()[1]
                configuration[option] = value
	
    configfile.close()
    return configuration
    
def createConfigFile():
	configuration = classes.configObject()
	
	outputFile = open("ucambuilder_auto.conf", 'w')
	lineString = "# This file has been automatically genereated by 'ultracamutils.py'.\n"
	outputFile.write(lineString)
	attributes = [a for a in dir(configuration) if not a.startswith('__') and not callable(getattr(configuration,a))]
	print attributes
	for key in attributes:
		lineString = str(key) + "\t" + str(configuration[key]) + "\n"
		outputFile.write(lineString)
	outputFile.close()
	
    
def getRunInfo(filename, runIdent):
	""" Loads some run info from the JSON file created by Tom Marsh. Places it into an object of class runInfoObject and returns the object
	"""
	JSONfile = open(filename, "r")

	allObjectsJSON = json.load(JSONfile)

	(runDate, runNumber) = separateRunNameAndDate(runIdent)

	run = classes.runObject(runDate, runNumber)
	
	runNumberStr = runNumber[3:]
	runNumber = int(runNumberStr)
	
	for object in allObjectsJSON:
		date = object['night']
		num = object['num']
		if ((date == runDate) & (runNumber == num)):
			run.updateRunInfo(object)
			
	return run
	
    
def getRunMetaData(runName):
    """ Opens an Ultracam raw file and gets some metadata from it
    """
    rdat  = ultracam.Rdata(runFilename, startFrame, server=False)
    numFrames = rdat.ntotal()

def getUniqueID(objects):
	""" Returns a unique ID (monotonically increasing) for a list of objects that have a property called 'id'
	""" 
	newID = 0 
	for o in objects:
		if (o.id >= newID): newID = o.id + 1	
	return newID

def getObjectByID(objects, id):
	""" Returns an object that matches a given ID, assuming the list of objects have an ID property
	""" 
	for o in objects:
		if (o.id == id): return o	
	return None


def percentiles(data, lo, hi):
    """ Returns a normalised array where lo percent of the pixels are 0 and hi percent of the pixels are 255
    """
    max = data.max()
    dataArray = data.flatten()
    pHi = numpy.percentile(dataArray, hi)
    pLo = numpy.percentile(dataArray, lo)
    range = pHi - pLo
    scale = range/255
    data = numpy.clip(data, pLo, pHi)
    data-= pLo
    data/=scale
    return data

def measureDistance(p1, p2):
	""" Measures the 2D pythogorean distance between to (x,y) tuples
    """
	(x1, y1) = p1
	(x2, y2) = p2
	
	distance = math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) )
	
	return distance
	
    
def buildObjectsFromJSON(filename):
	""" Reads a JSON file and re-constructs the object list... returns it as an array
	"""
	JSONfile = open(filename, "r")

	wholeFileString = JSONfile.read()

	allObjectsJSON = json.loads(wholeFileString)

	objects = []

	for i in allObjectsJSON:
		ob = json.loads(i)
		# Create a new instance of ObservedObject for this object 
		newObject = classes.ObservedObject(ob['id'])
		newObject.currentPosition = (ob['x'], ob['y'])
		dataArray = ob['data']
		for data in dataArray:
			newObject.addExposure(data[2],data[3], data[1], data[4], data[0], 0)
		objects.append(newObject)
	
	JSONfile.close()
		
	return objects
	
def runnameToUniqueID(runname):
	""" Converts a runname to a unique ID number.
	"""
	retString = runname
	
	for char in "-/run":
		retString = retString.replace(char,'')
	
	return retString
	
	
def createFITS(frameNumber, window, channel, imageData, runname):
    """ Writes a FITS file for the images in a frame of the CCD
    """
    uniqueID = runnameToUniqueID(runname)
    filename = "/tmp/ucamwin" + str(uniqueID) + str(frameNumber).zfill(5)+ "_" + str(window) + "_" + channel + ".fits"
    hdu = astropy.io.fits.PrimaryHDU(imageData)
    hdulist = astropy.io.fits.HDUList([hdu])
    hdulist.writeto(filename, clobber=True)
    return filename
    
def saveFITSImage(imageData, filename):
    """ Writes a FITS file for the images in a frame of the CCD
    """
    hdu = astropy.io.fits.PrimaryHDU(imageData)
    hdulist = astropy.io.fits.HDUList([hdu])
    hdulist.writeto(filename, clobber=True)
    
    
def writeFITSwithRunHeaders(imageData, filename, runInfo):
	
	ra = runInfo.ra * 15.  # Convert RA to degrees
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
	prihdr['CRPIX1'] = 512
	prihdr['CRPIX2'] = 512
	prihdr['CRVAL1'] = ra
	prihdr['CRVAL2'] = dec
	prihdr['CDELT1'] = fieldScaleX
	prihdr['CDELT2'] = fieldScaleY
	
	hdu = astropy.io.fits.PrimaryHDU(imageData, header=prihdr)
	hdulist = astropy.io.fits.HDUList([hdu])
	hdulist.writeto(filename, clobber=True)

	

def removeTMPFile(filename):
    """ Removes a temporary file once sextractor has finished with it
    """
    os.remove(filename)


def removeFITS(number):
    """ Removes a temporary FITS file once sextractor has finished with it
    """
    filename = "/tmp/frame" + str(number) + ".fits"
    os.remove(filename)

def removeCAT(number):
    """ Removes a temporary CAT file once sextractor has finished with it
    """
    filename = "/tmp/frame" + str(number) + ".cat"
    os.remove(filename)


def runSex(tmpFilename, apertures = 'False'):
	""" Runs a sextractor process for an image that has just been placed in a FITS file
	"""
	sexCommand = ["sex"]
	sexCommand.append(tmpFilename)
    
	catFilename = tmpFilename[:-5] + ".cat"
	catFileParameter = "-CATALOG_NAME " + str(catFilename)
	sexCommand.append(catFileParameter)
    
	if(apertures==True):
		apertureDirective = "-CHECKIMAGE_TYPE APERTURES"
		sexCommand.append(apertureDirective)
		apertureFiles =  "-CHECKIMAGE_NAME " + tmpFilename[:-5] + ".aperture.fits"
		sexCommand.append(apertureFiles)
		
	subprocess.call(sexCommand)
	return catFilename
    

def readSexObjects(catFilename, sexMagnitude):
	""" Reads a sextractor .cat file and returns a list of objects in that file
    """
	sexCatalog = astropy.io.fits.open(catFilename)
	headers = sexCatalog["LDAC_OBJECTS"].header
	data = sexCatalog["LDAC_OBJECTS"].data
	columns = sexCatalog["LDAC_OBJECTS"].columns
	objects = []
	
	for item in data:
		object = {}
		object['id']     = item[columns.names.index("NUMBER")]
		object['x']      = item[columns.names.index("X_IMAGE")]
		object['y']      = item[columns.names.index("Y_IMAGE")]
		#object['counts'] = item[columns.names.index("FLUX_APER")]
		object['counts'] = item[columns.names.index(sexMagnitude)]
		object['mag']    = item[columns.names.index("MAG_AUTO")]
		object['radius'] = item[columns.names.index("FLUX_RADIUS")]
		object['flags']  = item[columns.names.index("FLAGS")]
		objects.append(object)
	
	sexCatalog.close()
	
	return objects
	
def filterOutCosmicRays(objects):
	""" Returns a reduced list of objects with the cosmic rays removed
	"""
	filteredList = []
	for i in objects:
		if not i.isCosmicRay():
			filteredList.append(i)
	
	return filteredList

def filterOutPixels(objects, pixelThreshold):
	""" Returns a reduced list of objects with the 'single pixel' objects
	"""
	filteredList = []
	for i in objects:
		if i.meanFWHM > pixelThreshold:
			filteredList.append(i)
	
	return filteredList

		
def filterOutLowFrameCountObjects(objects, percentage):
	""" For a given list of objects, filters out those objects that do not appear on > percentage of frames
	"""
	maxFrames = 0
	for i in objects:
		if (i.numExposures>maxFrames): maxFrames = i.numExposures;
	
	filteredObjects = []
	threshold = maxFrames * percentage/100.
	for i in objects:
		if i.numExposures>=threshold: filteredObjects.append(i)
		
	return filteredObjects

def filterOutBadTimingFrames(objects):
	""" For a given list of objects, checks for frames with weird times (like the first few frames in a DRIFT mode run, or just some random jumps in time) and removes them from the data
	This approach is very inefficient since it checks all objects individually, when timing errors are more likely to be the same for all objects. Needs to be optimised.
	"""
	
	for o in objects:
		startTime = o.exposures[0].MJD
		endTime = o.exposures[-1].MJD
		
		newExposures = []
		for f, frame in enumerate(o.exposures):
			time = frame.MJD
			difference = time - startTime
			if (time>startTime) and (time<endTime): 
				newExposures.append(frame)
						
		o.resetExposures(newExposures)
				
	return objects		
		

def rejectBadObjects(objects):
    """ Returns a filtered set of objects without the sextractor bad objects
    """
    filteredObjects = []
    for ob in objects:
      if (ob['flags'] == 4): filteredObjects.append(ob)
      if (ob['flags'] == 0): filteredObjects.append(ob)
      if (ob['flags'] == 2): filteredObjects.append(ob)
    return filteredObjects

def getTopObjects(objectList, num):
	""" Returns the brightest n objects in the field
	"""
	brightest = []
	sortedList = sorted(objectList, key=lambda object: object.lastCounts, reverse = True)
	
	brightest = sortedList[0:num]
	
	return brightest
	
def calculateDistance(p1, p2):
	xd = (p1[0] - p2[0])
	yd = (p1[1] - p2[1])
	distance = math.sqrt(xd*xd + yd*yd)
	return distance 
	
def computeOffsetArray(arrayObject, dx, dy):
	""" Returns a new array with the values linearly interpolated by the offset amount
	"""
	# Start with dx
	row = arrayObject[0]
	print row
	outputRow = numpy.zeros(shape = len(row))
	
	for i in range(len(row)-1):
		newXvalue = i + dx
		print newXvalue, math.floor(newXvalue), math.ceil(newXvalue)
		y1 = row[math.floor(newXvalue)]
		y2 = row[math.ceil(newXvalue)]
		m = y2-y1
		c = y1 - m*i
		print "m,c",m,c, 
		ynew = (newXvalue) * m + c
		print "ynew", ynew
		outputRow[i] = ynew	
		
	print row
	print outputRow
	
def separateRunNameAndDate(name):
	""" Takes a string like 2013-07-21/run111 and returns a tuple with the date and the runName separated
	"""
	runs_re = re.compile(r'run[0-9]{3}')
	date_re = re.compile(r'20[0-9]{2}(-[0-9]{2}){2}')
		
	date = ""
	runName = ""
	
	d = date_re.search(name)
	if (d):
		date = d.group(0)

	r = runs_re.search(name)
	if (r):
		runName = r.group(0)
	
	return (date, runName) 
	
def addPaths(path1, path2):
	""" Adds two paths together inserting a '/' if required
	"""
	path = path1
	if (path[-1]!='/'): path+= '/';
	path+= path2
	return path
	
def createFolder(path):
	""" Creates a folder if it does not exist already
	"""
	if not os.path.exists(path):
		os.mkdir(path)
	
def writePNG(image, filename, caption = ""):
    """ Writes to a PNG file using the PIL library. Adds a caption if sent in the parameters. Also adds a .png extension if it isn't already there in 'filename'
	"""
    imageCopy = image.copy()    # We need to copy the image so we don't alter the original when adding the caption.
    config = readConfigFile()
    if (caption!=""): 
		font = ImageFont.truetype(config.FONT, 25) 
		draw = ImageDraw.Draw(imageCopy)
		if (imageCopy.mode == "L"):
			draw.text((0, 0), caption, 255, font = font)
		else: 
			draw.text((0, 0), caption, (255, 255, 255), font = font)
	
    if (filename[-4:]!=".png"): filename+= ".png"
    imageCopy.save(filename, "PNG")
	
def timedeltaTotalSeconds(td):
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / float(10**6)
    
def timedeltaHoursMinsSeconds(td):
	hours, remainder = divmod(td.seconds, 3600)
	minutes, seconds = divmod(remainder, 60)
	return (hours, minutes, seconds)
	
def writeFriendlyTimeMinutes(minutes):
	""" Writes a friendly time (in hours and/or minutes) based on an input of minutes
	"""
	timeStr = str(int(minutes)) + " minute"
	if minutes > 60: 
		hours = minutes/60.
		minutes = minutes - int(hours)*60
		timeStr= "%d hour"%(int(hours))
		if int(hours)!=1: timeStr+="s";
		timeStr+= ", %d minute"%(int(minutes))
	if int(minutes)!=1: timeStr+="s";
	
	return timeStr

def fromSexagesimal(raStr, decStr):
	""" Format for input ra and dec are 'HH:MM:SS.dd' and 'nDD:MM:SS.dd'
									or 	'HH MM SS.dd' and 'nDD MM SS.dd'
	"""
	separator = ':'
	if raStr.find(separator)==-1:
		separator = ' '
	raPieces = raStr.split(separator)
	raHours = int(raPieces[0])
	raMinutes = int(raPieces[1])
	raSeconds = float(raPieces[2])
	ra = 15 * (raHours + raMinutes/60.0 + raSeconds / 3600.0)
	
	decPieces = decStr.split(separator)
	if decPieces[0][0]=='-':
		south = True
	else:
		south = False
	decHours = int(decPieces[0])
	decMinutes = int(decPieces[1])
	decSeconds = float(decPieces[2])
	dec = decHours + decMinutes/60.0 + decSeconds / 3600.0
	if south:
		dec = decHours - decMinutes/60.0 - decSeconds / 3600.0
		
	return (ra, dec)

def toSexagesimal(world):
	raDeg = world[0]
	ra = raDeg/15.
	hours = int(ra)
	minutes = (ra - int(ra)) * 60
	seconds = (minutes - int(minutes)) * 60
				
	dec = world[1]
	decDegrees = int(dec)
	decMinutes = (dec - int(dec)) * 60
	decSeconds = (decMinutes - int(decMinutes)) * 60
		
	outString = "RA: %02d:%02d:%02.1f"%(hours, minutes, seconds)
	outString+= " DEC: %02d:%02d:%02.3f"%(dec, decMinutes, decSeconds)
	return outString

def writeFriendlyTimeSeconds(seconds):
	""" Writes a friendly time (in hours and/or minutes) based on an input of seconds
	"""
	minutes = seconds / 60.
	timeStr = str(int(minutes)) + " minutes"
	if minutes > 60: 
		hours = minutes/60.
		minutes = minutes - int(hours)*60
		timeStr= "%d hour"%(int(hours))
		if int(hours)!=1: timeStr+="s";
		timeStr+= ", %d minute"%(int(minutes))
	if int(minutes)!=1: timeStr+="s";
			
	return timeStr
