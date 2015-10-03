import math, json, datetime, os
import trm.ultracam
import rashley_utils as utils
import ultracamutils
import numpy

class dayObject:
	""" This class stores metadata for a particular night's observing
	"""
	def __init__(self, date):
		self.date = date
		self.runs = []
		self.totalTime = 0
		self.totalFrames = 0
        
	def addRun(self, runName):
		run = runObject(self.date, runName)
		self.runs.append(run)
        
	def getRuns(self):
		return self.runs
	
		
	def setRuns(self, runData):
		self.runs = runData

class runObject:
	def __init__(self, date, runName):
		self.runID = runName
		self.runDate = date
		config = utils.readConfigFile()
		runPath = utils.addPaths(config.ULTRACAMRAW, date)
		runPath = ultracamutils.addPaths(runPath, runName)
		self.totalTime = 0
		runMetaData = trm.ultracam.Rhead(runPath, server=False)
		self.mode = runMetaData.mode
		#self.userData = runMetaData.user
		try: 
			self.nblue = runMetaData.nblue
		except AttributeError:
			self.nblue = 1
		if (self.mode != "PONOFF"):
			runData =  trm.ultracam.Rdata(runPath, 1, server=False)
			self.numFrames = runData.ntotal()
		else: 
			self.numFrames = 0
		self.runClass = 0
		self.comment = ""
		self.ra = 0
		self.dec = 0
		self.objectID = "?"
		self.target = "?"
		self.expose = 0
		try: 
			self.exposeTime = runMetaData.exposeTime 
		except:
			self.exposeTime = -1
		self.num = 0
		self.dataProcessed = False
		self.numWindows = 0
		self.maxExtents = [0, 0, 0, 0]
		self.version = 'primary'
		
	def loadSelf(self, configData):
		""" Checks to see if a saved file exists that contains info for this run. Loads it into this object.
		"""
		filename = ultracamutils.addPaths(configData.SITE_PATH, self.runDate)
		filename = ultracamutils.addPaths(filename, self.runID)
		if (self.version!='primary'): filename+= "_" + str(self.version)
		filename+= "_info.json"
		print "Trying to load myself from a file.", filename
		if os.path.exists(filename):
			JSONfile = open(filename, "r")
			wholeFileString = JSONfile.read()
			parsedObject = json.loads(wholeFileString)
			for key in parsedObject.keys():
				print "Setting property:", key, parsedObject[key]
				setattr(self,key,parsedObject[key])
		else: 
			print "Not found... falling back to ultra.json"
			self.mergeULTRAJSON(configData.RUNINFO)
			
			
	def mergeULTRAJSON(self, ultrajsonFilename):
		""" Looks in Tom's ultra.json file and gets the data there. Merges it with this object
		"""
		JSONfile = open(ultrajsonFilename, "r")
		allObjectsJSON = json.load(JSONfile)
		run = {}
		runNumberStr = self.runID[3:]
		runNumber = int(runNumberStr)
	
		for object in allObjectsJSON:
			date = object['night']
			num = object['num']
			if ((date == self.runDate) & (runNumber == num)):
				self.comment = object["comment"]
				self.ra = object['ra']
				self.dec = object['dec']
				self.objectID = object['id']
				self.target = object['target']
				self.num = object['num']
				self.expose = object['expose']
			
	def checkForComments(self, rawDataPath):
		""" Does a check for comments in the DDDD-MM-YY.dat file in the raw data path
		"""
		filename = ultracamutils.addPaths(rawDataPath, self.runDate)
		filename = ultracamutils.addPaths(filename, self.runDate) + '.dat'
		dataFile = open(filename, 'r')
		
		for line in dataFile:
			runIdentifier = line[:6]
			if (runIdentifier==self.runID):
				self.comment = line[7:]
			
	def addSexInfo(self, config):
		""" Adds info about the sextractor parameters to this object
		"""
		self.sexMagnitude = config.SEX_MAGNITUDE;
		allSexOptions = {}
		self.sexOptions = {'ANALYSIS_THRESH': None, 'SATUR_LEVEL': None, 'PHOT_APERTURES': None, 'DETECT_MINAREA': None, 'DETECT_THRESH': None, \
		                   'PHOT_AUTOPARAMS': None, 'BACK_SIZE': None, 'BACK_FILTERSIZE': None}
		# Now read the default.sex file for more parameters
		try:
			sexConfigFile = open("default.sex", "r")
		except IOError:
			print "Warning: Cannot find the sextractor 'default.sex' file %s ... will ignore capturing of this info."%filename
			return

		for line in sexConfigFile:
			if(line[0]!="#"):      # Ignore lines that are comments
				if len(line.split())>=2:
					option = line.split()[0]
					value = line.split()[1]
					allSexOptions[option] = value
		
		sexConfigFile.close()
		
		# Filter out only the sextractor parameters that we care about
		for lookup in self.sexOptions:
			self.sexOptions[lookup] = allSexOptions[lookup]
			
		return self.sexOptions
	
		
	def updateRunInfo(self, object):
		self.comment = object['comment']
		self.ra = object['ra']
		self.dec = object['dec']
		self.objectID = object['id']
		self.target = object['target']
		self.num = object['num']
		self.expose = object['expose']

	def writeSelf(self, configData):
		""" Writes itself to a JSON file in the web folder 
		"""
		filename = ultracamutils.addPaths(configData.SITE_PATH, self.runDate)
		filename = ultracamutils.addPaths(filename, self.runID)
		if self.version!='primary': filename+= "_" + str(self.version);
		filename+= "_info.json"
		print "Updating runinfo in file: ", filename
		
		outputObject = {}
		outputObject["target"] = self.target
		outputObject["date"] = self.runDate
		outputObject["runID"] = self.runID
		outputObject["comment"] = self.comment
		outputObject["ra"] = self.ra
		outputObject["dec"] = self.dec
		outputObject["expose"] = self.expose
		outputObject["objectID"] = self.objectID
		outputObject["dataProcessed"] = self.dataProcessed
		outputObject['numWindows'] = self.numWindows
		outputObject['maxExtents'] = self.maxExtents
		
		if hasattr(self, 'sexMagnitude'):
			outputObject['sexMagnitude'] = self.sexMagnitude
			
		if hasattr(self, 'sexOptions'):
			outputObject['sexOptions'] = self.sexOptions
		
		JSONfile = open(filename, 'w')
		json.dump(outputObject, JSONfile)
		JSONfile.close()
		
		
	def __str__(self):
		outStr = "RunID: " + self.runID + "\n"
		outStr+= "Date: " + self.runDate + "\n"
		outStr+= "Target: " + self.target + " RA:" + str(self.ra) + " DEC: " + str(self.dec) + "\n"
		outStr+= "Frames: " + str(self.numFrames) + "\n"
		outStr+= "Expose time: " + str(self.exposeTime) + "\n"
		outStr+= "Mode: " + self.mode + "\n"
		outStr+= "nBlue: " + str(self.nblue) + "\n"
		outStr+= "Comments: " + str(self.comment) + "\n"
		outStr+= "Num: " + str(self.num) + "\n"
		outStr+= "Expose: " + str(self.expose) + "\n"
		
		return outStr
		
	def determineRunClass(self):
		""" Determines the type of run returns as an integer and label  (ntuple)
		0. Unknown
		1. Full field science run (a)... Many objects (>40), large windows (~ full fields), 100-2000 frames
		2. Intermediate science run (b)... Intermediate objects (4-40), intermediate windows (~ half fields), 100-2000 frames
		3. High speed science run (c)...  Few objects (2-4), short exposures, DRIFT mode, small windows (~ quarter size), many frames (2000-100000)
		4. Acquisition run. Science run (1, 2, 3) with few frames <50 and several glitches
		5. Bias run 
		6. Flat field. Full CCD exposed. Start or end of the evening. ~100 exposures
		7. Junk (unclassified)
		"""
		runClassLabels = ["Unknown", "Science A", "Science B", "Science C", "Acquisition", "Bias", "Flat", "Junk"]
		
		if (self.numFrames>2000):
			self.runClass = 1
			return 1, runClassLabels[1]
			
		self.runClass = 0	
		return 0, runClassLabels[0] 

	def setComment(self, comment):
		self.comment = comment


class configObject:
    """ This class stores all of the configuration for the run
	"""
    def __init__(self):
		self.ULTRACAMRAW = "/storage/astro1/phsaap/ultracam/raw_data"
		self.DEBUG = 1
		self.SITE_PATH = "/storage/astro2/phrnaw/ucamsite"
		self.TMP_PATH = "/tmp"
		self.WRITE_FITS = 0
		self.KEEP_TMP_FILES = 0
		self.WRITE_JSON = 0
		self.MOVIE_TMP_PATH = "/tmp/movie"
		self.MINPIXELDISTANCE = 5
		self.FONT = "/usr/share/fonts/truetype/ubuntu-font-family/Ubuntu-B.ttf"
		self.RUNTEMPLATE = "/home/rashley/astro/ucamsite/templates/runxxx.jinja"
		self.WORKINGDIR = "/storage/astro2/phrnaw/ucamsite"
		self.COMPARISON_THRESHOLD = 95.
		self.SEX_MAGNITUDE = "FLUX_AUTO"
		self.FRAME_STACK_BIN = 10
		self.FRAME_MINIMUM = 100
		self.SIGMA_THRESHOLD = 3.0
		self.APERTURE_RADIUS = 5.0
		self.INNER_SKY = 10.0
		self.OUTER_SKY = 15.0
		self.POLY_DEGREE = 5
		self.REF_APERTURES = 4

    def __setitem__(self, item, value):
		setattr(self, item, value)
	
    def __str__(self):
		out = ""
		for key, value in self.__dict__.items():
			out+= str(key) + " --> " + str(value) + "\n"
		return out

class ExposureObject:
	""" This is a class to encapsulate a single exposure
	"""
	def __init__(self):
		self.centroid = (0,0)     	# tuple containing the position of the object in this exposure
		self.counts = 0             	# counts measured by sextractor for this exposure
		self.exposureTime = 0       	# the exposure time for this exposure as given by Tom's ucam library
		self.MJD = 0                	# the calculated MJD for this exposure from Tom's ucam library
		self.FWHM = 0 			# the width of the image, as calculated by sextractor  
		
	def __str__(self):
		return "MJD: " + str(self.MJD) + " counts: " + str(self.counts)
	   
class FrameObject:
	""" This class contains the definitions of the windows for a specific run
	"""
	def __init__(self):
		self.numWindows = 0
		self.windows = []
		self.nxmax = 0
		self.nymax = 0
		self.minX = 0
		self.maxX = 1081
		self.minY = 0
		self.maxY = 1081
		
	def addWindow(self, xoffset, yoffset, xsize, ysize):
		window = WindowObject(xoffset, yoffset, xsize, ysize)
		self.windows.append(window)
		self.numWindows+= 1
	
	def calcMaxExtents(self):
		self.minX = 1085
		self.maxX = 0
		self.minY = 1085
		self.maxY = 0
		for w in self.windows:
			if w.xll<self.minX: self.minX = w.xll;
			if (w.xll+w.xsize)>self.maxX: self.maxX = w.xll+w.xsize - 1;
			if w.yll<self.minY: self.minY = w.yll;
			if (w.yll+w.ysize)>self.maxY: self.maxY = w.yll+w.ysize - 1;
		#print "Max extents: x: (%d - %d)  -  y: (%d - %d)"%(self.minX, self.maxX, self.minY, self.maxY)
		
	def getWindow(self, index):
		return self.windows[index]
		
	def getMaxExtents(self):
		return (self.minX, self.maxX, self.minY, self.maxY)
		
	def __str__(self):
		out = "NumWindows: " + str(self.numWindows) + " dimensions (" + str(self.nxmax) + ", " + str(self.nymax) + ")"
		out += "\n"
		for i in self.windows:
			out+= str(i) + "\n"
		return out

class stackedImage:
	def __init__(self):
		self.numWindows = 0
		self.windows = []
		self.windowImages = []
		self.exposuresCount = [] 
		
	def addWindow(self, windowData, windowImage):
		self.windowImages.append(windowImage)
		self.windows.append(windowData)
		self.exposuresCount.append(0)
		self.numWindows+= 1
		
	def getWindow(self, index):
		return self.windowImages[index]
		
	def addNewData(self, windowData, index):
		originalImage = self.windowImages[index]
		newImage = numpy.add(originalImage, windowData)
		self.windowImages[index] = newImage
		self.exposuresCount[index] = self.exposuresCount[index] + 1
		

class debugObject:
	def __init__(self, debugLevel):
		self.timeStamp = 0
		self.debugText = ""
		self.debugLevel = debugLevel
		self.timeLog = False
	
	def __str__(self):
		out = ""
		if(self.timeLog): 
			self.timeStamp = datetime.datetime.now().strftime("[%H:%M:%S]")
			out+= str(self.timeStamp) + " "
		out+= str(self.debugText)
		return out
		
	def setLevel(self, newLevel):
		self.debugLevel = newLevel
		
	def toggleTimeLog(self):
		if self.timeLog: 
			self.timeLog = False
		else:
			self.timeLog = True

	def error(self, debugText):
		self.debugText = "ERROR: " + debugText
		print str(self)
		
	def write(self, debugText, level = 3):
		self.debugText = debugText
		if (int(self.debugLevel)>=int(level)): print str(self)

				

class WindowObject:
	""" This class contains the definitions of a window within a run
	"""
	def __init__(self, xll, yll, xsize, ysize):
		self.xll = xll
		self.yll = yll
		self.xsize = xsize
		self.ysize = ysize
		
	def __str__(self):
		out = "lower left: (" + str(self.xll) + ", " + str(self.yll) + ") size (" + str(self.xsize) + ", " + str(self.ysize) + ")"
		return out

class combined3ColourObject:
	""" This is a class to encapsulate a single observed object with all three channel data combined
	"""
	colours = ['r', 'g', 'b']
	colourNames = ['Red', 'Green', 'Blue']
	
	def __init__(self, id):
		self.id = id
		self.colourIDs = {'r': -1, 'g': -1, 'b': -1}
		
	def setColourID(self, colour, ID):
		self.colourIDs[colour] = ID
		
	def getColourID(self, colourIndex):
		colourID = self.colourIDs[colourIndex]
		return colourID
		
	def __str__(self):
		outString = "ID: " + str(self.id) + " \n"
		
		for n, c in enumerate(combined3ColourObject.colours):
			outString+= "[%s: %s] \n"%(combined3ColourObject.colourNames[n], self.colourIDs[c])

		return outString
		
	def summaryString(self):
		outString = "ID: " + str(self.id) + " [R:%3d G:%3d B:%3d]"%(self.colourIDs['r'], self.colourIDs['g'], self.colourIDs['b']) 
		return outString
		
	def toJSON(self):
		outObject = {'id':0, 'r': -1, 'g':-1, 'b':-1}
		outObject['id'] = self.id
		for c in combined3ColourObject.colours: 
			outObject[c] = self.colourIDs[c]
		return json.dumps(outObject)

class FrameCatalogObject:
	""" This is a class to keep track of the catalogs we are going to follow. We will keep a catalog for each window. 
	    A catalog is an array of the positions (x, y) of all  the objects in the window.
	"""	
	def __init__(self, numWindows):
		self.numWindows = numWindows
		self.catalogs = []
		for i in range(numWindows):
			cat = []
			self.catalogs.append(cat)
		
	def setCatalog(self, windowIndex, catalog):
		self.catalogs[windowIndex] = catalog
		
	def getCatalog(self, windowIndex):
		return self.catalogs[windowIndex]

class ObservedObject:
	""" This is a class to encapsulate a single observed object in a single channel (colour), usually a star
	"""
	def __init__(self, id):
		self.id = id                # ID for this object, unique to this run
		self.CCDchannel = 0         # which channel am I on? r = 1, g = 2, b = 3, undefined = 0 
		self.windowIndex = 0           # which window am I on? ID from 0 to n 
		self.distanceThreshold = 3.0
		self.numExposures = 0       # The number of exposures in this run for this object
		self.exposures = []         # An array containing numExposures exposure objects
		self.currentPosition = (0,0) # The object's last known x, y position
		self.lastCounts = 0
		self.meanFlux = 0
		self.meanFWHM = 0
		self.ra = 0
		self.dec = 0
		
	def addExposure(self, x, y, counts, FWHM, MJD, exposureTime):
		""" Adds a new exposure object to this object
		"""
		exposure = ExposureObject()
		exposure.centroid = (x, y)
		exposure.counts = counts
		exposure.exposureTime = exposureTime
		exposure.MJD = MJD
		exposure.FWHM = FWHM
		self.exposures.append(exposure)
		self.numExposures+= 1
		self.currentPosition = (x, y)
		
	def resetExposures(self, exposures):
		self.exposures = exposures
		self.numExposures = len(exposures)
		
	def removeExposure(self, index):
		""" Removes an expsure by index. This is used to remove frames with bad timings (like at the start of a run in DRIFT mode)
		"""
		self.exposures.pop(index)
		print "Removing exposure at frame: ", index
		self.numExposures-= 1
		
	def calculateMeanFlux(self):
		meanFlux = 0
		for i in range(self.numExposures):
			flux = self.exposures[i].counts
			meanFlux += flux
		meanFlux /= self.numExposures
		
		self.meanFlux = meanFlux
		
		return self.meanFlux

	def calculateMeanFWHM(self):
		meanFWHM = 0
		for i in range(self.numExposures):
			fwhm = self.exposures[i].FWHM
			meanFWHM += fwhm
		meanFWHM /= self.numExposures
		
		self.meanFWHM = meanFWHM
		
		return self.meanFWHM


	def calculateMeanPosition(self):
		(mx, my) = (0, 0)
		for i in range(self.numExposures):
			(x, y) = self.exposures[i].centroid
			mx += x
			my += y
		mx /= self.numExposures
		my /= self.numExposures
		
		self.meanPosition = (mx, my)
		
		return self.meanPosition
		
	def setWindowIndex(self, windowIndex):
		self.windowIndex = windowIndex
		
	def setWorldPosition(self, ra, dec):
		self.ra = ra
		self.dec = dec
		
	def addExposureByObject(self, newValue, MJD):
		""" Adds a new exposure object to this object
		"""
		exposure = ExposureObject()
		exposure.centroid = ( newValue['absX'], newValue['absY'] ) 
		exposure.counts = newValue['counts']
		exposure.exposureTime = 0
		exposure.MJD = MJD
		exposure.FWHM = newValue['radius']
		self.exposures.append(exposure)
		self.numExposures+= 1
		self.currentPosition = (newValue['absX'], newValue['absY'] )
		
	def isDistanceMatch(self, object):
		""" Returns -1 if object is not a match, or the distance if it is closer than the distanceThreshold
		"""
		xo = object['absX'] - self.currentPosition[0]
		yo = object['absY'] - self.currentPosition[1]
		distance = math.sqrt(xo*xo + yo*yo)
		if (distance>self.distanceThreshold): return -1;
		return distance;
				

	def isCosmicRay(self):
		""" Returns true if this object thinks it is a cosmic ray, based on the simple test that it only appears in 1 frame
		"""
		if self.numExposures == 1: return True;
		return False
		
	def isInCircle(self, x, y, radius):
		""" Returns -1 if object is not in the circle, or the distance if it is inside the circle
		"""
		xo = x - self.currentPosition[0]
		yo = y - self.currentPosition[1]
		distance = math.sqrt(xo*xo + yo*yo)
		if (distance>radius): return -1;
		return distance;	
		
	def __repr__(self):
		self.lastCounts = self.exposures[-1].counts
		return repr((self.id, self.lastCounts))
		
	def __str__(self):
		""" Returns a nicely formatted description of this object (for debug purposes)
		"""
		out = ""
		out += "ID: " + str(self.id) + " (%.f,%.f)[%.2f]"%(self.currentPosition[0], self.currentPosition[1], self.meanFWHM)\
		        + " frames: " + str(self.numExposures) + \
				" mean counts: %.2f"%(self.meanFlux)
		return out
		
	def toJSON(self):
		testObject = {'id': 0, 'x':0, 'y':0, 'data':[]}
		testObject['id'] = self.id
		testObject['x'] = self.currentPosition[0]
		testObject['y'] = self.currentPosition[1]
		exposureDataArray = []
		for c in self.exposures:
			exposureData = (c.MJD, float(c.counts), c.centroid[0], c.centroid[1], float(c.FWHM))
			exposureDataArray.append(exposureData)
		testObject['data'] = exposureDataArray
		return json.dumps(testObject)
		
	def getDeltaXY(self):
		""" Returns the delta XY between the two most recent frames for this object or (0,0) if the object only has one frame
		"""
		if (self.numExposures<2): return (0, 0);
		deltaXY = (0, 0) 
		thisXY = self.exposures[-1].centroid
		previousXY = self.exposures[-2].centroid
		print thisXY, previousXY,
		deltaXY = (thisXY[0]-previousXY[0], thisXY[1]-previousXY[1]) 
		print deltaXY
		return deltaXY
		
	def getData(self):
		""" Returns a tuple of arrays with the MJD and counts for this object. Useful for plotting
		"""
		MJDArray = []
		CountsArray = []
		for i in self.exposures:
			MJDArray.append(i.MJD)
			CountsArray.append(i.counts)
		return (MJDArray, CountsArray) 
		
	def getMJDs(self):
		""" Returns an array with the MJD for this object. Useful for plotting
		"""
		MJDArray = []
		for i in self.exposures:
			MJDArray.append(i.MJD)
		return MJDArray 

	def getLastCounts(self):
		""" Returns an integer containing the last read counts figure for this object
		"""
		return self.exposures[-1].counts
		
	def getCountsForMJD(self, MJD):
		""" Returns counts for a certain MJD, returns 0 if not found
		"""
		counts = 0
		for e in self.exposures:
			 if e.MJD==MJD: counts = e.counts
		return counts
  
	def getCounts(self):
		""" Returns an array with the counts for this object. Useful for plotting
		"""
		CountsArray = []
		for i in self.exposures:
			CountsArray.append(i.counts)
		return CountsArray 

