import numpy, json
import xml.etree.ElementTree as ElementTree
from scipy.ndimage.filters import gaussian_filter
import ultracamutils
from astropy.table import Table


class sourceList:
	""" This is a class of a list of sources
	"""
	
	def __init__(self):
		self.sources = []

	def addSourceByPosition(self, position):
		newSource = source(position)
		self.sources.append(newSource)

	def addSource(self, source):
		if (source.id==0):
			# Generate a unique id for this new source
			newID = ultracamutils.getUniqueID(self.sources)
			source.id = newID
		self.sources.append(source)
		
	def getSources(self):
		return self.sources
		
	def getSourcesByWindow(self, index):
		sourcesInWindow = []
		for s in self.sources:
			if s.windowIndex == index:
				sourcesInWindow.append(s)
		return sourcesInWindow
		
	def writeToCSV(self, filename):
		outFile = file(filename, 'w')
		if len(self.sources)<1: return
		
		headerLine = "ID, window, x, y, abs_x, abs_y, flux, sharpness\n"
		outFile.write(headerLine)
		for s in self.sources:
			outLine = ""
			outLine+= str(s.id) + ", "  + str(s.windowIndex) + ", "
			outLine+= str(s.position[0]) + ", " + str(s.position[1]) + ", "
			outLine+= str(s.position[0] + s.xll) + ", " + str(s.position[1] + s.yll) + ", "
			outLine+= str(s.flux) + ", "
			outLine+= str(s.sharpness)
			outLine+= "\n"
			outFile.write(outLine)
		outFile.close()
		
	def getNumSources(self):
		return len(self.sources)
		
	def loadFromCSV(self, filename):
		try:
			inFile = file(filename, 'r')
		except IOError:
			return False
			
		
		for index, line in enumerate(inFile):
			#print line
			if index==0: continue
			valueList = line.split(",")
			id = int(valueList[0])
			window = int(valueList[1])
			position = ( float(valueList[2]), float(valueList[3]) )
			abs_position = ( float(valueList[4]), float(valueList[5]) )
			sourceObject = source(id, position, window)
			offset = numpy.subtract(abs_position, position)
			sourceObject.latestPosition = position
			sourceObject.setOffsets(offset)
			sourceObject.abs_position = abs_position
			flux = float(valueList[6])
			sharpness = float(valueList[7])
			sourceObject.flux = flux
			sourceObject.latestFlux = flux
			sourceObject.sharpness = sharpness
			self.sources.append(sourceObject)
			
		return True
		
	def keepTopSourcesOnly(self, n):
		self.sortByFlux()
		maxSources = n
		shortenedSourceList = []
		for index, s in enumerate(self.sources):
			if index>=maxSources: break
			self.shortenedSourceList.append(s)	
		print "Reduced number of sources is:", len(self.sources)

		
	def sortByFlux(self):
		sortedSources = sorted(self.sources, key=lambda object: object.flux, reverse = True)
		for index, s in enumerate(sortedSources):
			s.id = index
		self.sources = sortedSources
	
class referenceApertures:
	""" This class contains a list of sources that are the reference apertures for the run....
	"""
	def __init__(self):
		sources = []
		
	def addSource(self, source):
		self.sources.append(source)
		
	def getSources(self):
		return self.sources
		
	def initFromSourceList(self, sourceList, max=4):
		self.sources = []
		sourceList.sortByFlux()
		topPercent = 40.
		maxApertures = max
		
		for index, s in enumerate(sourceList.sources):
			if index>=maxApertures: break
			self.sources.append(s)
		#print "Number of reference apertures is:", len(self.sources)
		
	def calculateFrameCoverage(self, nFrames):
		for s in self.sources:
			numFrames = len(s.positionLog)
			coverage = float(numFrames) / float(nFrames) * 100.
			s.coverage = coverage
			print s, coverage
			
	def sortByCoverage(self):
		sortedSources = sorted(self.sources, key=lambda object: object.coverage, reverse = True)
		self.sources = sortedSources

	def sortByFlux(self):
		sortedSources = sorted(self.sources, key=lambda object: object.flux, reverse = True)
		self.sources = sortedSources

		
	def limitCoverage(self, minCoverage):
		newSources = []
			
		for s in self.sources:
			if s.coverage>=minCoverage: 
				newSources.append(s)
			else:
				print "Removing aperture: ", s
		self.sources = newSources

class source:
	""" This is a definition of a source 
	"""
	def __init__(self, id, position, windowIndex):
		self.id = id
		self.position = position
		self.windowIndex = windowIndex
		self.binningFactor = 0
		self.xll = 0
		self.yll = 0 
		self.positionLog = []
		self.latestPosition = (0, 0)
		self.fluxMeasurements = []
		self.latestFlux = 0
		self.polyFit = {}
		self.polyFit["x"] = {}
		self.polyFit["x"]["order"] = 0
		self.polyFit["y"] = {}
		self.polyFit["y"]["order"] = 0
		self.recentFail = False
		
	def setPolynomial(self, axis, order, parameters):
		if axis!="x" and axis!="y": return False
		self.polyFit[axis]["order"] = order
		self.polyFit[axis]["parameters"] = parameters
		print "Setting polynomial:", axis, order, parameters
		return True
		
		
	def setDAOPhotData(self, sharpness, roundness1, roundness2, npix, sky, peak, flux, mag):
		self.sharpness = sharpness
		self.roundness = (roundness1, roundness2)
		self.npix = npix
		self.sky = sky
		self.peak = peak
		self.flux = flux
		self.mag = mag
		
	def setOffsets(self, offset):
		self.xll, self.yll = offset
		
	def getPositionByFrame(self, frameNumber):
		for p in self.positionLog:
			if p['frameNumber']== frameNumber:
				return p['position']
		return (0, 0)
		
	def getData(self):
		dataPoints = []
		for f in self.fluxMeasurements:
			position = self.getPositionByFrame(f['frameNumber'])
			data = [ f['frameNumber'], f['flux'], position[0], position[1]]
			dataPoints.append(data)
		return dataPoints
		
	def setLatestPosition(self, frameNumber, newPosition, errors = (0, 0)):
		self.latestPosition = newPosition
		logEntry = {}
		logEntry['frameNumber'] = frameNumber
		logEntry['position'] = newPosition
		logEntry['positionError'] = errors
		self.recentFail = False
		self.positionLog.append(logEntry)
		
	def addFluxMeasurement(self, frameNumber, flux):
		self.latestFlux = flux
		logEntry = {}
		logEntry['frameNumber'] = frameNumber
		logEntry['flux'] = flux
		self.fluxMeasurements.append(logEntry)

	def __str__(self):
		outString = ""
		outString+= "ID: " + str(self.id) + " W[%d]"%self.windowIndex + " position (%3.2f, %3.2F)"%(self.latestPosition) + " Flux: %f"%(self.latestFlux)
		return outString 

class aperture:
	""" This is an aperture instance 
	"""
	def __init__(self, id, position, radius):
		self.position = position
		self.radius = radius
		self.hasSky = False
		self.skyInnerRadius = 0
		self.skyOuterRadius = 0
		self.id = id
		
	def pixelArea(self):
		return numpy.pi * self.radius * self.radius
	
	def fluxPerPixel(self, flux):
		return flux / self.pixelArea()	
		
	def setSkyAperture(self, skyInnerRadius, skyOuterRadius):
		self.hasSky = True
		self.skyInnerRadius = skyInnerRadius
		self.skyOuterRadius = skyOuterRadius
		
	def setRadius(self, radius):
		self.radius = radius
		
class window:
	""" This class defines a window for the ULTRASPEC camera
	"""
	def __init__(self):
		self.xbin = 0
		self.ybin = 0
		self.xll = 0
		self.yll = 0
		self.nx = 0
		self.ny = 0
		self.data = None
		self.stackedData = []
		self.BGSubtractedImage = []
		self.sources = None
		self.borderWidth = 10
		
	def setExtents(self, xll, yll, nx, ny):
		self.xll = xll
		self.yll = yll
		self.nx = nx
		self.ny = ny
		
	def setBinning(self, xbin, ybin):
		self.xbin = xbin
		self.ybin = ybin
		
	def setData(self, data):
		self.data = data
		if len(self.stackedData) == 0:
			self.stackedData = data
	
	def setBlankData(self, data):
		dataShape = numpy.shape(data)
		blanks = numpy.zeros(dataShape)
		if len(self.stackedData) == 0:
			self.stackedData = blanks

	def addData(self, data):
		self.stackedData = self.stackedData + data
		self.data = data
		
	def addToStack(self, data):
		self.stackedData = self.stackedData + data
		
	def filterBorderSources(self):
		currentNumSources = len(self.sources)
		if currentNumSources==0:
			return 0
		sourcesToRemove = []
		for index, s in enumerate(self.sources):
			x = s['xcentroid']
			y = s['ycentroid']
			if (x < self.borderWidth) or (x > (self.nx - self.borderWidth)): 
				sourcesToRemove.append(index)
				continue
			if (y < self.borderWidth) or (y > (self.ny - self.borderWidth)): 
				sourcesToRemove.append(index)
				continue
		self.sources.remove_rows(sourcesToRemove)
		return len(self.sources)
		
	def setSources(self, sources):
		self.sources = sources
		
	def setSourcesAvoidBorders(self, sources):
		self.sources = sources
		self.filterBorderSources()
		
	def getSources(self):
		return self.sources
		
	def __str__(self):
		retString = ""
		retString+= "xll, yll, nx, ny (%d, %d, %d, %d)  %dx%d"%(self.xll, self.yll, self.nx, self.ny, self.xbin, self.ybin)
		return retString

class sourceMap:
	""" This class defines a 2D histogram that is going to be used as a 'source map' or a 'heat map of sources'
	"""
	
	def __init__(self, dimensions):
		self.heatMap = numpy.zeros(dimensions)
		self.psize = 1.0
		self.fwhm = 4.0
		self.xsize = dimensions[1]
		self.ysize = dimensions[0]
		self.border = 10
		print "Intialised a source map with dimensions:", dimensions
		
	def updateMap(self, sources):
		for s in sources:
			j, i = int(s[0]), int(s[1])
			if (j<self.border) or (j>self.xsize-self.border): continue
			if (i<self.border) or (i>self.ysize-self.border): continue
			self.heatMap[i][j]+=1 
			
	def getSourceMap(self):
		return self.heatMap
		
	def getSmoothMap(self):
		return gaussian_filter(self.heatMap, self.fwhm/self.psize/2.3548, mode='constant')
		
class runInfo:
	""" This class is used to store the meta-data for the run
	"""
	def __init__(self, runPath):
		self.runPath = runPath
		self.runDate, self.runID = ultracamutils.separateRunNameAndDate(runPath)
		self.comment = ""
		self.ra = 0
		self.dec = 0
		self.target = "unknown"
		self.expose = 0
		self.num = 0
		self.comment = ""
		self.windowInfo = []
		self.powerCycle = False
		self.flags = ''
		
		
	def loadFromJSON(self, JSONFilename):
		JSONfile = open(JSONFilename, "r")
		allObjectsJSON = json.load(JSONfile)
		run = {}
		runNumberStr = self.runID[3:]
		runNumber = int(runNumberStr)
		
		objectNotFound = True
		
		for object in allObjectsJSON:
			date = object['night']
			num = object['num']
			if ((date == self.runDate) & (runNumber == num)):
				objectNotFound = False
				self.comment = object["comment"]
				self.ra = object['ra'] * 15. # Convert the ultra.json RA value to degrees
				self.dec = object['dec']
				self.objectID = object['id']
				self.target = object['target']
				self.num = object['num']
				self.expose = object['expose']
				#print object
				return True 
		return False
            
	def __str__(self):
		description = "%s/%s Target: %s  Comment:%s"%(self.runDate, self.runID, self.target, self.comment)
		return description

	def loadFromXML(self, rawDataDirectory):
		""" This method is called upon if the loadFromJSON method fails to return any data
		"""
		XMLFilename = ultracamutils.addPaths(rawDataDirectory, self.runPath) + ".xml"
		try:
			tree = ElementTree.parse(XMLFilename)
		except IOError as e:
			print "Could not get read the run's XML file at:", XMLFilename
			return False
		
		root = tree.getroot()
		# Read the window sizes and binning, etc
		filesave = root.find('filesave_status')
		instrument = filesave.find('instrument_status')
		detector = instrument.find('detector_status')
		
		powerToDetector = detector.attrib['power']
		if powerToDetector == 'OFF':
			self.powerCycle = True
			return True
		
		parameter = instrument.findall('parameter_status')
		parameters = {}
		for param in parameter:
			parameters[param.attrib['name']] = param.attrib['value']
		if parameters.has_key('COM_TAB_START'):
			self.powerCycle = True
			return True
			
		
		windows = []
		xbin = int(parameters['X_BIN'])
		ybin = int(parameters['Y_BIN'])
		for number in range(1, 5):
			try:
				xStartStr = "X%d_START"%number
				xStart = int(parameters[xStartStr]) 
				yStartStr = "Y%d_START"%number
				yStart = int(parameters[yStartStr]) 
				xSizeStr = "X%d_SIZE"%number
				xSize = int(parameters[xSizeStr]) 
				ySizeStr = "Y%d_SIZE"%number
				ySize = int(parameters[ySizeStr]) 
			except KeyError:
				continue
			if xSize!=0 and ySize!=0:
				# Add this to the list of windows
				windowDict = {}
				windowDict['xll'] = xStart
				windowDict['nx'] = xSize
				windowDict['yll'] = yStart
				windowDict['ny'] = ySize
				windowDict['xbin'] = xbin
				windowDict['ybin'] = ybin
				windows.append(windowDict)
		self.windowInfo = windows
				
		# Get the observation meta-data
		user = root.find('user')
		target = user.find('target').text
		flags = user.find('flags').text
		PI = user.find('PI').text
		ID = user.find('ID').text
		observers = user.find('Observers').text
		try:
			raStr = user.find('RA').text
			decStr = user.find('Dec').text
			ra, dec = ultracamutils.fromSexagesimal(raStr, decStr)
		except AttributeError:
			ra, dec = (0, 0)
		#print "Target:", target
		#print "PI: %s Observers: %s  Programme: %s"%(PI, observers, ID)
		#print "(ra, dec): (%f, %f)"%(ra, dec)
		
		self.ra = ra
		self.dec = dec
		self.target = target
		self.flags = flags
		
		return True
        
