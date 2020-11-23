from ij import IJ, ImagePlus, ImageStack, CompositeImage
from ij.process import FloatProcessor
from sys.float_info import max as MAX_FLOAT

#imp = IJ.getImage()
fname = "M:/tnw/bn/dm/Shared/Lukas/BEP/Experiments/WKS024/B02/B02_wts.tif"
#imp = IJ.openImage(fname)
imp = IJ.getImage()
ip = imp.getProcessor()
pixels = ip.getPixels()

'''
maximum = reduce(max, pixels)
network = {}
print maximum
'''
width = imp.width
print len(pixels)
maximum = 0


for iC in range(len(pixels)):
	xC = iC % width
	
	#if xC == 0:
	#	continue
	
	yC = iC / width

	iW = iC - 1
	iN = iC - width
	iNW = iC - 1 - width

	center = pixels[iC]
	
	if yC > maximum:
		maximum = yC

print maximum

'''

for i in range(1,int(maximum+1)):
	network[i] = []

for y in xrange(1,ip.getHeight()):
	for x in xrange(1,ip.getWidth()):
		center = ip.get(x,y)

		if center > maximum:
			print center, ' ', x, ' ', y

		if center > 0:
			north = ip.get(x,y-1)
			west = ip.get(x-1,y)
			northwest = ip.get(x-1,y-1)
			
			if north > 0 and center != north and north not in network[center]:
				network[center].append(north)
				network[north].append(center)
			
			if west > 0 and center != west and west not in network[center]:
				network[center].append(west)
				network[west].append(center)

			if northwest > 0 and center !=northwest and northwest not in network[center]:
				network[center].append(northwest)
				network[northwest].append(center)

print "Network is ready!!"
'''



'''
#imp = IJ.openImage("https://imagej.nih.gov/ij/images/flybrain.zip")
IJ.run(imp, "Properties...", "channels=3 slices=1 frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");

stack = imp.getImageStack()
print stack.getSize()

stack2 = ImageStack(imp.width, imp.height)

for i in range(1, imp.getNSlices()+1):
	cp = stack.getProcessor(i)
	red = cp.toFloat(0,None)
	blue = cp.toFloat(2,None)
	stack2.addSlice(None,red)
	stack2.addSlice(None,blue)
	print i

imp2 = ImagePlus("32-bit 2-channel composite", stack2)
imp2.setCalibration(imp.getCalibration().copy())

nChannels = 2
nSlices = stack.getSize()
nFrames = 1

imp2.setDimensions(nChannels, nSlices, nFrames)
comp = CompositeImage(imp2, CompositeImage.COLOR)
comp.show()
'''
'''
print 'Number of slices = ', imp.getNSlices()

ip = imp.getProcessor().convertToFloat()

ip0 = ip.toFloat(0,None)
pixels0 = ip0.getPixels()

ip2 = ip.toFloat(2,None)
pixels2 = ip2.getPixels()

minimum = reduce(max, pixels2)
print 'Minimum is ',minimum

'''
'''
# get minimum pixel
minimum = MAX_FLOAT

xc = 0
yc = 0
for pix in pixels:
	if pix < minimum:
		minimum = pix

print 'Minimum is ',minimum

# subtract minimum pixel
pixels2 = map(lambda x: x - minimum, pixels)
ip2 = FloatProcessor(ip.width, ip.height, pixels2, None)
imp2 = ImagePlus(imp.title, ip2)

imp2.show()

#imp.show()
'''