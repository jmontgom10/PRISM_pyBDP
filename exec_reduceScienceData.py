#Apply the calibration fields to the RAW science data to produce
#reduced science data.

#Import whatever modules will be used
import os
import sys
import pdb
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits, ascii
from scipy import stats
from pyBDP import Image, Bias, Dark, Flat


#Setup the path delimeter for this operating system
delim = os.path.sep

#==============================================================================
# *********************** CUSTOM USER CODE ************************************
# this is where the user specifies where the raw data is stored
# and some of the subdirectory structure to find the actual .FITS images
#==============================================================================
# This is the location of the raw data for the observing run
rawPath = '/home/jordan/ThesisData/PRISM_Data/Raw_data/'
# This is a list of strings containing the subdirectory structure for each night
subDirs = ['20150117/PRISM_Images/', \
           '20150118/PRISM_Images/', \
           '20150119/PRISM_Images/']

# This is a list of strings containing the prefix of the filenames
#filePrefix = ['20150117', '20150118', '20150119']

# This line prepends the rawPath variable to each element in the subDirs list
rawDirs = [rawPath + subDir for subDir in subDirs]

calibrationDir = '/home/jordan/ThesisData/PRISM_Data/Calibration/'

#Loop through each night and build a list of all the files in observing run
fileList = []
for night in subDirs:
    nightPath = rawPath + night
    for file in os.listdir(nightPath):
        fileList.extend([os.path.join(nightPath, file)])

#Sort the fileList
fileNums = [''.join((file.split(delim).pop().split('.'))[0:2]) for file in fileList]
sortInds = np.argsort(np.array(fileNums, dtype = np.int))
fileList = [fileList[ind] for ind in sortInds]

# These are the overscan regions for all PRISM frames at 1x1 binning
#                       ((x1,y1), (x2, y2))
overscanPos = np.array([[2110, 8],[2177, 2059]], dtype = np.int32)
sciencePos  = np.array([[70,  32],[2070, 2032]], dtype = np.int32)

#==============================================================================
# ***************************** INDEX *****************************************
# Build an index of the file type and binning, and write it to disk
#==============================================================================
# Check if a file index already exists... if it does then just read it in
indexFile = 'fileIndex.dat'

# Read the fileIndex back in as an astropy Table
print('\nReading file index from disk')
fileIndex = ascii.read(indexFile)

biasBool = (fileIndex['Data'] == 'BIAS')
darkBool = (fileIndex['Data'] == 'DARK')
flatBool = (fileIndex['Data'] == 'FLAT')
sciBool  = (fileIndex['Data'] == 'OBJECT')
waveBand = fileIndex['Waveband']
polAng   = fileIndex['Polaroid Angle']
binType  = fileIndex['Binning']
lights   = fileIndex['Lights']

#==============================================================================
# ***************************** BIAS *****************************************
# Setup the paths to the bias images and compute the bias map
#==============================================================================

# Find the number of unique binnings used in biases
uniqBins = np.unique(binType[biasBool])

masterBiases = {}
binPolyDegrees = []
for thisBin in uniqBins:
    # Construct the filename for this bias.
    filename = calibrationDir + 'MasterBias{0:g}.fits'.format(thisBin)

    # Read in the masterBias file
    print('\nLoading file into masterBias list')
    print(filename)
    masterBiases.update({thisBin:
                         Bias(filename)})

# Read in the overscan polynomial degrees
overscanPolyDegrees = ascii.read('overscanPolynomials.dat')


#==============================================================================
# ***************************** DARKS *****************************************
# Setup the paths to the dark images and compute the dark current map
#==============================================================================

# Find the number of unique binnings used in darks
uniqBins = np.unique(binType[darkBool])

masterDarks = {}

for thisBin in uniqBins:
    # Construct filename for this dark
    filename = calibrationDir + 'MasterDark{0:g}.fits'.format(thisBin)
    
    # Read in the file
    print('\nLoading file into masterDark list')
    print(filename)
    masterDarks.update({thisBin:
                        Dark(filename)})

#==============================================================================
# ***************************** FLATS *****************************************
# Setup the paths to the flat images and compute the flat map
#==============================================================================

# Find the number of unique wavebands used in flats
uniqBands = np.unique(waveBand[flatBool])

masterFlats = {}

#Loop through each waveband
for thisBand in uniqBands:
    # Compute the unique values for the polaroid rotation angle
    uniqPolAngs = np.unique(polAng[flatBool &
                            (waveBand == thisBand)])
    for thisAng in uniqPolAngs:
        # Compute the unique values for the binning level
        uniqBins = np.unique(binType[flatBool &
                             (waveBand == thisBand) &
                             (polAng == thisAng)])
        # TODO replace binNumber syntax with dictionary syntax
        for thisBin in uniqBins:
            # Construct the keyname and filename for this image
            keyname  = '{0:s}_{1:g}_{2:g}'.format(thisBand, thisAng, thisBin)
            filename = calibrationDir + 'MasterFlat{0:s}_{1:g}_{2:g}.fits'.format(
              thisBand, thisAng, thisBin)
            
            # Read in the masterFlat file
            print('\nLoading file into masterFlat list')
            print(filename)
            masterFlats.update({keyname:
                                Flat(filename)})

#==============================================================================
# **************************** SCIENCE ****************************************
# Setup the paths to the
#==============================================================================
reducedDir      = '/home/jordan/ThesisData/PRISM_Data/Reduced_data'
scienceImgFiles = fileIndex['Filename'][sciBool]
print('\nBeginning to reduce science data.')
for file in scienceImgFiles:
    # Read in the science image from disk
    thisImg = Image(file)
    # Perform the overscan correction appropriate for this binning
    polyInd            = (overscanPolyDegrees['Binning'] == thisImg.binning)    
    overscanPolyDegree = int((overscanPolyDegrees[polyInd])['Polynomial'])
    thisImg.overscan_correction(overscanPos, sciencePos,
                                overscanPolyDegree)
    
    # Correct the 2D bias structure
    thisImg.arr = thisImg.arr - masterBiases[thisImg.binning].arr
    
    # Subtract dark current (not necessary in this case!)
#    avg, sig = np.median(thisImg.arr), thisImg.arr.std()
#    plt.imshow(thisImg.arr, vmin = avg-0.5*sig, vmax = avg+0.5*sig)
#    plt.show()
    
    # Corect the flat-field structure
    thisBand = thisImg.header['FILTNME3']
    thisAng  = thisImg.header['POLPOS']
    keyname  = '{0:s}_{1:g}_{2:g}'.format(thisBand, thisAng, thisImg.binning)
    thisImg.arr = thisImg.arr/masterFlats[keyname].arr
    
    outFile = reducedDir + delim + file.split(delim).pop()
    
    # TODO should this be changed to use the Image.write() method?
    fits.writeto(outFile,
                 thisImg.arr,
                 thisImg.header,
                 clobber = True)

# Let the user know everything is finished
print('\n..........')
print('Finished reducing science data.')
#    avg, sig = np.median(thisImg.arr), thisImg.arr.std()
#    plt.imshow(thisImg.arr, vmin = avg-0.5*sig, vmax = avg+0.5*sig)
#    plt.show()
