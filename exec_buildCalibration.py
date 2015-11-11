#This scirpt will build the master calibration fields
#==========
#MasterBias
#MasterDark
#MasterFlat

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
if not os.path.isfile(indexFile):
    # Loop through each night and test for image type
    print('\nCategorizing files into "BIAS", "DARK", "FLAT", "OBJECT"\n')
    startTime = os.times().elapsed
    # Begin by initalizing some arrays to store the image classifications
    obsType  = []
    group    = []
    binType  = []
    polAng   = []
    waveBand = []
    lights   = []
    fileCounter = 0
    percentage  = 0
    
    #Loop through each file in the fileList variable
    for file in fileList:
        
        # Classify each file type and binning
        obsType.append(fits.getval(file, 'OBSTYPE'))
        tmpGroup = fits.getval(file, 'OBJECT')
        if len(tmpGroup) < 1:
            tmpGroup = 'blank'
        group.append(tmpGroup)
        polAng.append(fits.getval(file, 'POLPOS'))
        waveBand.append(fits.getval(file, 'FILTNME3'))
        
        # Try to grab the lights on/off status from the comments
        try:
            comments = fits.getval(file, 'COMMENT')
            for comment in comments:
                valueWritten = False
                if 'lights on' in comment.lower():
                    lights.append('lights on')
                    valueWritten = True                    
                elif 'lights off' in comment.lower():
                    lights.append('lights off')
                    valueWritten = True
            
            # If no comments on the lights were found, write 'N/A'
            if not valueWritten:
                lights.append('N/A')
        
        # If there are no comments, then the lights are not relevant
        except KeyError:
            lights.append('N/A')

        binTest  = fits.getval(file, 'CRDELT*')
        if binTest[0] == binTest[1]:
            binType.append(int(binTest[0]))
    
        # Count the files completed and print update progress message
        fileCounter += 1
        percentage1  = np.floor(fileCounter/len(fileList)*100)
    #    if (percentage1 % 5) == 0 and percentage1 != percentage:
        if percentage1 != percentage:
            print('completed {0:3g}%'.format(percentage1), end="\r")
        percentage = percentage1
        
        if len(lights) != fileCounter:
            pdb.set_trace()
     
    endTime = os.times().elapsed
    print('\nFile processing completed in {0:g} seconds'.format(endTime -startTime))
    
    # Write the file index to disk
    fileIndex = Table([fileList, obsType, waveBand, target, polAng, binType, lights],
                      names = ['Filename', 'Data',  'Waveband', 'Group',
                      'Polaroid Angle', 'Binning', 'Lights'])
    ascii.write(fileIndex, indexFile)
else:
    # Read the fileIndex back in as an astropy Table
    print('\nReading file index from disk')
    fileIndex = ascii.read(indexFile)#, guess=False, delimiter=' ')

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

    # Test if there is a MasterBias image for this binnng level.
    if os.path.isfile(filename):
        # Read in the file if it exists
        print('\nLoading file into masterBias list')
        print(filename)
        masterBiases.update({thisBin:
                             Bias(filename)})
    else:
        print('\nProcessing biases with ({0:g}x{0:g}) binning'.format(thisBin))
        
        # Construct the key name for this bias image
        biasKey = '{0:g}'.format(thisBin)
        
        # Select the bias images with the correct binning level
        biasImgFiles = fileIndex['Filename'][biasBool & (binType == thisBin)]
        
        # Loop through each of the files and add them to the biasImgList list
        biasImgList  = []
        for file in biasImgFiles:
            biasImgList.append(Bias(file))
        
        #=========================================================================
        # **************************** OVERSCAN **********************************
        # Use overscan regions from the bias images and use them to figure 
        # out what overscan polynomial should be fit to both PRESCAN and 
        # POSTSCAN regions.
        #=========================================================================
        # Find the best fitting polynomial to the prescan and postscan regions.
        overscanPolyDegree = Bias.overscan_polynomial(biasImgList, overscanPos)
        binPolyDegrees.append(overscanPolyDegree)
        print('\nOverscan region fit with {0:g}-order polynomial'.
          format(overscanPolyDegree))
        print('\nProceeding to remove overscans from the bias images')
        for bias in biasImgList:
            bias.overscan_correction(overscanPos, sciencePos,
                                     overscanPolyDegree)
        
        # Perform the Bias.master_bias() method to compute the master bias map
        #
        #
        # TODO write a 'build master bias header' component
        # of the 'master_bias' method.
        #
        # Generate masterBias Image object to save to disk
        masterBias = Bias()
        masterBias.arr = Bias.master_bias(biasImgList)
        masterBias.header = biasImgList[0].header.copy()
        masterBiases.update({thisBin: masterBias})
        
        # Write masterBias object to disk
        fits.writeto(filename,
                     (masterBiases[thisBin].arr).astype(np.float32),
                     biasImgList[0].header,
                     clobber = True)
        
        print('\nThe mean bias level is {0:g} counts\n'.
          format(np.mean(masterBiases[thisBin].arr)))
        
        # Do a quick cleanup to make sure that memory survives
        del biasImgList
        del masterBias

# Check if the polynomial degrees for overscan correction have been saved
if not os.path.isfile('overscanPolynomials.dat'):
    overscanPolyDegrees = Table([uniqBins, binPolyDegrees],
                                 names = ['Binning', 'Polynomial'])
    ascii.write(overscanPolyDegrees, 'overscanPolynomials.dat')
else:
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
    # Test if darks have been computed for this binning level.
    if os.path.isfile(filename):
        print('\nLoading file into masterDark list')
        print(filename)
        masterDarks.update({thisBin:
                            Dark(filename)})
    else:
        print('\nProcessing darks with ({0:g}x{0:g}) binning'.format(thisBin))
        
        # Select the bias images with the correct binning level
        darkImgFiles = fileIndex['Filename'][darkBool & (binType == thisBin)]
        
        print('\nProceeding to remove overscans and bias from the dark images')
        
        # Select the correct polynomial order
        polyInd            = (overscanPolyDegrees['Binning'] == thisBin)
        overscanPolyDegree = int((overscanPolyDegrees[polyInd])['Polynomial'])
        
        # Loop through each of the files and add them to the biasImgList list
        darkImgList  = []
        for file in darkImgFiles:
            thisDark = Dark(file)
            thisDark.overscan_correction(overscanPos, sciencePos,
                                         overscanPolyDegree)
            thisDark.arr = thisDark.arr - masterBiases[thisBin].arr
            darkImgList.append(thisDark)
    
        # Generate a Dark(Image) object to store the final dark current map
        darkCurrent = Dark()
        # Perform the Dark.dark_current() method to compute the dark current map
        darkCurrent.arr = Dark.dark_current(darkImgList)
        darkCurrent.header = darkImgList[0].header
        masterDarks.update({thisBin: darkCurrent})
        print('\nThe mean dark current is {0:g} counts/second'.
          format(np.mean(masterDarks[thisBin].arr)))
        
        # Write darkCurrent Image object to disk
        outFile = calibrationDir + 'MasterDark{0:g}.fits'.format(thisBin)
        fits.writeto(outFile,
                     masterDarks[thisBin].arr.astype(np.float32),
                     darkImgList[0].header,
                     clobber = True)
        
        # Do a quick cleanup to make sure that memory survives
        del darkImgList
        del darkCurrent

#==============================================================================
# ***************************** FLATS *****************************************
# Setup the paths to the flat images and compute the flat map
#==============================================================================
# TODO loop through each waveband (exterior loop)
# TODO compute average lights on flat
# TODO compute average lights off flat
# TODO define magic_method "sub" for image class
# TODO compute masterFlat = lightsOn - lightsOff

# Find the number of unique wavebands used in flats
uniqBands = np.unique(waveBand[flatBool])

# Create an empty dictionary to store the masterFlats,
# keyed to each band/binning pair
# TODO change ALL "master" lists to "master" dictionaries
# TODO make sure all image lists get deleted from memory
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
            
            # Test if the file exists
            if os.path.isfile(filename):
                # Read in the file if it exists
                print('\nLoading file into masterFlat list')
                print(filename)
                masterFlats.update({keyname:
                                    Flat(filename)})
            else:
                # Create the file if it does not exist
                print('\nProcessing flats for')
                print('band    = {0:s}'.format(thisBand))
                print('polAng  = {0:g}'.format(thisAng))
                print('binning = ({0:g}x{0:g})'.format(thisBin))
                
                flatOnImgFiles = fileIndex['Filename'][flatBool &
                                                       (polAng == thisAng) &
                                                       (binType == thisBin) &
                                                       (lights == 'lights on')]
                
                flatOffImgFiles = fileIndex['Filename'][flatBool &
                                                        (polAng == thisAng) &
                                                        (binType == thisBin) &
                                                        (lights == 'lights off')]
                
                print('\nProceeding to remove overscans and bias from the flat images')

                # Select the correct polynomial order for this binning 
                polyInd            = (overscanPolyDegrees['Binning'] == thisBin)
                overscanPolyDegree = int((overscanPolyDegrees[polyInd])['Polynomial'])
                
                # Loop through each of the files and add them to the biasImgList list
                flatOnImgList  = []
                for file in flatOnImgFiles:
                    thisFlat = Flat(file)
                    thisFlat.overscan_correction(overscanPos, sciencePos,
                                                 overscanPolyDegree)
                    # TODO define a subtraction method for image array
                    # test for if it's another image or a similarly sized array
                    thisFlat.arr = thisFlat.arr - masterBiases[thisBin].arr
                    flatOnImgList.append(thisFlat)
                
                flatOffImgList = []
                for file in flatOffImgFiles:
                    thisFlat = Flat(file)
                    thisFlat.overscan_correction(overscanPos, sciencePos,
                                                 overscanPolyDegree)
                    thisFlat.arr = thisFlat.arr - masterBiases[thisBin].arr
                    flatOffImgList.append(thisFlat)
                
                # Compute AVERAGE flat 'lights on' image
                flatImg = Flat.master_flat(flatOnImgList)

                # Compute an off flat to subtract
                if len(flatOffImgList) > 0:
                    flatOffImg = Flat.master_flat(flatOffImgList)
                    flatImg   -= flatOffImg
                
                # Normalize the flat field by the mode
                # Compute the positions of the non-overscan region
                sciencePos1 = sciencePos/thisBin
                flatModeImg = flatImg[sciencePos1[0][1]:sciencePos1[1][1], \
                                      sciencePos1[0][0]:sciencePos1[1][0]]
                
                # Compute the number of bins that will be needed to find mode
                flatBins = np.ceil(0.1*(np.max(flatModeImg) -
                                         np.min(flatModeImg)))
                
                # Generate a histogram of the flat field
                hist, flatBins = np.histogram(flatModeImg.flatten(), flatBins)

                # Locate the histogram maximum
                maxInd = np.where(hist == np.max(hist))
                if len(maxInd) == 1:
                    maxInd = maxInd[0]
                else:
                    print('The mode is ambiguous')
                    pdb.set_trace()
                
                # Estimate flatMode from histogram maximum
                flatMode = np.mean(flatBins[maxInd:maxInd+2])
                
                # Grab the data within 100 counts of the estimated mode
                # Use kernel density estimator to find a more accurate mode
                print('Evaluating gaussian kernel density estimator')
                flatModeData = flatModeImg[np.where(
                  np.abs(flatModeImg - flatMode) < 100)]
                kernel = stats.gaussian_kde(flatModeData)
                
                # Establish bins for evaluating the resultant
                xmin     = np.floor(flatMode - 30)
                xmax     = np.ceil(flatMode + 10)
                flatBins = np.mgrid[xmin:xmax:0.5]
                density  = kernel.evaluate(flatBins)
                flatMode = flatBins[np.where(density == density.max())]
                
                # Normalize flatImg
                flatImg /= flatMode
                
                # TODO replace zeros with average of surrounding non-zero values
                # Replace zero values in overscan regions
                # to prevent division problems later.
                zeroInds = np.where(flatImg == 0)
                flatImg[zeroInds] = 1

                masterFlats.update({keyname:flatImg})
                
                # Write normalized flat to disk
                fits.writeto(filename,
                             (masterFlats[keyname]).astype(np.float32),
                             flatOnImgList[0].header,
                             clobber = True)

                # Do a quick cleanup to make sure that memory survives
                del flatOnImgList
                del flatOffImgList

# Let the user know everything completed
print('\n..........')
print('Finished producing master calibration fields')
