#This scirpt will build the master calibration fields
#==========
#MasterBias
#MasterDark
#MasterFlat

#Import whatever modules will be used
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.table import Column
from astropy.io import fits, ascii
from scipy import stats

# Import AstroImage
import astroimage as ai

# Add the header handler to the BaseImage class
from PRISM_header_handler import PRISM_header_handler
ai.BaseImage.set_header_handler(PRISM_header_handler)

#Setup the path delimeter for this operating system
delim = os.path.sep

#==============================================================================
# *********************** CUSTOM USER CODE ************************************
# this is where the user specifies where the raw data is stored
# and some of the subdirectory structure to find the actual .FITS images
#==============================================================================
# This is the location of the raw data for the observing run
rawDir = 'C:\\Users\\Jordan\\FITS Data\\PRISM_data\\raw_data\\201612\\'

# Define the path to the parent directory for all pyBDP products
pyBDP_data = 'C:\\Users\\Jordan\\FITS_data\\PRISM_data\\pyBDP_data\\201612\\'

# Define the directory into which the average calibration images will be placed
calibrationDir = os.path.join(pyBDP_data, 'master_calibration_images')
if (not os.path.isdir(calibrationDir)):
    os.mkdir(calibrationDir, 0o755)

# Read the fileIndex back in as an astropy Table
print('\nReading file index from disk')
indexFile = os.path.join(pyBDP_data, 'rawFileIndex.csv')
fileIndex = Table.read(indexFile, format='csv')

# Locate where the bias, dark, flat, and science images are in the index
biasBool = (fileIndex['Data'] == 'BIAS')
darkBool = (fileIndex['Data'] == 'DARK')
flatBool = (fileIndex['Data'] == 'FLAT')
sciBool  = (fileIndex['Data'] == 'OBJECT')

# Extract lists of the waveband, polaroid angle, binning, and lights-on/off
waveBand = fileIndex['Waveband']
polAng   = fileIndex['Polaroid Angle']
binning  = fileIndex['Binning']
lights   = fileIndex['Lights']

#==============================================================================
# ***************************** BIAS *****************************************
# Setup the paths to the bias images and compute the bias map
#==============================================================================

# Find the number of unique binnings used in biases
uniqBins = np.unique(binning).astype(int)

masterBiasDict = {}
for thisBin in uniqBins:
    # Construct the filename for this bias.
    masterBiasFilename = 'MasterBias{0:g}.fits'.format(thisBin)
    masterBiasFilename = os.path.join(calibrationDir, masterBiasFilename)

    # Test if there is a MasterBias image for this binnng level.
    if os.path.isfile(masterBiasFilename):
        # Read in the file if it exists
        print('\nLoading file into masterBias list')
        print(masterBiasFilename)
        masterBias = ai.MasterBias.read(masterBiasFilename)
        masterBias = masterBias.astype(np.float32)
        masterBiasDict.update({thisBin: masterBias})

        # If a masterBias was found, then berak out of the loop
        continue

    # Locate the raw bias images with this binning
    thisBiasBinBool = np.logical_and(biasBool, (binning == thisBin))
    thisBiasBinInds = np.where(thisBiasBinBool)
    biasImgFiles    = fileIndex['Filename'][thisBiasBinInds]

    # Otherwise proceed to read in and process a masterBias image
    print('\nProcessing {0} biases with ({1:g}x{1:g}) binning'.format(
        biasImgFiles.size, thisBin))

    # Loop through each of the files and add them to the biasImgList list
    biasImgList  = []
    for filename in biasImgFiles:
        # Read the raw bias image from the disk
        rawBias = ai.RawBias.read(filename)

        # Append the raw bias (overscan corrected) image to the list of biases
        biasImgList.append(rawBias)

    # Construct an ImageStack out of the bias image list
    biasStack = ai.ImageStack(biasImgList)

    # Compute the master bias
    masterBias = biasStack.combine_images()

    # Store the master bias in the dictionary of colibration data
    masterBiasDict.update({thisBin: masterBias})

    # Write masterBias object to disk
    masterBias.write(masterBiasFilename, dtype=np.float32, clobber=True)

    print('\nThe mean bias level is {0:g} counts\n'.
      format(masterBias.data.mean()))

    # Do a quick cleanup to make sure that memory survives
    del biasStack
    del masterBias

#==============================================================================
# ***************************** DARKS *****************************************
# Setup the paths to the dark images and compute the dark current map
#==============================================================================

masterDarkDict = {}
for thisBin in uniqBins:
    # Construct the filename for this dark.
    masterDarkFilename = 'MasterDark{0:g}.fits'.format(thisBin)
    masterDarkFilename = os.path.join(calibrationDir, masterDarkFilename)

    # Test if there is a MasterDark image for this binnng level.
    if os.path.isfile(masterDarkFilename):
        # Read in the file if it exists
        print('\nLoading file into masterDark list')
        print(masterDarkFilename)
        masterDark = ai.MasterDark.read(masterDarkFilename)
        masterDark = masterDark.astype(np.float32)
        masterDarkDict.update({thisBin: masterDark})

        # If a masterDark was found, then continue out of the loop
        continue

    # Locate the raw dark images with this binning
    thisDarkBinBool = np.logical_and(darkBool, (binning == thisBin))
    thisDarkBinInds = np.where(thisDarkBinBool)
    darkImgFiles    = fileIndex['Filename'][thisDarkBinInds]

    # Otherwise continue to read in and process a masterDark image
    print('\nProcessing {0} darks with ({1:g}x{1:g}) binning'.format(
        darkImgFiles.size, thisBin))

    # Loop through each of the files and add them to the darkImgList list
    darkImgList  = []
    for filename in darkImgFiles:
        # Read the raw dark image from disk
        rawDark = ai.RawDark.read(filename)

        # Apply the basic data processing
        reducedDark = rawDark.process_image(
            bias=masterBiasDict[thisBin]
        )

        # Divide the reduced dark by its own exposure time
        # (just to make sure that all darks have the same properties)
        reducedDark = reducedDark.divide_by_expTime()

        # Append the reduced dark to the list of dark images
        darkImgList.append(reducedDark)

    # Construct an ImageStack out of the dark image list
    darkStack = ai.ImageStack(darkImgList)

    # Compute the master dark
    masterDark = darkStack.combine_images()

    # Convert to a 32 bit float
    masterDark = masterDark.astype(np.float32)

    # Store the master dark in the dictionary of colibration data
    masterDarkDict.update({thisBin: masterDark})

    # Write masterDark object to disk
    masterDark.write(masterDarkFilename, dtype=np.float32, clobber=True)

    print('\nThe mean dark level is {0:g} counts\n'.
      format(masterDark.data.mean()))

    # Do a quick cleanup to make sure that memory survives
    del darkStack
    del masterDark

#==============================================================================
# ***************************** FLATS *****************************************
# Setup the paths to the flat images and compute the flat map
#==============================================================================
# Find the number of unique wavebands used in flats
uniqBands = np.unique(waveBand)

# Create an empty dictionary to store the masterFlatDict,
# keyed to each band/polAng/binning combination
masterFlatDict = {}

#Loop through each waveband
for thisBand in uniqBands:
    # Compute the unique values for the polaroid rotation angle
    thisFlatWaveBool = np.logical_and(flatBool, (waveBand == thisBand))
    thisFlatWaveInds = np.where(thisFlatWaveBool)
    uniqPolAngs      = np.unique(polAng[thisFlatWaveInds])

    for thisAng in uniqPolAngs:
        # Compute the unique values for the binning level
        thisFlatAngBool = np.logical_and(thisFlatWaveBool, (polAng == thisAng))
        thisFlatAngInds = np.where(thisFlatAngBool)
        uniqBins        = np.unique(binning[thisFlatAngInds]).astype(int)

        for thisBin in uniqBins:
            # Construct the flatKey and filename for this image
            flatKey            = (thisBand, thisAng, thisBin)
            flatKeyStr         = '{0:s}_{1:g}_{2:g}'.format(*flatKey)
            masterFlatFilename = 'MasterFlat' + flatKeyStr + '.fits'
            masterFlatFilename = os.path.join(calibrationDir, masterFlatFilename)

            # Test if the file exists
            if os.path.isfile(masterFlatFilename):
                # Read in the file if it exists
                print('\nLoading file into masterFlat list')
                print(masterFlatFilename)
                masterFlat = ai.MasterFlat.read(masterFlatFilename)
                masterFlat = masterFlat.astype(np.float32)
                masterFlatDict.update({flatKey: masterFlat})

                # If a master flat was found, then continue out of the loop...
                continue

            # Locate the raw bias images with this binning
            thisFlatBinBool = np.logical_and(thisFlatAngBool, (binning == thisBin))
            thisFlatBinInds = np.where(thisFlatBinBool)
            flatImgFiles    = fileIndex['Filename'][thisFlatBinInds]

            # Otherwise continue to read in and process a masterFlat image
            # Create the file if it does not exist
            print('\nProcessing {0} flats for'.format(flatImgFiles.size))
            print('band    = {0:s}'.format(thisBand))
            print('polAng  = {0:g}'.format(thisAng))
            print('binning = ({0:g}x{0:g})'.format(thisBin))

            # Loop through each of the files and add them to the biasImgList list
            flatImgList  = []
            for filename in flatImgFiles:
                # Read the raw flat image from disk
                rawFlat = ai.RawFlat.read(filename)

                # Apply basic data processing
                reducedFlat = rawFlat.process_image(
                    bias=masterBiasDict[thisBin],
                    dark=masterDarkDict[thisBin]
                )
                # Divide the flat image by its own mode
                reducedFlat = reducedFlat/reducedFlat.mode

                # Add the mode-normalized flat to the list of flat images
                flatImgList.append(reducedFlat)

            # Loop through all the flats and determine which ones have
            # lights-on vs. lights-off
            # Start by grabbing all the stats for each flat.
            flatStats = [img.sigma_clipped_stats() for img in flatImgList]
            flatStats = np.array(flatStats)

            # Now compute the SNR level in each flat
            flatSNRs = flatStats[:,1]/flatStats[:,2]

            # Estimate that the images with SNR > 1.5 are lights-on images
            thisFlatOnBool  = (flatSNRs > 1.5)
            thisFlatOnInds  = np.where(thisFlatOnBool)
            thisFlatOffBool = np.logical_not(thisFlatOnBool)
            thisFlatOffInds = np.where(thisFlatOffBool)

            # Convert the list of images into an array for indexing
            flatImgList = np.array(flatImgList)

            # Grab the lights-on flats and compute the average image
            if np.sum(thisFlatOnBool.astype(int)) > 0:
                # Grab the lights-on images
                flatOnImgList = flatImgList[thisFlatOnInds]

                # Construct an ImageStack out of the dark image list
                flatStack = ai.ImageStack(flatOnImgList)

                # Compute the master flat
                masterOnFlat = flatStack.combine_images()

                # Divide the master flat by its own mode (re-normalize)
                masterOnFlat = masterOnFlat/masterOnFlat.mode

                # Convert to a 32 bit float
                masterOnFlat = masterOnFlat.astype(np.float32)
            else:
                raise RuntimeError('There are no flat images. This is a major problem!')

            # Grab the lights-off flats and compute the average image
            if np.sum(thisFlatOffBool.astype(int)) > 0:
                # Grab the lights-on images
                flatOffImgList = flatImgList[thisFlatOnInds]

                # Construct an ImageStack out of the dark image list
                flatStack = ai.ImageStack(flatOffImgList)

                # Compute the master flat
                masterOffFlat = flatStack.combine_images()

                # Divide the master flat by its own mode (re-normalize)
                masterOffFlat = masterOffFlat/masterOffFlat.mode

                # Convert to a 32 bit float
                masterOffFlat = masterOffFlat.astype(np.float32)
            else:
                # Not having any lights-off flats probabbly isn't such a big deal
                flatOffImgList = []
                masterOffFlat = 0

            # Compute the difference between the lights-on and lights-off flats
            masterFlat = masterOnFlat - masterOffFlat

            # Make sure the header is fully up to date
            masterFlat._properties_to_header()

            # Catch any instances of where the flat is zero and set to one.
            zeroPix = (masterFlat.data == 0)
            if np.sum(zeroPix) > 0:
                # Copy the data array
                fixedArray = masterFlat.data.copy()

                # Find the pixels where there are zeros
                zeroInds   = np.where(zeroPix)
                fixedArray[zeroInds] = 1

                # Replace the master flat data array
                masterFlat.data = fixedArray

            # Store the master dark in the dictionary of colibration data
            masterFlatDict.update({flatKey: masterFlat})

            # Write masterFlat object to disk
            masterFlat.write(masterFlatFilename, dtype=np.float32, clobber=True)

            # Do a quick cleanup to make sure that memory survives
            del flatOnImgList
            del flatOffImgList
            del flatStack

# Let the user know everything completed
print('\n..........')
print('Finished producing master calibration fields')
