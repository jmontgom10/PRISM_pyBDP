#Apply the calibration fields to the RAW science data to produce
#reduced science data.

#Import whatever modules will be used
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
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

# Reduced directory (for saving the final images)
reducedDir = os.path.join(pyBDP_data, 'pyBDP_reduced_images')
if (not os.path.isdir(reducedDir)):
    os.mkdir(reducedDir, 0o755)

# Read the fileIndex back in as an astropy Table
print('\nReading file index from disk')
indexFile = os.path.join(pyBDP_data, 'rawFileIndex.csv')
fileIndex = Table.read(indexFile, format='csv')

biasBool = (fileIndex['Data'] == 'BIAS')
darkBool = (fileIndex['Data'] == 'DARK')
flatBool = (fileIndex['Data'] == 'FLAT')
sciBool  = (fileIndex['Data'] == 'OBJECT')
waveBand = fileIndex['Waveband']
polAng   = fileIndex['Polaroid Angle']
binning  = fileIndex['Binning']
lights   = fileIndex['Lights']

#==============================================================================
# ***************************** BIAS *****************************************
# Setup the paths to the bias images and compute the bias map
#==============================================================================

# Find the number of unique binnings used in biases
uniqBins = np.unique(binning)

masterBiasDict = {}
for thisBin in uniqBins:
    # Construct the filename for this bias.
    masterBiasFilename = 'MasterBias{0:g}.fits'.format(thisBin)
    masterBiasFilename = os.path.join(calibrationDir, masterBiasFilename)

    # Read in the masterBias file
    print('\nLoading file into masterBias list')
    print(masterBiasFilename)

    masterBias = ai.MasterBias.read(masterBiasFilename)
    masterBias = masterBias.astype(np.float32)
    masterBiasDict.update({thisBin: masterBias})


#==============================================================================
# ***************************** DARKS *****************************************
# Setup the paths to the dark images and compute the dark current map
#==============================================================================

masterDarkDict = {}
for thisBin in uniqBins:
    # Construct the filename for this dark.
    masterDarkFilename = 'MasterDark{0:g}.fits'.format(thisBin)
    masterDarkFilename = os.path.join(calibrationDir, masterDarkFilename)

    # Read in the file
    print('\nLoading file into masterDark list')
    print(masterDarkFilename)

    masterDark = ai.MasterDark.read(masterDarkFilename)
    masterDark = masterDark.astype(np.float32)
    masterDarkDict.update({thisBin: masterDark})

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

            # Read in the masterFlat file
            print('\nLoading file into masterFlat list')
            print(masterFlatFilename)
            masterFlat = ai.MasterFlat.read(masterFlatFilename)
            masterFlat = masterFlat.astype(np.float32)
            masterFlatDict.update({flatKey: masterFlat})

#==============================================================================
# **************************** SCIENCE ****************************************
# Setup the paths to the
#==============================================================================
# Grab the indices of the science images
scienceInds     = np.where(sciBool)
scienceImgFiles = fileIndex['Filename'][scienceInds]

print('\nBeginning to reduce science data.')
for filename in scienceImgFiles:
    # Read in the raw science image from disk
    rawScience = ai.RawScience.read(filename)

    # Extract the information on this file
    thisBand  = rawScience.filter
    thisAng   = rawScience.header['POLPOS']
    thisBin   = np.int(np.unique(rawScience.binning))

    # Construct the key for the flat dictionary
    thisFlatKey = (thisBand, thisAng, thisBin)

    # Process this raw science image using the relevant calibration files
    reducedScience = rawScience.process_image(
        bias=masterBiasDict[thisBin],
        dark=masterDarkDict[thisBin],
        flat=masterFlatDict[thisFlatKey]
    )

    # Write the file to disk
    reducedFileName = os.path.join(
        reducedDir,
        os.path.basename(rawScience.filename)
    )
    reducedScience.write(reducedFileName, dtype=np.float32, clobber=True)

# Let the user know everything is finished
print('\n..........')
print('Finished reducing science data.')
