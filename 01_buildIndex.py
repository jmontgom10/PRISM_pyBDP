# This script will generate an index of the raw files

#Import whatever modules will be used
import os
import sys
import time
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

################################################################################
# Define a recursive file search which takes a parent directory and returns all
# the FILES (not DIRECTORIES) beneath that node.
def recursive_file_search(parentDir, exten='', fileList=[]):
    # Query the elements in the directory
    subNodes = os.listdir(parentDir)

    # Loop through the nodes...
    for node in subNodes:
        # If this node is a directory,
        thisPath = os.path.join(parentDir, node)
        if os.path.isdir(thisPath):
            # then drop down recurse the function
            recursive_file_search(thisPath, exten, fileList)
        else:
            # otherwise test the extension,
            # and append the node to the fileList
            if len(exten) > 0:
                # If an extension was defined,
                # then test if this file is the right extension
                exten1 = (exten[::-1]).upper()
                if (thisPath[::-1][0:len(exten1)]).upper() == exten1:
                    fileList.append(thisPath)
            else:
                fileList.append(thisPath)

    # Return the final list to the user
    return fileList
################################################################################

#Setup the path delimeter for this operating system
delim = os.path.sep

#==============================================================================
# *********************** CUSTOM USER CODE ************************************
# this is where the user specifies where the raw data is stored
# and some of the subdirectory structure to find the actual .FITS images
#==============================================================================
# This is the location of the raw data for the observing run
rawDir   = 'C:\\Users\\Jordan\\FITS_data\\PRISM_data\\raw_data\\201612\\'
fileList = recursive_file_search(rawDir, exten='.fits')

#Sort the fileList
fileNums = [''.join((file.split(delim).pop().split('.'))[0:2]) for file in fileList]
fileNums = [file.split('_')[0] for file in fileNums]
sortInds = np.argsort(np.array(fileNums, dtype = np.int64))
fileList = [fileList[ind] for ind in sortInds]

# Define the path to the parent directory for all pyBDP products
pyBDP_data = 'C:\\Users\\Jordan\\FITS_data\\PRISM_data\\pyBDP_data\\201612\\'

# Setup new directory for reduced data if it does not already exist
if (not os.path.isdir(pyBDP_data)):
    os.mkdir(pyBDP_data, 0o755)

# Define the directory into which the average calibration images will be placed
calibrationDir = os.path.join(pyBDP_data, 'master_calibration_images')

# Reduced directory (for saving the final images)
reducedDir = os.path.join(pyBDP_data, 'pyBDP_reduced_images')

#==============================================================================
# ***************************** INDEX *****************************************
# Build an index of the file type and binning, and write it to disk
#==============================================================================

# Define the path to the rawFileIndex.csv file
indexFile = os.path.join(pyBDP_data, 'rawFileIndex.csv')
# Test for image type
print('\nCategorizing files into "BIAS", "DARK", "FLAT", "OBJECT"\n')
startTime = time.time()
# Begin by initalizing some arrays to store the image classifications
name     = []
data     = []
binType  = []
polAng   = []
waveBand = []
lights   = []
fileCounter = 0
percentage  = 0

#Loop through each file in the fileList variable
for file in fileList:

    # Classify each file type and binning
    data.append(fits.getval(file, 'OBSTYPE'))
    tmpName = fits.getval(file, 'OBJECT')
    if len(tmpName) < 1:
        tmpName = 'blank'
    name.append(tmpName)
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

endTime = time.time()
numFiles = len(fileList)
print('\nFile processing completed in {0:g} seconds'.format(endTime -startTime))

# Query the user about the targets of each group...
# Write the file index to disk
fileIndex = Table([fileList, name, data, waveBand, polAng, binType, lights],
                  names = ['Filename', 'Name', 'Data', 'Waveband', 'Polaroid Angle', 'Binning', 'Lights'])

# Write file to disk
fileIndex.write(indexFile, format='csv')
