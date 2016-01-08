# This script will generate an index of the raw files

#Import whatever modules will be used
import os
import sys
import pdb
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.table import Column
from astropy.io import fits, ascii
from scipy import stats

# Add the AstroImage class
sys.path.append("C:\\Users\\Jordan\\Libraries\\python\\AstroImage")
from AstroImage import AstroImage, Bias, Dark, Flat


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
rawDir   = 'C:\\Users\\Jordan\\FITS Data\\PRISM_Data\\Raw_data\\'
fileList = recursive_file_search(rawDir, exten='.fits')

#Sort the fileList
fileNums = [''.join((file.split(delim).pop().split('.'))[0:2]) for file in fileList]
fileNums = [file.split('_')[0] for file in fileNums]
sortInds = np.argsort(np.array(fileNums, dtype = np.int64))
fileList = [fileList[ind] for ind in sortInds]

# Define the directory into which the average calibration images will be placed
calibrationDir = 'C:\\Users\\Jordan\\FITS Data\\PRISM_Data\\Calibration'

# Reduced directory (for saving the final images)
# This is where the "fileIndex.csv" file will be located
reducedDir = 'C:\\Users\\Jordan\\FITS Data\\PRISM_Data\\Reduced_data'

# These are the overscan regions for all PRISM frames at 1x1 binning
#                       ((x1,y1), (x2, y2))
overscanPos = np.array([[2110, 8],[2177, 2059]], dtype = np.int32)
sciencePos  = np.array([[70,  32],[2070, 2032]], dtype = np.int32)

#==============================================================================
# ***************************** INDEX *****************************************
# Build an index of the file type and binning, and write it to disk
#==============================================================================
# Check if a file index already exists... if it does then just read it in
indexFile = rawDir + delim + 'fileIndex.csv'

# Test for image type
print('\nCategorizing files into "BIAS", "DARK", "FLAT", "OBJECT"\n')
startTime = os.times().elapsed
# Begin by initalizing some arrays to store the image classifications
obsType  = []
name     = []
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

endTime = os.times().elapsed
numFiles = len(fileList)
print('\nFile processing completed in {0:g} seconds'.format(endTime -startTime))

# Query the user about the targets of each group...
# Write the file index to disk
fileIndex = Table([fileList, name, waveBand, polAng, binType, lights],
                  names = ['Filename', 'Name', 'Waveband', 'Polaroid Angle', 'Binning', 'Lights'])

# Re-sort by file-number
fileSortInds = np.argsort(fileIndex['Filename'])
fileIndex1   = fileIndex[fileSortInds]

# Write file to disk
fileIndex1.write(indexFile, format='csv')
