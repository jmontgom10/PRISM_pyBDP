"""
Generates an index of the raw files.
"""

#Import whatever modules will be used
import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from astropy.io import fits

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
#==============================================================================
# *********************** CUSTOM USER CODE ************************************
# this is where the user specifies where the raw data is stored
# and some of the subdirectory structure to find the actual .FITS images
#==============================================================================
# This is the location of the raw data for the observing run
rawDir   = 'C:\\Users\\Jordan\\FITS_data\\PRISM_data\\raw_data\\201612\\'
fileList = np.array(recursive_file_search(rawDir, exten='.fits'))

#Sort the fileList
fileNums = [''.join((os.path.basename(f).split('.'))[0:2]) for f in fileList]
fileNums = np.array([f.split('_')[0] for f in fileNums], dtype=np.int64)
sortInds = fileNums.argsort()
fileList = fileList[sortInds]

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
OBJECT   = []
OBSTYPE  = []
TELRA    = []
TELDEC   = []
FILTER   = []
POLPOS   = []
BINNING  = []
NIGHT    = []
fileCounter = 0
percentage  = 0

#Loop through each file in the fileList variable
for filename in fileList:
    # Grab the OBJECT header value
    tmpOBJECT = fits.getval(filename, 'OBJECT')
    if len(tmpOBJECT) < 1:
        tmpOBJECT = 'blank'
    OBJECT.append(tmpOBJECT)

    # Grab the TELRA header value
    try:
        TELRA.append(fits.getval(filename, 'TELRA'))
    except:
        TELRA.append(0)

    try:
        TELDEC.append(fits.getval(filename, 'TELDEC'))
    except:
        TELDEC.append(0)

    # Grab the OBSTYPE header value
    OBSTYPE.append(fits.getval(filename, 'OBSTYPE'))

    # Grab the FILTNME3 header value
    FILTER.append(fits.getval(filename, 'FILTNME3'))

    # Grab the POLPOS header value
    POLPOS.append(fits.getval(filename, 'POLPOS'))

    # Grab the binning from the CRDELT values
    tmpBINNING = fits.getval(filename, 'CRDELT*')
    if tmpBINNING[0] == tmpBINNING[1]:
        BINNING.append(int(tmpBINNING[0]))
    else:
        raise ValueError('Binning is illogical')

    # Grab the EXPTIME value from the

    # Assign a NIGHT value for this image
    NIGHT.append(''.join((os.path.basename(filename).split('.'))[0]))

    # Count the files completed and print update progress message
    fileCounter += 1
    percentage1  = np.floor(fileCounter/len(fileList)*100)
    if percentage1 != percentage:
        print('completed {0:3g}%'.format(percentage1), end="\r")
    percentage = percentage1

endTime = time.time()
numFiles = len(fileList)
print('\nFile processing completed in {0:g} seconds'.format(endTime - startTime))

# Query the user about the targets of each group...
# Write the file index to disk
fileIndex = Table([fileList,   NIGHT,    OBJECT,   TELRA,   TELDEC,   OBSTYPE,   FILTER,   POLPOS,   BINNING],
          names = ['FILENAME', 'NIGHT', 'OBJECT', 'TELRA', 'TELDEC', 'OBSTYPE', 'FILTER', 'POLPOS', 'BINNING'])

# Write file to disk
fileIndex.write(indexFile, format='ascii.csv', overwrite=True)
