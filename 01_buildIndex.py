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


#Setup the path delimeter for this operating system
delim = os.path.sep

#==============================================================================
# *********************** CUSTOM USER CODE ************************************
# this is where the user specifies where the raw data is stored
# and some of the subdirectory structure to find the actual .FITS images
#==============================================================================
# This is the location of the raw data for the observing run
rawPath = 'C:\\Users\\Jordan\\FITS Data\\PRISM_Data\\Raw_data\\'
# This is a list of strings containing the subdirectory structure for each night
subDirs = ['20150117\\PRISM_Images\\', \
           '20150118\\PRISM_Images\\', \
           '20150119\\PRISM_Images\\']

# This is a list of strings containing the prefix of the filenames
#filePrefix = ['20150117', '20150118', '20150119']

# This line prepends the rawPath variable to each element in the subDirs list
rawDirs = [rawPath + subDir for subDir in subDirs]

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
indexFile = reducedDir + delim + 'fileIndex.dat'
if not os.path.isfile(indexFile):
    #Loop through each night and build a list of all the files in observing run
    fileList = []
    for night in subDirs:
        nightPath = rawPath + night
        for file in os.listdir(nightPath):
            fileList.extend([os.path.join(nightPath, file)])

    #Sort the fileList
    fileNums = [''.join((file.split(delim).pop().split('.'))[0:2]) for file in fileList]
    sortInds = np.argsort(np.array(fileNums, dtype = np.int64))
    fileList = [fileList[ind] for ind in sortInds]

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
    fileIndex = Table([fileList, name, waveBand, polAng, binType],
                      names = ['Filename', 'Name', 'Waveband', 'Polaroid Angle', 'Binning'])
    fileIndex.add_column(Column(name='Use',
                                data=np.ones((numFiles)),
                                dtype=np.int),
                                index = 0)

    # Group by "Name"
    groupFileIndex = fileIndex.group_by('Name')

    # Grab the file-number orderd indices for the groupFileIndex
    fileIndices = np.argsort(groupFileIndex['Filename'])

    # Loop through each "Name" and assign it a "Target" value
    targetList = []
    ditherList = []
    for group in groupFileIndex.groups:
        # Select this groups properties
        thisName = np.unique(group['Name'])

        # Test if the group name truely is unique
        if len(thisName) == 1:
            thisName = thisName[0]
        else:
            print('There is more than one name in this group!')
            pdb.set_trace()

        # Count the number of elements in this group
        groupLen = len(group)

        # Add the "Target" column to the fileIndex
        thisTarget = input('\nEnter the target for group "{0}": '.format(thisName))
        thisTarget = [thisTarget]*groupLen

        # Ask the user to supply the dither pattern for this group
        thisDitherEntered = False
        while not thisDitherEntered:
            # Have the user select option 1 or 2
            print('\nEnter the dither patttern for group "{0}": '.format(thisName))
            thisDither = input('[1: ABBA, 2: HEX]')

            # Test if the numbers 1 or 2 were entered
            try:
                thisDither = np.int(thisDither)
                if (thisDither == 1) or (thisDither == 2):
                    # If so, then reassign as a string
                    thisDither = ['ABBA', 'HEX'][(thisDither-1)]
                    thisDitherEntered = True
            except:
                print('Response not recognized')

        # Create a list of "thisDither" entries
        thisDither = [thisDither]*groupLen

        # Add these elements to the target list
        targetList.extend(thisTarget)
        ditherList.extend(thisDither)

    pdb.set_trace()
    # Add the "Target" and "Dither columns"
    groupFileIndex.add_column(Column(name='Target',
                                data=np.array(targetList)),
                                index = 2)
    groupFileIndex.add_column(Column(name='Dither',
                                data=np.array(ditherList)),
                                index = 7)

    # Re-sort by file-number
    fileSortInds = np.argsort(groupFileIndex['Filename'])
    fileIndex1   = groupFileIndex[fileSortInds]

    # Write file to disk
    fileIndex1.write(indexFile, format='csv')

    # Write the file index to disk
    fileIndex = Table([fileList, obsType, waveBand, target, polAng, binType, lights],
                      names = ['Filename', 'Data',  'Waveband', 'Group',
                      'Polaroid Angle', 'Binning', 'Lights'])
    # ascii.write(fileIndex, indexFile)
    fileIndex.write(indexFile, format='csv')
