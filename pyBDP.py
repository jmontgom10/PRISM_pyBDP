import pdb
import psutil
import time
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from astropy.io import fits

class Image(object):
    """An object which stores an image array and header and provides a
    read method and an overscan correction method.
    """
    
    def __init__(self, filename=''):
        if len(filename) > 0:
            try:
                HDUlist     = fits.open(filename)
                self.header = HDUlist[0].header.copy()
                floatFlag   = self.header['BITPIX'] < 0
                numBits     = np.abs(self.header['BITPIX'])
                
                # Determine the appropriate data type for the array.
                if floatFlag:
                    if numBits >= 64:
                        dataType = np.float64
                    else:
                        dataType = np.float32
                else:
                    if numBits >= 64:
                        dataType = np.int64
                    elif numBits >= 32:
                        dataType = np.int32
                    else:
                        dataType = np.int16
                
                # Store the data as the correct data type
                self.arr      = HDUlist[0].data.astype(dataType, copy = False)

                # Check that binning makes sense and store it if it does
                uniqBins = np.unique((self.header['CRDELT1'],
                                     self.header['CRDELT2']))
                if len(uniqBins) > 1:
                    raise ValueError('Binning is different for each axis')
                else:
                    self.binning  = int(uniqBins[0])

                self.filename = filename
                self.dtype    = dataType
                HDUlist.close()

            except:
                print('File not found')
    
#        def read_file(self, filename):
    
    def copy(self):
        output = Image()
        output.arr = self.arr.copy()
        output.binning = self.binning
        output.header = self.header.copy()
        output.filename = self.filename.copy()
        output.dtype = self.dtype
        return output
    
    # TODO define a subtraction method for image array
    # test for if it's another image or a similarly sized array

    
    def overscan_correction(self, overscanPos, sciencePos,
                            overscanPolyOrder = 3):
        """Fits a polynomial to the overscan column and subtracts and extended
        polynomial from the entire image array.
        Note: The prescan region is not an accurate representation of the bias,
        so it is completely ignored.
        """
        # Check the binning
#        binning = np.unique((self.header['CRDELT1'], self.header['CRDELT2']))

        
        # Compute the rebinned locations of the pre-scan and post-scan regions.
        overscanPos1 = overscanPos/self.binning
        sciencePos1  = sciencePos/self.binning
#        overscanRegion = self.arr[overscanPos1[0][1]:overscanPos1[1][1], \
#                                  overscanPos1[0][0]:overscanPos1[1][0]]
        
        # Grab the pre-scan and post-scan regions of the array
        overscanRegion  = self.arr[:, overscanPos1[0][0]:overscanPos1[1][0]]
        
        # Make sure I know what is "PRESCAN" (right?) and "POSTSCAN" (left?)
        # and which one corresponds to the "overscan region"
        
        # Grab the shape of the array for future use
        ny, nx = self.arr.shape
        
        # Fit the correct order polynomial to the overscan columns
        overscanRowValues  = np.mean(overscanRegion, axis = 1)
        rowInds            = range(ny)
        overscanPolyCoeffs = np.polyfit(rowInds, overscanRowValues, overscanPolyOrder)
        overscanPolyValues = np.polyval(overscanPolyCoeffs, rowInds)
        
        # Expand the overscan along the horizontal axis and subtract.
        overscan  = np.tile(overscanPolyValues, (self.arr.shape[1],1)).T
        self.arr -= overscan
        
        # Trim the array to include only the science data
        self.arr  = self.arr[:, sciencePos1[0][0]:sciencePos1[1][0]]
    
    def stacked_average(imgList, clipSigma = 3.0):
        """Compute the median filtered mean of a stack of images.
        Standard deviation is computed from the variance of the stack of
        pixels.
        
        parameters:
        imgList   -- a list containing Image class objects.
        clipSigma -- the level at which to trim outliers (default = 3)
        """
        numImg = len(imgList)
        print('\nEntered averaging method')
        if numImg > 1:
            # Test for the correct number of bits in each pixel
            dataType    = imgList[0].dtype
            if dataType == np.int16:
                numBits = 16
            elif (dataType == np.int32) or (dataType == np.float32):
                numBits = 32
            elif (dataType == np.int64) or (dataType == np.float64):
                numBits = 64

            # Compute the number of pixels that fit under the memory limit.
            memLimit    = (psutil.virtual_memory().available/
                          (numBits*(1024**2)))
            memLimit    = int(10*np.floor(memLimit/10.0))
            numStackPix = memLimit*(1024**2)*8/numBits
            ny, nx      = imgList[0].arr.shape
            numRows     = int(np.floor(numStackPix/(numImg*nx)))
            if numRows > ny: numRows = ny
            numSections = int(np.ceil(ny/numRows))
            
            # Compute the number of subsections and display stats to user
            print('\nAiming to fit each stack into {0:g}MB of memory'.format(memLimit))
            print('\nBreaking stack of {0:g} images into {1:g} sections of {2:g} rows'
              .format(numImg, numSections, numRows))
            
            # Initalize an array to store the final averaged image
            finalImg = np.zeros((ny,nx))
            
            # Compute the stacked average of each section
            #
            #
            # TODO Check that this section averaging is working correctly!!!
            #
            #
            for thisSec in range(numSections):
                # Calculate the row numbers for this section
                thisRows = (thisSec*numRows,
                            min([(thisSec + 1)*numRows, ny]))
                
                # Stack the selected region of the images.
                secRows = thisRows[1] - thisRows[0]
                stack   = np.ma.zeros((numImg, secRows, nx), dtype = dataType)
                flipCount = np.zeros((numImg, secRows, nx), dtype = np.int8)
                for i in range(numImg):
                    stack[i,:,:] = imgList[i].arr[thisRows[0]:thisRows[1],:]

                print('\nAveraging rows {0[0]:g} through {0[1]:g}'.format(thisRows))
            
                # Iteratively clip outliers until answer converges.
                # Use the stacked median for first image estimate.
                imgEstimate = np.median(stack, axis = 0)
                stackSigma  = stack.std(axis = 0)
                outliers    = np.ndarray(stack.shape, dtype = bool)
                outliers1   = outliers.copy()
               
                # Setup vairables to track if the mask has converged
                # and how many times the iteration 
                numMaskChange = 1
                iterCounter   = 0
                
                # This loop will iterate until the mask converges to an
                # unchanging state, or until 40 iterations  have completed.
                while (numMaskChange != 0) and (iterCounter < 12):
                    print('\n\tBeginning iteration {0:g}'.format(iterCounter+1))
                    
                    # Loop through the stack, and find the outliers.
                    for j in range(numImg):
                        deviation        = np.absolute(stack[j,:,:] - imgEstimate)
                        outliers1[j,:,:] = deviation > clipSigma*stackSigma
                    
                    # Count the number of changed mask elements,
                    # and set the new mask to match these outliers.
                    maskChange    = outliers != outliers1
                    flipCount    += maskChange
                    numMaskChange = np.sum(maskChange)
                    stack.mask    = outliers1
                    print('\t{0:g} mask elements have changed'.format(numMaskChange))
                    
                    #Compute new stack statistics
                    stackSigma   = stack.std(axis = 0)
                    imgEstimate1 = stack.mean(axis = 0)
                    
                    #Save old variables for use in the next iteration
                    outliers     = outliers1.copy()
                    imgEstimate  = imgEstimate1.copy()
                    iterCounter += 1
                
                # Now that this section has been averaged, store it in output.
                finalImg[thisRows[0]:thisRows[1],:] = imgEstimate
                
                # Check for pixels where the clipping eliminated entire stack.
                blockedPix = (np.sum(outliers, axis = 0) == numImg)
                if np.sum(blockedPix) > 0:
                    print('\n\tReplacing blocked pixels with median estimator')
                    medianImg            = np.median(stack, axis = 0)
                    finalImg[blockedPix] = medianImg[blockedPix]
                    del medianImg
            return finalImg
        else:
            return imgList[0].arr

            
class Bias(Image):
    """A subclass of the "Image" class: stores bias images and provides some
    methods for bias type operations.
    """
    
    def __init__(self, filename = ''):
        super(Bias, self).__init__(filename)
    
    def average(self):
        return np.mean(self.arr)
    
    def master_bias(biasList, clipSigma = 3.0):
        return Image.stacked_average(biasList, clipSigma = clipSigma)

    def overscan_polynomial(biasList, overscanPos):
        ## Loop through biasList and build a stacked average of the
        ## overscanRegion
        numBias      = len(biasList)
#        overscanList = biasList.copy()
        binning      = np.unique([(biasList[i].header['CRDELT1'], 
                                   biasList[i].header['CRDELT2'])
                                   for i in range(numBias)])
        if len(binning) > 1:
            raise ValueError('Binning is different for each axis')
        overscanPos1   = overscanPos/binning

        ## Trim the overscan list to contain ONLY the overscan.
        overscanList = []
        for i in range(len(biasList)):
            overscanList.append(biasList[i].copy())
            overscanList[i].arr = \
              biasList[i].arr[overscanPos1[0][1]:overscanPos1[1][1], \
                                  overscanPos1[0][0]:overscanPos1[1][0]]
        masterOverscan = Image.stacked_average(overscanList)

        ## Average across each row.
        overscanRowValues = np.mean(masterOverscan, axis = 1)
        rowInds           = range(len(overscanRowValues))

        ## Compute statistics for a zeroth order polynomial.
        polyCoeffs = np.polyfit(rowInds, overscanRowValues, 0)
        polyValues = np.polyval(polyCoeffs, rowInds)
        chi2_m     = np.sum((overscanRowValues - polyValues)**2)

        ## Loop through polynomial degrees and store Ftest result
#        alpha  = 0.05   #Use a 5% "random probability" as a cutoff
        sigma  = 7.0    #Use a 7 sigma requirement for so much data
        alpha  = (1 - scipy.stats.norm.cdf(sigma))
        Ftests = []
        coeffs = []
        for deg in range(1,10):
            dof        = len(overscanRowValues) - deg - 1
            polyCoeffs = np.polyfit(rowInds, overscanRowValues, deg)
            coeffs.append(polyCoeffs)
            polyValues = np.polyval(polyCoeffs, rowInds)
            chi2_m1    = np.sum((overscanRowValues - polyValues)**2)
#            print('reduced Chi2 = ' + str(chi2_m1/dof))
            Fchi       = (chi2_m - chi2_m1)/(chi2_m1/dof)
            prob       = 1 - scipy.stats.f.cdf(Fchi, 1, dof)
            Ftests.append(prob < alpha)
            
            ## Store chi2_m1 in chi2 for use in next iteration
            #***** NOTE HERE *****
            # there are so many degrees of freedom (v2 ~ 20000 pixels!)
            # that Fchi(v1=1, v2) = Fchi(v1=1, v2=Inf), so almost all
            # polynomials will pass the test (until we get to Fchi < 3.84
            # for alpha = 0.05). We can decrease alpha (stricter test), or
            # we can rebin along the column
            chi2_m     = chi2_m1
        
        
        ## Find the LOWEST ORDER FAILED F-test and return the degree
        ## of the best fitting polynomial
        bestDegree = np.min(np.where([not test for test in Ftests]))
        
#        # Plot all of the fits...
#        for plotNum in range(len(coeffs)):
#            cf = coeffs[plotNum]
#            plt.subplot(3, 3, plotNum + 1)
#            plt.plot(rowInds, overscanRowValues)
#            polyValues = np.polyval(cf, rowInds)
#            plt.plot(rowInds, polyValues, linewidth = 4, color = 'r')
#            plt.text(1250, 1029, str(Ftests[plotNum]))
#        plt.show()
        return bestDegree
            
#    @property
#    def read_noise(self):
#        return np.std(self.arr)

class Flat(Image):
    """A subclass of the "Image" class: stores flat frames and provides some
    methods for flat type operations.
    """
    
    def __init__(self, filename = ''):
        super(Flat, self).__init__(filename)
    
    def master_flat(flatList, clipSigma = 3.0):
        return Image.stacked_average(flatList, clipSigma = clipSigma)

class Dark(Image):
    """A subclass of the "Image" class: stores dark frames and provides some
    methods for dark type operations.
    """
    
    def __init(self, filename = ''):
        super(Dark, self).__init__(filename)
    
    def dark_time(darkList):
        
        # Build a list of all the exposure times
        expTimeList = [dark.header['EXPTIME'] for dark in darkList]

        # Check if all the exposure times are the same...
        expTimes = np.unique(expTimeList)
        if expTimes.size > 1:
            raise ValueError('More than one exposure time found in list')
        else:
            return expTimes[0]
    
    def dark_current(darkList, clipSigma = 3.0):
        avgImg   = Image.stacked_average(darkList, clipSigma = clipSigma)
        darkTime = Dark.dark_time(darkList)
        return avgImg/darkTime
