# pyBDP

Python PRISM reduction pipeline

# Dependencies and Setup

pyBDP requires installation of the "AstroImage" class. This can be downloaded
from the jmontgom10 GitHub repository and has its own set of dependencies, which
the user should take care to install properly (see the AstroImage README File a
guide on install those packages).

I recommend using the Anaconda environment, as that comes with numpy, scipy,
astropy, and matplotlib preinstalled. If you elect not to use Anaconda, then
make sure to get those packages properly installed before proceeding to install
the AstroImage dependencies.

# Procedure

This package contains three scripts for reducing CCD data from the PRISM
instrument on the 1.8 m Perkins telescope. Each script is numbered, and the
scripts should be executed in their numbered order, as in
```$ python 01_buildIndex.py
```

Each script contains some code which needs to be modified by the user before
execution. Most of this is simply pointing the computer to the correct
directories to find the raw data as well as supplying the directory where the
reduced data should be stored. To help keep the directory structure clean, I
recommend using a directory structure similar to the default values provided in
these scripts.

## 01_buildIndex.py

This script loops through each of the \*.fits files in the raw data directory
and builds an inventory of Bias, Dark, Flat, and Science images as well as any
binning information.

If the directories you provided are well formulated, then this script should
simply execute and take about 1 minute to run (modulo Moore's law).

## 02_buildCalibration.py

*Wait!* Before you execute this script, be sure to update the directory
structure in the first 40 lines of code. I know, I know... I should have put all
that information in a "configuration file". Maybe I will do that in the next
version.

This script, when properly configured, will compute **average** bias, dark, and
wavelength and polaroid angle specific flat-field images. This procedure takes
advantage of the "stacked_average" method in the AstroImage class, which in turn
uses the vectorized masked numpy arrays and averages along the "z-axis". This
takes a few minutes but should run smoothly with no user interaction.

## 03_reduceScienceData

Begin by once more configuring the directory structure at the top of the script.
This script will apply a basic de-biasing and flattening using the calibration
images computed by the previous script. The reduced science image files are
written to the corresponding directory where they can be addressed by the pyPol
scripts.

The current version of the script *does not subtract dark current*, as the
dark-current images obtained indicated that the dark current was less than the
read-noise. Thus, subtracting the dark-current would simply introduce extraneous
noise with no improvement in the accuracy of the measurements.
