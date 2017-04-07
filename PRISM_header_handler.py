# This file provides a single place to store the appropriate header_handler for
# the PRISM instrument

def PRISM_header_handler(header):
    """Makes some small modifications to the header as its read in"""
    # Copy the header for manipulation
    outHeader = header.copy()

    # Change the location of the binning values
    outHeader['ADELX_01'] = header['CRDELT1']
    outHeader['ADELY_01'] = header['CRDELT2']

    # Make the ADU units lower case so that astropy.units can handle them.
    outHeader['BUNIT'] = header['BUNIT'].strip().lower()

    # Set the gain for the instrument (since it wasn't included in the header)
    # PRISM got a new CCD with a gain of 2.63 e-/ADU and a 7.30 e- read noise
    # as of January 28, 2015.
    # Retrieved from
    #
    # https://jumar.lowell.edu/confluence/display/Perkins/New+PRISM+CCD
    #
    # on April 7, 2017
    outHeader['AGAIN_01'] = 2.68
    outHeader['ARDNS_01'] = 7.30

    # # Before Janury 28, 2015, use this gain and read noise.
    # # Retrieved from old web-page... cannot find again...
    # outHeader['AGAIN_01'] = 3.3
    # outHeader['ARDNS_01'] = 13.0

    return outHeader
