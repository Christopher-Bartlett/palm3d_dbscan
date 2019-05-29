"""

getpalm3ddata.py

Extract x-y-z localisations coordinates from palm3d output file(s)
    Requires palm_3d.py from York et al. Nat. Meth. 2011
Save as NumPy array AllCorrectedLocalizations.npy

Written by Alistair Curd and Christopher Bartlett
University of Leeds 11-Sept-2017

"""

import cPickle
import shelve
from Tkinter import Tk
from tkFileDialog import askopenfilename
from tkFileDialog import askdirectory
import numpy as np

def getcorrecteddata(pkldata, locsfile, offset,
                     qualmin,
                     x_to_nm, y_to_nm, z_to_nm):
    """
    Get co-ordinates for selected data.

    z extremes truncated to avoid these gluts of poor data.
    Co-ordinates corrected for drift and converted to nm.

    Args:
        palm3d .pkl file
        palm3d locsfile
        x-y-z offset
        minimum correlation
        x-y-z nm conversion factors

    Returns:
       This data as a numpy array with shape (N, 3),
       where N is the number of localisations.
    """

    # Separate out offset components for use when corrected x, y and z, below.

    (xoffset, yoffset, zoffset) = offset

    # Get raw localisations (not drift-corrected).
    # locsimages is a shelve, or uber-list of images, called '0' '1', '2', etc.
    # Each image - locsimages['0'] or locsimages[repr(0)] -
    # is a list of localisations within that image
    # Each localisation - e.g. locsimages[repr(1)][1] -
    # is a dictionary of properties
    # Each localisation contains position, e.g. 'x', and other information

    locsimages = shelve.open(locsfile, protocol=2)

    # Initialise co-ordinates and data quality/correlation arrays
    # Over all of the images, find drift-corrected coordinates and 'qual'

    xs = []
    ys = []
    zs = []
    qual = []

    for imagecount in range(len(pkldata.images)):
        for loc in locsimages[repr(imagecount)]:
            # Exclude poor correlations with calibration stack
            if loc['qual'] > qualmin:
                # Exclude the extremes in z
                # (where incorrect localisations in z end up)
                if loc['z'] > 5 and loc['z'] < 75:
                    xs.append(loc['x'] - pkldata.drift['x'][imagecount] - xoffset)
                    ys.append(loc['y'] - pkldata.drift['y'][imagecount] - yoffset)
                    zs.append(loc['z'] - pkldata.drift['z'][imagecount] - zoffset)
                    qual.append(loc['qual']) # May use it later

    locsimages.close()

    # Numpy array of coordinates...with each (x,y,z) along one row

    xyz = np.array([xs, ys, zs])
    xyz = np.transpose(xyz)

    # Convert pixels to nm

    conv_to_nm = np.array([x_to_nm, y_to_nm, z_to_nm])
    xyz = np.multiply(xyz, conv_to_nm)

    # Save a little memory/time
    # We could usually get away with 16-bit
    # (350 pixels * 167 nm = 58450 nm < 65536)
    # But do want to keep 1 or decimal places as well.

    xyz = xyz.astype(np.float32)

    print '\n%i localisations found before subsetting.' % len(xyz)

    return xyz

def combineacquisitions(savefolder):
    """
    Combine data from multiple acquisitions.
    User input required to find acquistion data.

    Args:
        Folder in which to save data

    Returns:
       Numpy array as in getcorrecteddata().
    """

    print '\nPlease find me the first pickle (pkl)!'
    pklfile = askopenfilename()
    print '%s\n' % pklfile

    print 'Please find me the first localizations file.'
    print '(Usually in raw data folder, e.g. z=...)'
    print 'Use _palm_localizations for unlinked localisations'
    print 'or _palm_particles for localisations after linking.'
    locsfile = askopenfilename()
    print '%s' % locsfile

    # Get acquisition information, including drift.

    with open(pklfile, 'rb') as f:
        pkldata = cPickle.load(f)

    # Initial fiducial position, from which subsequent offsets are calculated

    veryfirstfidpos = pkldata.drift['initial_xyz']

    # Get pixel to distance conversion factors
    # and minimum correlation with calibration stack

    x_to_nm = float(raw_input('\nHow many nm per pixel in x? '))
    y_to_nm = float(raw_input('How many nm per pixel in y? '))
    z_to_nm = float(raw_input('How many nm per pixel in z? '))
    qualmin = float(raw_input('\nChoose minimum correlation for the data: '))

    print '\nReading these localisations. Please wait...'

    xyz = getcorrecteddata(pkldata=pkldata, locsfile=locsfile,
                           offset=np.array([0., 0., 0.]),
                           qualmin=qualmin,
                           x_to_nm=x_to_nm, y_to_nm=y_to_nm, z_to_nm=z_to_nm)

    while(raw_input(
            '\nDo you want to include another acqn. run? [y]/n: ') != 'n'):

        # Find the data, as before
        print '\nPlease find me that pickle then (pkl).'
        pklfile = askopenfilename()
        print '%s\n' % pklfile

        print 'Please find me relevant localizations file too.'
        print 'Use _palm_localizations for unlinked localisations'
        print 'or _palm_particles for localisations after linking.'
        locsfile = askopenfilename()
        print '%s' % locsfile

        with open(pklfile, 'rb') as f:
            pkldata = cPickle.load(f)

        print '\nReading these localisations. Please wait...'

        # Get the offset relative to the first run

        firstfidpos = np.array(pkldata.drift['initial_xyz'])
        offset = firstfidpos - veryfirstfidpos

        xyz_next = getcorrecteddata(pkldata, locsfile, offset,
                                    qualmin,
                                    x_to_nm, y_to_nm, z_to_nm)
        xyz = np.append(xyz, xyz_next, axis=0)

    np.save('%s\\AllCorrectedLocalizations' % savefolder, xyz)

    return None

def main():
    """
    Call getcorrecteddata() for each acqn. run,
    using combineacquisitions()
    """

    Tk().withdraw()

    print 'Where do you want to save the results?'
    savefolder = askdirectory()
    print '%s\n' % savefolder

    combineacquisitions(savefolder)

    return None

if __name__ == '__main__':
    main()
    print '\nHit Enter to exit'
    raw_input()
