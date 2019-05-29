"""

make_palm3d_histogram.py

Generate 3D histograms for visualisation
    Requires AllCorrectedLocalisations.npy
    Output of getpalm3ddata.py

Written by Christopher Bartlett and Alistair Curd
University of Leeds 11-Sept-2017

"""

from Tkinter import Tk
from tkFileDialog import askopenfilename
from tkFileDialog import askdirectory
import numpy as np

def make_histogram(savefolder, subset,
                   xmin_nm, xmax_nm,
                   ymin_nm, ymax_nm,
                   zmin_nm, zmax_nm,
                   title):
    """
    Save ROI as 3D histogram

    User input for bin size

    Args:
        Folder to save data
        ROI to save as histogram
        ROI x-y-z min/max positions
        Output file name
    """

    binsize = float(input('What is the histogram bin size? '))

    # Bins for histogram
    # Bins are centred on integer numbers of binsize (factor of 1.5)
    xmin_hist = int(xmin_nm / binsize) * binsize - 1.5 * binsize
    xmax_hist = int(xmax_nm / binsize) * binsize + 1.5 * binsize
    ymin_hist = int(ymin_nm / binsize) * binsize - 1.5 * binsize
    ymax_hist = int(ymax_nm / binsize) * binsize + 1.5 * binsize
    zmin_hist = int(zmin_nm / binsize) * binsize - 1.5 * binsize
    zmax_hist = int(zmax_nm / binsize) * binsize + 1.5 * binsize

    edgex = np.arange(xmin_hist, xmax_hist, step=binsize)
    edgey = np.arange(ymin_hist, ymax_hist, step=binsize)
    edgez = np.arange(zmin_hist, zmax_hist, step=binsize)
    hist, _ = np.histogramdd(subset, bins=(edgex, edgey, edgez))

    fname = '%s\\%s_HistogramAbsolutePositions-%i-%i-%i.raw' % (
        savefolder, title, hist.shape[2], hist.shape[1], hist.shape[0])
    hist.tofile(fname)
    print '\nHistogram saved (.raw)'

    return None

def crop(savefolder):
    """
    Crop region of interest

    User input required to select AllCorrectedLocalizations.npy file
    User input for min/max coordinate values from 350x350 image

    Args:
        Folder to save data
    """

    print '\nWhere is the AllCorrectedLocalisations.npy file?'
    xyzfile = askopenfilename()
    xyz = np.load(xyzfile)

    x_to_nm = float(input('\nHow many nm per pixel in x? '))
    y_to_nm = float(input('How many nm per pixel in y? '))
    z_to_nm = float(input('How many nm per pixel in z? '))
    conv_to_nm = np.array([x_to_nm, y_to_nm, z_to_nm])

    print '\nTo subset an area for analysis: '
    print 'use pixel positions defined in ImageJ'
    print 'from uncorrected reconstruction of the full FOV images'
    print 'e.g. nm/px in x = 100 and nm/px in y = 100'

    while True:
        # Coordinates for ROI
        # ----- REMEMBER ----- #
        # palm3d x-y coordinates are opposite to ImageJ
        # e.g x in palm3d = y in ImageJ
        ymin = float(input('\nMinimum x pixel (0 for full FOV): '))
        ymax = float(input('Maximum x pixel (350 for full FOV): '))
        xmin = float(input('\nMinimum y pixel (0 for full FOV): '))
        xmax = float(input('Maximum y pixel (350 for full FOV): '))
        zmin = float(input('\nMinimum z pixel (0 for full FOV): '))
        zmax = float(input('Maximum z pixel (80 for full FOV): '))

        title = ('Region_X' + str(int(ymin)) + '-' + str(int(ymax)) +
                 '_Y' + str(int(xmin)) + '-' + str(int(xmax)) +
                 '_Z' + str(int(zmin)) + '-' + str(int(zmax)))

        # Convert x-y-z to nm
        xmin_nm, ymin_nm, zmin_nm = np.multiply([xmin, ymin, zmin], conv_to_nm)
        xmax_nm, ymax_nm, zmax_nm = np.multiply([xmax, ymax, zmax], conv_to_nm)

        # Subset FOV into ROI
        xyzmin = np.array([xmin_nm, ymin_nm, zmin_nm])
        xyzmax = np.array([xmax_nm, ymax_nm, zmax_nm])
        subsetfilter = np.logical_and(xyz >= xyzmin, xyz <= xyzmax)
        chosenlocs = np.all(subsetfilter, axis=1)
        subset = xyz[chosenlocs]
        np.save('%s\\%s.npy' %(savefolder, title), subset)
        print '\nROI saved (.npy)'

        ans = raw_input('Do you want to save the ROI as a histogram? [y]/n: ')
        if ans != 'n':
            # Save histogram of ROI
            make_histogram(savefolder, subset,
                           xmin_nm, xmax_nm,
                           ymin_nm, ymax_nm,
                           zmin_nm, zmax_nm,
                           title)

        # Repeat cropping?
        ans = raw_input('\nDo you want to crop another region? [y]/n: ')
        if ans == 'n':
            break

    return None

def main():
    """
    Ask user for directory to save files

    Call crop() function
    """

    Tk().withdraw()

    print '\nWhere do you want to save the results?'
    savefolder = askdirectory()

    crop(savefolder)

    return None

if __name__ == '__main__':
    main()
    print '\nFinished!'
