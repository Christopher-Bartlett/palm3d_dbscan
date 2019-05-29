"""

run_dbscan.py

DBSCAN cluster analysis on ROIs
Remove fiducial locs by filtering out locs
    > 500 neighbours within 30 nm

INPUT:
        Directory containing region(s) of interest .npy file(s)

MEASUREMENTS:
        num:      number of locs per cluster
        xymedian: median pairwise distance between locs,
                  within a cluster in x-y

OUTPUT:
        Cluster measurements saved as single file:
            AllMeasurements.csv

        Average cluster from nomalised x-y-z coordinates as .dat file:
            AverageCluster-%i-%i-%i.dat
            %i-%i-%i are histogram shape dimensions in x-y-z

        Number of clusters per ROI as single file:
            AllROI_NumberOfClusters.csv

Written by Christopher Bartlett and Alistair Curd
University of Leeds 11-Sept-2017

"""

import time
import os
from Tkinter import Tk
from tkFileDialog import askdirectory
import numpy as np
from sklearn.cluster import DBSCAN
from scipy.spatial import distance

def removefiducial(data):
    """
    Remove fiducial localisations from palm3d files

        Search filter applied to localisations
        If localisation has > 1000 locs within [40,40,60] nm cuboid
            - considered part of fiducial
        Save all localisations other than fiducial for DBSCAN analysis

    Args:
        x-y-z palm3d localisation data as NumPy array

    Returns:
        Array of same shape as input without fiducial localisations

    """

    # Distance to search for localisations
    filterarray = np.array([40, 40, 60])

    x_coords = []
    y_coords = []
    z_coords = []
    loc = 0

    while loc < len(data):

        loc_a = data[loc]

        # A filter to search for localisation coordinates within filter distance
        # Returns an np array of True / False
        testfilter = np.logical_and(data > loc_a - filterarray,
                                    data < loc_a + filterarray)

        chosenlocs = np.all(testfilter, axis=1)
        subxyz = data[chosenlocs]

        if len(subxyz) < 1000:
            x_coords = np.append(x_coords, loc_a[0])
            y_coords = np.append(y_coords, loc_a[1])
            z_coords = np.append(z_coords, loc_a[2])

        loc += 1

    xyz = np.array([x_coords, y_coords, z_coords])
    xyz = np.transpose(xyz)
    print 'Locs classified as fiducial and removed: %i' %(len(data) - len(xyz))

    return xyz

def run_dbscan(data, neighbourhood, nlocs):

    """
    DBSCAN clustering on ROI

    Args:
        x-y-z palm3d localisation data as NumPy array
        Search radius for DBSCAN
        Min points for DBSCAN

    Returns:
        Arrays containing measurements for ROI
    """

    # DBSCAN clustering
    clusterer = DBSCAN(eps=neighbourhood, min_samples=nlocs).fit(data)

    # Combine XYZ coordinates with cluster identity from DBSCAN
    # Coordinates identified as "noise" are labelled with cluster id of -1
    data_dbscan = np.c_[data, clusterer.labels_]

    id_ = []
    num = []
    xymedian = []
    xavg_cluster = []
    yavg_cluster = []
    zavg_cluster = []

    for cluster_id in range(clusterer.labels_.max() + 1):

        xyzclust = data_dbscan[data_dbscan[:, 3] == cluster_id]

        loc_a = 0
        while loc_a < len(xyzclust):

            # Register localisation coordinates to cluster centroid position
            xnorm = np.subtract(xyzclust[loc_a, 0], np.mean(xyzclust[:, 0]))
            ynorm = np.subtract(xyzclust[loc_a, 1], np.mean(xyzclust[:, 1]))
            znorm = np.subtract(xyzclust[loc_a, 2], np.mean(xyzclust[:, 2]))

            # Only save coordinates if fall inside 500 nm cuboid around centroid
            if xnorm > -250 and xnorm < 250:
                if ynorm > -250 and ynorm < 250:
                    if znorm > -250 and znorm < 250:

                        xavg_cluster = np.append(xavg_cluster, xnorm)
                        yavg_cluster = np.append(yavg_cluster, ynorm)
                        zavg_cluster = np.append(zavg_cluster, znorm)

            loc_a += 1

        # Calculate pairwise Euclidean distances between localisations in XY
        xydist = distance.pdist(xyzclust[:, 0:2], 'euclidean')

        id_ = np.append(id_, cluster_id)
        num = np.append(num, len(xyzclust))
        xymedian = np.append(xymedian, np.median(xydist))

    total_clusters = clusterer.labels_.max() + 1
    print 'Clusters in this ROI: %i' % total_clusters

    return (id_, num, xymedian,
            xavg_cluster, yavg_cluster, zavg_cluster,
            total_clusters)

def main():
    """
    Calls removefiducial function
    Calls run_dbscan function
    
    User input to specify DBSCAN clustering parameters
    """
    Tk().withdraw()

    file_list = []
    filename_total = []
    id_total = []
    num_total = []
    xymedian_total = []
    xavg_cluster_total = []
    yavg_cluster_total = []
    zavg_cluster_total = []
    filename_roi = []
    num_clusters = []

    print '\nWhere are the .npy region(s) of interest file(s)?'
    directory = askdirectory()

    directory_list = os.listdir(directory)
    for item in directory_list:
        if item.startswith('Region') and item.endswith('.npy'):
            file_list.append(item)

    print '%i file(s) in this folder to to be analysed' % len(file_list)

    # Parameters for DBSCAN clustering
    search_radius = int(raw_input('\nWhat is the DBSCAN search radius?: '))
    min_cluster_size = int(raw_input('How many locs within this radius?: '))

    start_time = time.time()
    for item in file_list:
        print '\nAnalysing region %i/%i: %s' %(file_list.index(item) + 1,
                                               len(file_list),
                                               item[7:-4])
        xyz = np.load(directory + '/' + item)

        # Remove fiducial locs
        xyz = removefiducial(data=xyz)

        # Run DBSCAN analysis
        (id_, num, xymedian,
         xavg_cluster, yavg_cluster, zavg_cluster,
         total_clusters) = run_dbscan(data=xyz,
                                      neighbourhood=search_radius,
                                      nlocs=min_cluster_size)

        # Keep output from run_dbscan for each region of interest
        filename_total = np.append(filename_total, np.full(len(id_), item[7:-4]))
        id_total = np.append(id_total, id_)
        num_total = np.append(num_total, num)
        xymedian_total = np.append(xymedian_total, xymedian)

        xavg_cluster_total = np.append(xavg_cluster_total, xavg_cluster)
        yavg_cluster_total = np.append(yavg_cluster_total, yavg_cluster)
        zavg_cluster_total = np.append(zavg_cluster_total, zavg_cluster)

        filename_roi = np.append(filename_roi, item[7:-4])
        num_clusters = np.append(num_clusters, total_clusters)

    print '\nFinished DBSCAN analysis'
    print '%.2f sec to analyse all ROIs' % (time.time() - start_time)

    # Save all measurements as AllMeasurements.csv
    parameters = np.transpose(np.array([filename_total,
                                        id_total,
                                        num_total,
                                        xymedian_total]))

    np.savetxt('%s\\AllMeasurements.csv' % directory, parameters,
               fmt='%s', delimiter=',',
               header='ROI,cluster ID,number of locs,xymedian')

    # Save avg cluster as 3D histogram
    # Histogram - 4 nm px bins within 500 nm bounding box
    # %i-%i-%i = histogram shape for ImageJ .dat import
    avg_cluster_total = np.transpose(np.array([xavg_cluster_total,
                                               yavg_cluster_total,
                                               zavg_cluster_total]))

    edge = np.arange(-250, 254, step=4)
    hist, _ = np.histogramdd(avg_cluster_total, bins=(edge, edge, edge))

    hist.tofile('%s\\AverageCluster-%i-%i-%i.dat' % (directory,
                                                     hist.shape[2],
                                                     hist.shape[1],
                                                     hist.shape[0]))

    # Save number of clusters / region
    all_clusters_roi = np.transpose(np.array([filename_roi, num_clusters]))

    np.savetxt('%s\\AllROI_NumberOfClusters.csv' % directory, all_clusters_roi,
               fmt='%s', delimiter=',',
               header='ROI, Number of Clusters')

    return None

if __name__ == '__main__':
    main()
