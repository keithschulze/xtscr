#
#  Spinal Cord Recovery XTension
#
#  Copyright (c) 2014 Keith Schulze (keith.schulze@monash.edu)
#  MIT-style copyright and disclaimer apply
#
#    <CustomTools>
#      <Menu>
#       <Item name="SCR" icon="Python" tooltip="Spinal Cord Recovery">
#         <Command>PythonXT::xtscr(%i)</Command>
#       </Item>
#      </Menu>
#    </CustomTools>
#

""" Imaris Xtension for the analysis of Spinal Cord Recovery in live
confocal images.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d
import seaborn as sns
from coord import coord
from numpy import linalg
from pIceImarisConnector import pIceImarisConnector as ice
from scipy.spatial import distance


__author__ = "Keith Schulze"
__copyright__ = "Copyright 2015, Monash Micro Imaging"
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "keith.schulze@monash.edu"
__status__ = "development"


def get_lesion_coords(mp):
    """Get lesion coords from an Imaris measurement points object.

    Parameters
    ----------
    mp: MeasurementPoints object
        Measurement points object from Imaris

    Return
    ------
    out: ndarray
         Numpy array in format [[x, y, z]]
    """
    return np.array(mp.GetPositionsXYZ(), np.float)


def get_spots_odr(spot_df):
    """Determine central axis of the spinal coord by fitting a ODR.

    Parameters
    ----------
    spot_df: pandas.DataFrame
             DataFrame with the columns for x, y and z coordinate for each
             spot. Also needs a column, time, specifying the time points of
             each coordinate.

    Return
    ------
    line: numpy.ndarray
          A 3 x 2 numpy array containing 2 3D coords specifying the start and
          end coords of a ODR line through the spots coordinates.
    """
    time0 = spot_df.groupby('time')[["x", "y", "z"]].get_group(0.0).values
    coord_mean = time0.mean(axis=0)

    # Use SVD to determine first principal component through the mean of the
    # spot data.
    _, _, vv = np.linalg.svd(time0 - coord_mean)

    # Estimate the max distance between the two furtherest coords
    xmin = np.min(spot_df.x)
    xmax = np.max(spot_df.x)
    ymin = np.min(spot_df.y)
    ymax = np.max(spot_df.y)
    zmin = np.min(spot_df.z)
    zmax = np.max(spot_df.z)
    max_dist = distance.cdist(np.array([[xmin, ymin, zmin]]),
                              np.array([[xmax, ymax, zmax]]))
    hmd = max_dist[0, 0]/2

    # Get begin and end coords of the line
    linepts = vv[0] * np.mgrid[-hmd:hmd:2j][:, np.newaxis]
    linepts += coord_mean

    return linepts


def project_lc(lesion, line):
    """Project the estimated lesion coord onto the closest point of a line.

    Parameters
    ----------
    lesion: numpy.ndarray
            A 3 x 1 (x, y, z) numpy array specifying the estimated lesion
            coordinate.
    line: numpy.ndarray
          A 3 x 2 numpy array containing 2 3D coords specifying the start and
          end coords of a ODR line through the spots coordinates.

    Return
    ------
    p_lesion: numpy.ndarray
              A 3 x 1 (x, y, z) numpy array specifying the nearest
              coordinate along the line to the given lesion coordinate
              i.e., the projected coordinate.
    """
    line_v = np.subtract(line[1], line[0])
    lesion_v = np.subtract(lesion[0], line[0])
    proj_lc = (np.dot(lesion_v, line_v)/np.dot(line_v, line_v)) * line_v

    return proj_lc + line[0]


def xtscr(appid):
    """Entry function for the Imaris XTension."""
    conn = ice(appid)

    # Get the lesion coords from the measurement points object
    mp = conn.getAllSurpassChildren(False, "MeasurementPoints")[0]
    if not mp:
        print "No measurement points marking lesion site specified"
        return
    ref_points = get_lesion_coords(mp)

    # Get spot coords from the selected spots object
    spots = conn.getSurpassSelection("Spots")
    if not spots:
        print "No spots were found. Please select spots object you want to"\
              " analyse."
        return
    spot_coords = np.array(spots.GetPositionsXYZ(), dtype=np.float)

    # Create a pandas dataframe to summarise data
    spot_df = pd.DataFrame(spot_coords.view(dtype=[('x', np.float),
                                                   ('y', np.float),
                                                   ('z', np.float)]).ravel())

    # Add tracks to the dataframe
    spot_df["track"] = np.zeros(len(spot_df["x"]))
    spot_track_edges = spots.GetTrackEdges()
    spot_track_ids = spots.GetTrackIds()

    for i, e in zip(spot_track_ids, spot_track_edges):
        if spot_df.track[e[0]] == 0:
            spot_df.track[e[0]] = i
        if spot_df.track[e[1]] == 0:
            spot_df.track[e[1]] = i

    spot_df.track = spot_df.track.astype(np.int).astype(np.str)

    spot_df['time'] = np.array(spots.GetIndicesT(), np.float) *\
        conn.mImarisApplication.GetDataSet().GetTimePointsDelta()

    # Get coordinates of line along central axis.
    line_coords = get_spots_odr(spot_df)

    # Project estimated lesion coords onto the central axis line.
    lesion_coords = project_lc(ref_points, line_coords)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(*spot_coords.T)
    ax.scatter3D(*ref_points.T, s=300, c=u'r')
    ax.scatter3D(*lesion_coords, s=300, c=u'g')
    ax.scatter3D(*line_coords.T, s=50, c=u'g')
    ax.plot3D(*line_coords.T)
    ax.set_xlim(0, 260)
    ax.set_ylim(0, 260)
    ax.set_zlim(-0, -200)
    fig.set_facecolor('white')
    fig.savefig("/Users/skeith/Desktop/cart3d.png", bbox_inches='tight',
                facecolor=fig.get_facecolor(), transparent=True)

    # Adjust origin of spot, line and lesion coords to be the start line coord.
    spot_coords = np.subtract(spot_coords, line_coords[0])
    lesion_coords = np.subtract(np.array([lesion_coords]), line_coords[0])
    ref_points = np.subtract(ref_points, line_coords[0])
    line_coords = np.subtract(line_coords, line_coords[0])

    # Switch to spherical coordinate system to estimate the angle of elevation
    # of the central axis
    spot_coords = coord.cartesian2spherical(spot_coords)
    lesion_coords = coord.cartesian2spherical(lesion_coords)
    ref_points = coord.cartesian2spherical(ref_points)
    line_coords = coord.cartesian2spherical(line_coords)

    # Adjust all coords by substracting angle of elevation of the best fit line,
    # thus making central axis run along Z-axis
    spot_coords[:, 2] = spot_coords[:, 2] - line_coords[1, 2]
    lesion_coords[:, 2] = lesion_coords[:, 2] - line_coords[1, 2]
    ref_points[:, 2] = ref_points[:, 2] - line_coords[1, 2]
    line_coords[:, 2] = line_coords[:, 2] - line_coords[1, 2]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(*coord.spherical2cartesian(spot_coords).T)
    ax.scatter3D(*coord.spherical2cartesian(lesion_coords).T, s=50, c=u'r')
    ax.scatter3D(*coord.spherical2cartesian(ref_points).T, s=50, c=u'r')
    ax.scatter3D(*coord.spherical2cartesian(line_coords).T, s=50, c=u'g')
    ax.plot3D(*coord.spherical2cartesian(line_coords).T)
    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    ax.set_zlim(0, 300)
    fig.set_facecolor('white')
    fig.savefig("/Users/skeith/Desktop/sph3d.png", bbox_inches='tight',
                facecolor=fig.get_facecolor(), transparent=True)


    # Switch to cylindrical coordinate system.
    spot_coords = coord.spherical2cylindrical(spot_coords)
    lesion_coords = coord.spherical2cylindrical(lesion_coords)
    ref_points = coord.spherical2cylindrical(ref_points)
    line_coords = coord.spherical2cylindrical(line_coords)

    # Reset origin to lesion coord.
    spot_coords[:, 2] = spot_coords[:, 2] - lesion_coords[:, 2]
    line_coords[:, 2] = line_coords[:, 2] - lesion_coords[:, 2]
    ref_points[:, 2] = ref_points[:, 2] - lesion_coords[:, 2]
    lesion_coords[:, 2] = lesion_coords[:, 2] - lesion_coords[:, 2]

    # Reset azimuth angle to that of second reference point
    spot_coords[:, 1] = spot_coords[:, 1] - ref_points[1, 1]
    line_coords[:, 1] = line_coords[:, 1] - ref_points[1, 1]
    lesion_coords[:, 1] = lesion_coords[:, 1] - ref_points[1, 1]
    ref_points[:, 1] = ref_points[:, 1] - ref_points[1, 1]

    # Output 3D plot of cylindrical registered scr data
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter3D(*coord.cylindrical2cartesian(spot_coords).T)
    ax.scatter3D(*coord.cylindrical2cartesian(lesion_coords).T, s=50, c=u'r')
    ax.scatter3D(*coord.cylindrical2cartesian(ref_points).T, s=50, c=u'r')
    ax.scatter3D(*coord.cylindrical2cartesian(line_coords).T, s=50, c=u'g')
    ax.plot3D(*coord.cylindrical2cartesian(line_coords).T)
    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    ax.set_zlim(-150, 150)
    fig.set_facecolor('white')
    fig.savefig("/Users/skeith/Desktop/cyl3d.png", bbox_inches='tight',
                facecolor=fig.get_facecolor(), transparent=True)
