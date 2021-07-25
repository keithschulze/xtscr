'''
SCR module for the analysis of Spinal Cord Recovery from live imaging
of the Zebrafish spinal cord.
'''
import math
import operator
import numpy as np
import pandas as pd
from scr import coord
from scipy import optimize
from skimage.filters import threshold_otsu, gaussian as gaussian_filter
from skimage.measure import label, regionprops
from skimage.morphology import closing, disk


def project_coord_onto_line(coordinate, line):
    """Project a coord onto the closest point of a line.

    Parameters
    ----------
    coord: numpy.ndarray
            A 3 x 1 (x, y, z) numpy array specifying a coordinate.
    line: numpy.ndarray
          A 3 x 2 numpy array containing 2 3D coords specifying the start and
          end coords of a ODR line through the spots coordinates.

    Return
    ------
    proj_coord: numpy.ndarray
                A 3 x 1 (x, y, z) numpy array specifying the nearest
                coordinate along the line to the given lesion coordinate
                i.e., the projected coordinate.
    """
    line_v = np.subtract(line[1], line[0])
    coord_v = np.subtract(coordinate, line[0])
    proj_coord = (np.dot(coord_v, line_v) / np.dot(line_v, line_v)) * line_v

    return proj_coord + line[0]


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters.
    Adapted from http://wiki.scipy.org/Cookbook/FittingData

    Parameters
    ----------
    height: float
            Peak height of Gaussian curve
    center_x: float
              X coordinate of center of the Gaussian curve
    center_y: float
              Y coordinate of the center of the Gaussian curve
    width_x: float
             Width in X of the Gaussian curve
    width_y: float
             Width in Y of the Gaussian curve

    Returns
    -------
    gauss: function
           Gaussian function that return the amplitude for a given x and y
           coordinate
    """
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x, y: height*np.exp(-(((center_x-x)/width_x)**2 +
                                      ((center_y-y)/width_y)**2)/2)


def moments(data):
    """Estimates the height, center_x, center_y, width_x, width_y
    gaussian parameters of a 2D distribution by calculating its
    moments.

    Parameters
    ----------
    data: ndarray
          2D array of image intensity data

    Returns
    -------
    height, height, center_x, center_y, width_x, width_y: float
    """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y


def fitgaussian(data, params=None):
    """Estimates height, center_x, center_y, width_x, width_y
    gaussian parameters of a 2D Gaussian distribution by
    least squares fit.

    Parameters
    ----------
    data: ndarray
          2D array of image intensity data
    params: (float, float, float, float, float)
            (Optional) Tuple of best guess starting params for
            height, center_x, center_y, width_x, width_y to best
            optimized. If not provided, these will be estimated from
            the image data.
    """
    if not np.all(params):
        params = moments(data)

    def errorfunction(p):
        return np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return p


def spinal_cord_boundingbox(img):
    """Determine the coordinates x_min, y_min, x_max and y_max coordinates of a
    bounding box around the intensity data of cross-sections of the end of the
    spinal cord. If more that one region is detected, then the largest region
    is selected.

    Parameters
    ----------
    img: np.array
         2D array of image intensity data

    Returns
    -------
    bbox: tuple
          (x_min, y_min, x_max, y_max)
    """
    thresh = threshold_otsu(img)
    bw = closing(img > thresh, disk(5))
    label_img = label(bw)

    regions = [region for region in regionprops(label_img, img)
               if region.area > 200 and region.max_intensity > thresh]

    max_idx = 0
    if len(regions) < 1:
        raise Exception("no region for spinal cord detected")
    elif len(regions) > 1:
        max_idx, _ = max([(i, r.area) for i, r in enumerate(regions)],
                         key=operator.itemgetter(1))

    return regions[max_idx].bbox


def extract_sc_centres(vol, axis=0, sigma=10, slices=10):
    """Extracts the end coordinates of central axis (i.e., a line) of a spinal
    cord running along the specified axis. A specifed number of slices/planes
    perpendicular to the specified axis are taken from each end of the volumme,
    and max projected into a single plane representing the shape of the
    intensity data at each end. These smoothed and then fitted with a 2D
    Gaussian curve. The center coordinates of the fitted Gaussian are returned.

    Parameters
    ----------
    vol: array
         Z x Y x X 3D array of spinal cord intensity data
         representing a volume.
    axis: int
          Axis along which the spinal cord runs.
    sigma: int
           Width of gaussian filter to use for smoothing.
    slices: int
            Number of slices to max project at each end of the spinal cord.
    """
    axis_size = vol.shape[axis]
    if slices > axis_size:
        raise Exception("Number of specified slices: %s exceeds the dimension"
                        "along axis: %s" % (slices, axis))

    vol_swap = np.swapaxes(vol, 0, axis)
    top = 0
    bottom = axis_size - 1

    # Max projection
    top_max = np.amax(vol_swap[top:top+slices,:,:], axis=0)
    bottom_max = np.amax(vol_swap[bottom-slices:bottom+1,:,:], axis=0)

    # Smooth using gaussian
    tef = gaussian_filter(top_max, sigma)
    bef = gaussian_filter(bottom_max, sigma)

    # Dtermine bounding box
    teminz, teminx, temaxz, temaxx = spinal_cord_boundingbox(tef)
    beminz, beminx, bemaxz, bemaxx = spinal_cord_boundingbox(bef)
    tef_crop = tef[teminz:temaxz, teminx:temaxx]
    bef_crop = bef[beminz:bemaxz, beminx:bemaxx]

    # Fit gaussian and extract center coords
    paramt = fitgaussian(tef_crop)
    paramt[1] = paramt[1] + teminz
    paramt[2] = paramt[2] + teminx
    top_coord = [paramt[2], paramt[1]]
    top_coord.insert(axis, 0)
    paramb = fitgaussian(bef_crop)
    paramb[1] = paramb[1] + beminz
    paramb[2] = paramb[2] + beminx
    bottom_coord = [paramb[2], paramb[1]]
    bottom_coord.insert(axis, axis_size)
    return [tuple(top_coord), tuple(bottom_coord)]


def vector_rotation_matrix(axis, angle):
    '''Calculate rotation matrix to rotate vector in 3D around a specified
    origin and by a specified angle via Euler-Rodrigues formula.

    Parameters
    ----------
    axis: ndarray, float
          3 x 1 Numpy array describing the axis (origin) of rotation
    angle: float
           rotation angle in radian

    Returns
    -------
    rot: ndarray, float
         3x3 array representing the rotation matrix
    '''
    u = axis/math.sqrt(np.dot(axis, axis))
    a = np.cos(angle/2)
    b, c, d = -u * np.sin(angle/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    return np.array([[aa + bb - cc - dd, 2*(b*c - a*d), 2*(b*d + a*c)],
                     [2*(b*c + a*d), aa + cc - bb - dd, 2*(c*d - a*b)],
                     [2*(b*d - a*c), 2*(c*d + a*b), aa + dd - bb - cc]])


def register_cylinder(spots, axes, lesion):
    """ Register cell coordinates relative to central axis and lesion coords at
    each time point, and convert to cylindrical coordinate system.

    Parameters
    ----------
    spots: pandas.DataFrame
           DataFrame with x, y, z, tindex, trackid (columns) for each spot.
    axes:  numpy.ndarray
           n x 2 x 3 array representing the top and bottom coords of the spinal
           at each timepoint (n).
    lesion: numpy.ndarray
            1 x 3 array representing the estimate lesion coord at each
            timepoint 1.

    Returns
    -------
    spots: pandas.DataFrame, float
           DataFrame containing registered cyclindrical
           coordinates (r, l, theta) for each spot.
    """
    # Adjust axes coords to be relative to the base axes coord
    ascc = np.zeros_like(axes, dtype=np.float)
    ascc[:,1] = np.subtract(axes[:,1], axes[:,0])
    ascc[:,0] = np.subtract(axes[:,0], axes[:,0])

    ## Determine origin and zenith angle of central axis at each time point.
    # Switch to spherical coord system.
    ascc_sph = np.zeros_like(ascc)
    ascc_sph[:,0] = coord.cartesian2spherical(ascc[:,0])
    ascc_sph[:,1] = coord.cartesian2spherical(ascc[:,1])

    # Axis origin at each time
    rot_axis = np.ones_like(ascc_sph[:,1])
    rot_axis[:,1] = ascc_sph[:,1,1] - np.pi/2
    rot_axis[:,2] = np.pi/2
    rot_axis = coord.spherical2cartesian(rot_axis)

    # Zenith angle of axis at each time
    zen_angle = ascc_sph[:,1,2]

    # Adjust lesion coord and project onto central axis at each timepoint
    alc = np.subtract(lesion, axes[0,0])
    palc = np.array([project_coord_onto_line(alc[0], a) for a in ascc])

    # Rotate transform lesion coords to align with Z axis
    palc = np.array([np.dot(lc, vector_rotation_matrix(rot_axis[i], zen_angle[i]))
                     for i, lc in enumerate(palc)])

    # Adjust spots to be relative to base of central axis at each timepoint
    spots[['x', 'y', 'z']] = spots.groupby("time").\
        apply(lambda df: np.subtract(df[['x', 'y', 'z']],
                                     axes[df.tindex.iloc[0],0]))

    # Rotate transform spot coords to align with Z axis
    def rotate(df):
        i = df.tindex.iloc[0]
        coords = df[['x', 'y', 'z']]
        out = np.dot(coords, vector_rotation_matrix(rot_axis[i], zen_angle[i]))
        return pd.DataFrame({'time': df.time,
                             'tindex': df.tindex,
                             'x': out[:,0],
                             'y': out[:,1],
                             'z': out[:,2]})

    spots[['time', 'tindex', 'x', 'y', 'z']] = spots.groupby("time").apply(rotate)

    # Adjust spots coords to be relative to lesion coord (origin) at each time
    # point.
    spots[['x', 'y', 'z']] = spots.groupby("time").\
        apply(lambda df: np.subtract(df[['x', 'y', 'z']],
                                     palc[df.tindex.iloc[0]]))

    # Calculate cylindrical coordinates
    spots['r'], spots['theta'], spots['l'] = \
        zip(*coord.cartesian2cylindrical(spots[['x', 'y', 'z']].values))

    return spots[['r', 'l', 'theta']]
