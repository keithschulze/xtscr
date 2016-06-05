#
#  Spinal Cord Recovery XTension
#
#  Copyright (c) 2016 Keith Schulze (keith.schulze@monash.edu)
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
import math
import os
import time
import threading
import Queue
import numpy as np
import pandas as pd
from pIceImarisConnector import pIceImarisConnector as ice
from scr import scr

import Tkinter
import ttk
import tkFileDialog
import tkMessageBox


__author__ = "Keith Schulze"
__copyright__ = "Copyright 2015, Monash Micro Imaging"
__license__ = "MIT"
__version__ = "0.1.0"
__email__ = "keith.schulze@monash.edu"
__status__ = "development"



class SCRGUI(object):

    """docstring for SCRGUI"""
    def __init__(self, master, queue):
        self.queue = queue
        self.labeltext = Tkinter.StringVar()
        self.labeltext.set("SCR Processing")
        label = Tkinter.Label(master, width=50, textvariable=self.labeltext)
        self.prog = ttk.Progressbar(master, orient='horizontal', length=300,
                                    mode='determinate')
        label.pack(padx=10, pady=10)
        self.prog.pack(padx=10, pady=10)

    def processincoming(self):
        while self.queue.qsize():
            try:
                msg, prog = self.queue.get(0)
                self.labeltext.set(msg)
                self.prog.step(prog)
            except Queue.Empty:
                pass


class SCRProcessor(object):
    """docstring for SCRProcessor"""

    def __init__(self, master, conn):
        super(SCRProcessor, self).__init__()
        self.master = master
        self.conn = conn
        self.queue = Queue.Queue()
        self.gui = SCRGUI(master, self.queue)
        self.running = True

        self.dataset_name, _ = os.path.splitext(
            os.path.basename(conn.mImarisApplication.GetCurrentFileName())
        )

        dir_options = {'mustexist': True,
                       'title': 'Please select an output directory'}
        self.output_dir = tkFileDialog.askdirectory(**dir_options)
        if len(self.output_dir) == 0:
            self.master.quit()
            raise Exception("Open folder cancelled")

        threading.Thread(target=self.workerthread).start()
        self.periodiccall()

    def periodiccall(self):
        self.gui.processincoming()
        if self.running:
            self.master.after(10, self.periodiccall)

    def extract_sc_centres(self, t):
        """Extract central axes for a specified time point.

        Parameters
        ----------
        t: int
           timepoint

        Returns
        -------
        axes: np.ndarray (float)
              2 x 3 array storing coords for the start and end coords of
              the central axis for the given timepoint.
        """
        self.queue.put(("Extracting central axis for timepoint: " + str(t), 1))
        vol = self.conn.getDataVolume(0, t)
        return scr.extract_sc_centres(vol, 1, slices=50)


    def workerthread(self):
        if not os.path.exists(self.output_dir):
            print "Folder path does not exist"
            tkMessageBox.showwarning(
                "nonexistent folder",
                "Folder does not exist!"
            )
            time.sleep(2)
            self.master.quit()
            raise Exception("Folder not found.")

        self.queue.put(("Getting dataset dimensions", 1))
        # Get dataset parameters
        xsize, ysize, zsize, csize, tsize = self.conn.getSizes()
        xvox, yvox, zvox = self.conn.getVoxelSizes()
        xmin, _, ymin, _, zmin, _ = self.conn.getExtends()

        # Determine central axis in each frame
        central_axis_cache = os.path.join(self.output_dir,
                                          self.dataset_name+".npy")

        if os.path.exists(central_axis_cache):
            self.queue.put(("Extracting central axis for each time point.", 90))
            sc_centre_coords = np.load(central_axis_cache)
        else:
            sc_centre_coords = np.array(
                [self.extract_sc_centres(t) for t in range(tsize)])

            # Convert to physical coords using voxel sizes
            sc_centre_coords[:,:,0] = sc_centre_coords[:,:,0] * xvox + xmin
            sc_centre_coords[:,:,1] = sc_centre_coords[:,:,1] * yvox + ymin
            sc_centre_coords[:,:,2] = sc_centre_coords[:,:,2] * zvox + zmin

            np.save(central_axis_cache, sc_centre_coords)

        # Get the lesion coords from the measurement points object
        self.queue.put(("Retrieve lesion coordinates", 1))
        mp = self.conn.getAllSurpassChildren(False, "MeasurementPoints")[0]
        if not mp:
            print "No measurement points marking lesion site specified"
            tkMessageBox.showerror("No lesion site",
                                   "Please use measurement points to indicate"
                                   " lesion site and centre of the notochord"
                                   " perpendicular to lesion.")
            self.master.quit()
            return

        lesion_coords = SCRProcessor.get_lesion_coordinates(mp)

        # Get spot coords from the selected spots object
        self.queue.put(("Retrieve spot coordinates", 1))

        spots = None

        try:
            spots = self.conn.getSurpassSelection("Spots")
        except Exception, e:
            tkMessageBox.showerror("No Spots object selected",
                                   "No spots object found or selected.")
            self.master.quit()
            return

        if not spots:
            print "No spots were found. Please select spots object you want"\
                  "to analyse."
            tkMessageBox.showerror("No Spots",
                                   "No spots object founds or not selected."
                                   " Please select the spots object to want"
                                   " analyse")
            self.master.quit()
            return

        spot_coords = np.array(spots.GetPositionsXYZ(), dtype=np.float)

        # Create a pandas dataframe to summarise data
        self.queue.put(("Creating dataframe", 1))
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
            self.conn.mImarisApplication.GetDataSet().GetTimePointsDelta()

        # Add time index
        spot_df['tindex'] = np.array(spots.GetIndicesT())

        # Register to cylindrical coordinates
        self.queue.put(("Registering on cylindrical coordinate system", 1))
        spot_df[['r', 'l', 'theta']] = scr.register_cylinder(spot_df, sc_centre_coords,
                                                             lesion_coords)
        spot_df['abs_r'] = np.abs(spot_df['r'])
        spot_df['theta_deg'] = np.rad2deg(np.add(spot_df.theta, np.pi))

        self.queue.put(("Saving dataframe", 1))
        spot_df.to_csv(os.path.join(self.output_dir, self.dataset_name+".csv"))
        self.running = False
        self.master.quit()

    @classmethod
    def get_lesion_coordinates(cls, measurement_points):
        """ Get Lesion coordinates from a measurement points object.

        Parameters
        ----------
        measurement_points: MeasurementPoints object
                            MeasurementPoints object from Imaris

        Returns
        -------
        lesion_coords: numpy.ndarray
                       2 x 3 ndarray array representing coordinates for
                       lesion site and centre of the notochord perpendicular
                       to the lesion site.
        """
        return np.array(measurement_points.GetPositionsXYZ(), np.float)


def xtscr(appid):
    """Entry function for the Imaris XTension."""
    conn = ice(appid)

    if not conn.isAlive():
        print "Could not connect to Imaris"
        tkMessageBox.showwarning("Connection failed",
                                 "Could not connect to Imaris!")
        time.sleep(2)
        return

    root = Tkinter.Tk()
    client = SCRProcessor(root, conn)
    root.mainloop()
