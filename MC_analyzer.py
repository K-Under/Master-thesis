import numpy as np
import xml.etree.ElementTree as ET
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

"""
This script contains functions to analyze the .mifx and .raw files from the MCS.
"""


# needed for the conversion
MIFXTYPES2NPTYPES = {
    'uint8': np.uint8,
    'int8': np.int8,
    'int16': np.int16,
    'uint16': np.uint16,
    'int32': np.int32,
    'uint32': np.uint32,
    'real32': np.float32,
    'real64': np.float64
}


# loading the file
def load_mifx(pathmifx):
    """
    Loads the .raw file as a numpy array given a .mifx file. Naming needs to be consistent amongst files
    :param pathmifx: path to .mifx file
    :return: data stored in .raw
    """
    with open(pathmifx.with_suffix('.raw'), 'rb') as f_in:
        tree = ET.parse(pathmifx)
        root = tree.getroot()
        shape_x = int(root.find('NumberOfXVoxels').text)
        shape_y = int(root.find('NumberOfYVoxels').text)
        shape_z = int(root.find('NumberOfZVoxels').text)
        datatype = root.find('BinaryType').text
        data = np.frombuffer(f_in.read(), dtype=MIFXTYPES2NPTYPES[datatype])
        data = data.reshape((shape_z, shape_y, shape_x))
        return data


################################
#    big / useful function and class
################################

# We define a class to store data more easily, reuse functions and stored variables and for easier user access
class BreastAnalyzer:
    def __init__(self, data, voxel_size=0.3, standard_current=32):
        self.data = data  # pure data as a numpy array
        self.voxel_size = voxel_size  # in mm
        self.standard_current = standard_current  # in mA

        # dimensions
        self.x_length = self.data.shape[2]
        self.y_length = self.data.shape[1]
        self.z_length = self.data.shape[0]

        # breast data and geometry
        self.breast_data = self.ring_remover()  # breast_data (data without ring)
        self.diameter_lyst = self.breast_volume_diameter_estimate()  # breast diameter lyst
        self.breast_length = self.breast_length_estimate()  # length of the breast (mm)
        self.breast_midpoint = self.reference_point_find()  # midpoint index of breast (rounded down)

        #self.breast_mean_dose = self.dose_mean()  # filtered somewhat




    def ring_remover(self):  # removes the ring from the data set (nr.9) and turns them to air
        data_copy = self.data.copy()
        data_copy[data_copy == 9] = 0  # removes the ring
        return data_copy


    # breast diameter and length estimation, slice or breast based

    def breast_slyce_diameter_estimate(self, slyce):
        # estimates the breast width given a xy slice of a segmentation file
        # needs a slice with no ring in it! (from self.breast_data)

        if np.all(slyce == 0):
            return 0

        boundary_points = np.argwhere(slyce)
        hull = ConvexHull(boundary_points)  # gets the surface of the shape as a scipy ConvexHull object (not indices themselves)
        hull_points = boundary_points[hull.vertices]  # selects only the point needed to create the surface (shape (points, 2))


        # calculating the diameter
        diff = hull_points[:, None, :] - hull_points[None, :, :]  # use numpy broadcasting to take the difference of all possibilities
        dist_sq = np.sum(diff ** 2, axis=2)  # s quare each possibility
        diameter = np.sqrt(dist_sq.max())  # take the max value and take the root
        diameter = diameter*self.voxel_size  # convert to mm

        return diameter


    def breast_volume_diameter_estimate(self):  # applies the breast_slyce_diameter_estimate function to the whole z-axis of the breast
        width_lyst = []
        for z in range(self.z_length):
            width_lyst.append(self.breast_slyce_diameter_estimate(self.breast_data[z]))
        return width_lyst


    def breast_length_estimate(self):  # estimates the breast length based on the breast volume
        length = 0
        for z in range(self.z_length):  # go slice by slice
            slyce = self.breast_data[z]
            if not np.all(slyce == 0):  # check if the slice contains breast tissue (includes all tissue)
                length += self.voxel_size
        return length


    # breast dose estimation

    @staticmethod
    def bi_exp(a, b, k1, k2, c, x):  # bi-exponential from Shin et al., x can be anything here
        return (a*np.exp(k1*x))+(b*np.exp(k2*(x**2)))+c


    def dose_slyce_estimate(self, diameter):  # estimates the MGD for a diameter (other metrics possible)
        return self.bi_exp(9.22, 0, -0.076, 0, 3.1, diameter)  # constants from Shin et all.


    def dose_volume_estimate(self):
        dose_lyst = []
        for z in range(self.z_length):  # go slice by slice
            diameter = self.diameter_lyst[z]
            if diameter != 0:
                dose_lyst.append(self.dose_slyce_estimate(diameter))
            else:
                dose_lyst.append(0)
        return dose_lyst


    # reference point finding and dose reduction value / tube current modulation

    def reference_point_find(self):  # finds the mid-breast reference point for a slyce
        # checks where slices are nonzero and returns the indices of where this is the case
        occupied_slices = np.where(
            np.any(self.breast_data != 0, axis=(1, 2))
        )[0]

        if len(occupied_slices) == 0:
            return None

        first_slice = occupied_slices[0]
        last_slice = occupied_slices[-1]

        return (first_slice + last_slice) // 2


    def dose_volume_reduction_factor(self):  # calculates the dose reduction factor for a slyce
        volume_dose_est = np.array(self.dose_volume_estimate(), dtype=float)

        midpoint = self.breast_midpoint
        midpoint_dose = volume_dose_est[midpoint]

        volume_dose_est = volume_dose_est / midpoint_dose

        red_fac_lyst = [
            self.standard_current/dose if dose != 0 else 0
            for dose in volume_dose_est]

        return red_fac_lyst


    def circle_adder(self, circle_bool=True):  # adds the 100 (circle) or 50 (half circle) current modulations
        # for the MCS for the final circle at the start of the breast
        current_mod_lyst = self.dose_volume_reduction_factor()
        i = len(current_mod_lyst)
        while current_mod_lyst[i - 1] == 0:  # going to the last nonzero value
            i -= 1

        # the idea is to split the final value of the current onto the whole circle/half-circle
        if circle_bool:  # for whole circle
            circle_current = current_mod_lyst[i-1] / 100

            circle_current_lyst = [circle_current for i in range(100)]
            current_mod_lyst[i:i] = circle_current_lyst

            current_mod_lyst[i] = circle_current
            return np.flip(current_mod_lyst)

        # for half circle
        circle_current = current_mod_lyst[i] / 50

        circle_current_lyst = [circle_current for i in range(50)]
        current_mod_lyst[i:i] = circle_current_lyst

        circle_current_lyst[i] = circle_current
        return np.flip(current_mod_lyst)


    def tcm_file_write(self, filename, dirname=None):  # writes the .tcm file with the modulated current
        tube_current_lyst = self.circle_adder()

        if dirname is not None:
            filename = filename.parent / dirname / filename.name
            filename.parent.mkdir(parents=True, exist_ok=True)

        with open(filename, "w") as tcm:
            for fac in tube_current_lyst:
                tcm.write(f"{fac}\n")
            print(f"Wrote {filename} with {len(tube_current_lyst)} lines")


    def tcm_file_geometry(self, filename, dirname=None, linescount=0, geom_bool=True):  # writes the .tcm file with the modulated current
        if linescount == 0:
            linescount = self.z_length + 100 # 100 for the full circle
        tube_current_lyst = np.linspace(0, 32, linescount)
        if not geom_bool: # second option: have a step tcm function
            for i in range(len(tube_current_lyst)):
                tube_current_lyst[i] = 16 * (i * 4 // len(tube_current_lyst)) # * 8 due to current standard
        #tube_current_lyst = np.flip(tube_current_lyst)

        if dirname is not None:
            filename = filename.parent / dirname / filename.name
            filename.parent.mkdir(parents=True, exist_ok=True)

        with open(filename, "w") as tcm:
            for fac in tube_current_lyst:
                tcm.write(f"{fac}\n")
            print(f"Wrote {filename} with {len(tube_current_lyst)} lines")



    # plotting stuff
    @staticmethod
    def plot_data_simple(data):  # simple plot
        x_axis = np.linspace(0, len(data), len(data))
        plt.plot(x_axis, data)
        plt.grid()
        plt.show()



################################
#  small / temporary functions
################################


    @staticmethod
    def showim(data):
        fig, ax = plt.subplots()
        image = ax.imshow(data[0], cmap="cool", vmin=data.min(), vmax=data.max())

        for z in range(data.shape[0]):
            image.set_data(data[z])
            ax.set_title(f"Z = {z}")
            plt.pause(0.1)

        plt.show()

    @staticmethod
    def showimslice(data):
        fig, ax = plt.subplots()
        image = ax.imshow(data, cmap="cool", vmin=data.min(), vmax=data.max())
        plt.show()

# this function is kind of useless since we have segmentation data not dose!!
def dose_mean(data, filter_bool=True):  # filter_bool: filters out data lower than 0.01 or higher than 1000
    if filter_bool:  # filtering out high and low values
        mean = np.mean(data[(data > 0.01) & (data < 1000)])
        print("Average of voxels above 0.01:", mean)
        return mean
    mean = np.mean(data)
    print("Average of all voxels:", mean)
    return mean