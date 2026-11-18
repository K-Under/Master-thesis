import numpy as np
from pathlib import Path

import pydicom

# importing auxillary scripts
import geometry_funcs as gf

"""
This script contains functions to analyze a dicom topogram and generate a TCM file
"""


# loading the file
def load_sim(pathrecon):
    """
    Loads dicom data based on an input folder
    :param pathrecon: path to the recon folder
    :return: data stored in recon folder
    """
    dicom_lyst = []
    for dcm in pathrecon.iterdir():
        # Load DICOM file
        ds = pydicom.dcmread(dcm)
        dicom_lyst.append(ds.pixel_array)

    return np.stack(dicom_lyst)


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
        # self.breast_data = self.ring_remover()  # breast_data (data without ring)
        self.diameter_lyst = self.breast_volume_diameter_estimate()  # breast diameter lyst
        self.breast_length = self.breast_length_estimate()  # length of the breast (mm)
        self.breast_midpoint = self.reference_point_find()  # midpoint index of breast (rounded down)

        # self.breast_mean_dose = self.dose_mean()  # filtered somewhat







    def breast_volume_diameter_estimate(
            self):  # applies the breast_slyce_diameter_estimate function to the whole z-axis of the breast
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
        return (a * np.exp(k1 * x)) + (b * np.exp(k2 * (x ** 2))) + c

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
            self.standard_current / dose if dose != 0 else 0
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
            circle_current = current_mod_lyst[i - 1] / 100

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

    def tcm_file_geometry(self, filename, dirname=None, linescount=0,
                          geom_bool=True):  # writes the .tcm file with the modulated current
        if linescount == 0:
            linescount = self.z_length + 100  # 100 for the full circle
        tube_current_lyst = np.linspace(0, 32, linescount)
        if not geom_bool:  # second option: have a step tcm function
            for i in range(len(tube_current_lyst)):
                tube_current_lyst[i] = 16 * (i * 4 // len(tube_current_lyst))  # * 8 due to current standard
        # tube_current_lyst = np.flip(tube_current_lyst)

        if dirname is not None:
            filename = filename.parent / dirname / filename.name
            filename.parent.mkdir(parents=True, exist_ok=True)

        with open(filename, "w") as tcm:
            for fac in tube_current_lyst:
                tcm.write(f"{fac}\n")
            print(f"Wrote {filename} with {len(tube_current_lyst)} lines")


################################
#  small / temporary functions
################################

# plotting stuff

