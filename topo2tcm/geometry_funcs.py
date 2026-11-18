import numpy as np

"""
This script handles the geometry functions for analyzing the topogram.
I.e. width and length estimation
"""


def slyce2histogram(slyce_data, cutoff=-200, axis=-1, sumtype=0):
    """
    This function creates two histograms from a dicom slice and a cutoff value for HU.
    :param slyce_data: data of the slice
    :param cutoff: cutoff value for HU
    :param axis: -1 both, 0 y, 1 x
    :param sumtype: 0 sums up number of pixels that are valid, 1 sums up HU
    :return: histogram of the summed axes
    """
    if axis == -1: # incase both axis should be summed up, calculate both with the same function
        hist1 = slyce2histogram(slyce_data, axis=0, sumtype=sumtype)
        hist2 = slyce2histogram(slyce_data, axis=1, sumtype=sumtype)
        return hist1, hist2

    mask = slyce_data > cutoff # defining a binary mask by removing all of the values below cutoff
    if sumtype == 0: # summing amount of pixels
        hist = np.sum(mask, axis=axis)
        return hist

    hist = np.sum(np.where(mask, slyce_data, 0), axis=axis) # summing pixel values
    return hist


def width_estimation_hist(hist):
    """
    This function estimates the width of a dicom histogram and a cutoff.
    :param hist: histogram
    :return: width
    """

    hist_max = hist.max()


    return hist_max

