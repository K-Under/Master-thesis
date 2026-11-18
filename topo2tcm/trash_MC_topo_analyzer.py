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
    hull = ConvexHull(
        boundary_points)  # gets the surface of the shape as a scipy ConvexHull object (not indices themselves)
    hull_points = boundary_points[
        hull.vertices]  # selects only the point needed to create the surface (shape (points, 2))

    # calculating the diameter
    diff = hull_points[:, None, :] - hull_points[None, :,
                                     :]  # use numpy broadcasting to take the difference of all possibilities
    dist_sq = np.sum(diff ** 2, axis=2)  # s quare each possibility
    diameter = np.sqrt(dist_sq.max())  # take the max value and take the root
    diameter = diameter * self.voxel_size  # convert to mm

    return diameter