import numpy as np

THR_BREAST = -500
SEGM_CYL = 9


def identify_ring(image_voxels):
    def mean_top_ten(a):
        return np.mean(sorted(a, reverse=True)[:10])

    ring_mask = np.zeros_like(image_voxels)

    num_slices = image_voxels.shape[0]
    # are we highres?
    is_highres = image_voxels.shape[-2] > 1000  # why -2: if 2D equals to [0], if 3D equals to middle axis which is never z
    ring_thickness_px = 24 if is_highres else 14  # TODO improve values
    bkg_thr = -800

    for s in range(num_slices):
        slice = image_voxels[s]
        x_l, x_r, y_t, y_b = None, None, None, None
        # find the outer size of the ring in the x direction
        for i in range(slice.shape[0]):
            # print(s, i, 'mean', mean_top_ten(slice[i, :]))
            if mean_top_ten(slice[i, :]) > bkg_thr:
                x_l = i
                break
        for i in range(slice.shape[0]):
            if mean_top_ten(slice[-i, :]) > bkg_thr:
                x_r = slice.shape[0] - i
                break
        # find the outer size of the ring in the y direction
        for i in range(slice.shape[1]):
            if mean_top_ten(slice[:, i]) > bkg_thr:
                y_t = i
                break
        for i in range(slice.shape[1]):
            if mean_top_ten(slice[:, -i]) > bkg_thr:
                y_b = slice.shape[1] - i
                break
        # calculating center and radius
        if None not in [x_l, x_r, y_t, y_b]:
            center_x, center_y = int((x_r + x_l) / 2), int((y_b + y_t) / 2)
            radius_max = int((x_r - x_l) / 2)
            radius_min = radius_max - ring_thickness_px
            grid_x, grid_y = np.ogrid[:slice.shape[0], :slice.shape[1]]
            r_2d = np.hypot(grid_x - center_x, grid_y - center_y)
            ring_mask[s][(r_2d > radius_min) & (image_voxels[s] > THR_BREAST)] = SEGM_CYL
    return ring_mask

scan = None  # 'scan' is a 3D numpy array

# identify and hide ring
print('Identifying ring...')
ring_mask = identify_ring(scan)
print('Hiding ring...')
scan[ring_mask == SEGM_CYL] = -1024  # replacing ring with air

