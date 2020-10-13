import os
from helpers.projections import compute_projections_given_geometry
from helpers.system_matrix import get_position_of_point_in_HSI_vector, get_position_of_point_in_CTIS_vector
import cv2
import csv
from scipy import sparse
import numpy as np

######################
# File of the main cube>CTIS conversion method
######################

class BadHsiShapeException(Exception):
    pass

class OverlappingProjectionsException(Exception):
    pass

# Spectral sensitivity of the FPA is simply a value by which each cube wavelength will be multiplied.
# A csv specifying these values will be read with this function
def get_spectral_sensitivity(path_sensi, L):
    spectral_sensitivity = []
    try:
        with open(path_sensi) as f:
            reader = csv.reader(f, delimiter=";")
            for line in reader:
                spectral_sensitivity.append(float(line[0]))
    except FileNotFoundError:
        raise Exception("Could not find ",path_sensi)

    # If the csv's number of lines does not correspond to the cube's number of wavelengths, the spectral sensitivity is
    # resized
    if L != len(spectral_sensitivity):
        spectral_sensitivity = [x[0] for x in cv2.resize(np.array(spectral_sensitivity), (1, L))]

    spectral_sensitivity = [x / sum(spectral_sensitivity) for x in spectral_sensitivity]
    return spectral_sensitivity

# Method to convert a hyperspectral cube to a CTIS image
def create_CTIS_from_HSI(HSI, CTIS_geom_params, o0_attenuator, gain, wl_range, L_CTIS, path_save,
                         compute_system_matrix = False, simulate_field_stop_effect = False,
                         path_spectral_sensitivity = None, verbose = False, check_no_projection_overlap = False):
    verboseprint = print if verbose else lambda *a, **k: None

    # Get the cube
    H, W, L = HSI.shape
    mid_H = int(H/2)
    mid_W = int(W/2)

    # Create empty CTIS image
    CTIS = np.zeros((L_CTIS,L_CTIS))
    nb_pixels = L_CTIS**2

    # A CTIS can have different geometries of projections (hexagonal, etc), and different number of orders.
    geometry = CTIS_geom_params[0]
    nb_orders = CTIS_geom_params[1]
    # d0 is the size of the 0th order, which defines the spatio/spectral resolution ratio.
    d0 = CTIS_geom_params[2]

    half_d0 = int(d0 / 2)
    m = int(L_CTIS / 200)  # margins between orders are actually quite reduced. This value is arbitrary.

    measures = [L_CTIS, d0, wl_range, m]

    # Relative energy in each order (based on diffraction theory)
    if nb_orders == 2:
        diffraction_brightnesses = [0.84, 0.13, 0.03]
    else:
        diffraction_brightnesses = [0.84, 0.16]

    # Get spectral sesnitivity if there is one
    if path_spectral_sensitivity is not None:
        spectral_sensitivity = get_spectral_sensitivity(path_spectral_sensitivity,L)
    else:
        spectral_sensitivity = [1 / L for i in range(L)]

    # Get all projections' positions for the given geometry
    verboseprint("Computing projections")
    Projections = compute_projections_given_geometry(geometry, nb_orders, measures, L_CTIS)
    nb_proj = len(Projections)

    # The CTIS action can be seen as a linear transformation between the HSI at full stop (so reduced HSI) and the
    # raw output. It can be represented as a matrix, where each pixel in the raw output is a CL of the reduced HSI pixels.
    # This is called the "system matrix", noted S.
    # THis matrix can be used for cube reconstruction afterwards, so it can be interesting to compute. We've put a flag
    # to do so, because it takes quite a bit of time, and does not need to be computed for each image (once by CTIS
    # geometry is enough).
    # This matrix can be quite big, (nb_voxels in HSI * nb_pixels in CTIS image), so we spent quite a bit of time
    # to find the right tricks to compute and store it: in short, use sparse matrices.
    if compute_system_matrix:
        # initialize empty sparse matrix
        nb_voxels = d0 * d0 * L
        verboseprint("Creating the sparse matrix of size ",nb_pixels," x ", nb_voxels," ie ",
                  nb_pixels * nb_voxels," elements.")
        S = sparse.lil_matrix((nb_pixels, nb_voxels), dtype=np.float32)
    else:
        S = None

   # Actually filing up the projections with cube spectral slices (css)

    # If we are in "field stop" mode, we resize the cube to a reduced version (simulating the optical lens and field stop)
    # where the spatial size of the cube is the one of the 0th order
    if simulate_field_stop_effect:
        HSI_fs = np.stack([cv2.resize(HSI[:,:,x], (d0,d0)) for x in range(L)])
        HSI_fs = np.transpose(HSI_fs, (1, 2, 0))
    # If not, we just keep the original cube. But it must be of the 0th order size, otw an Exeption is raised.
    else:
        if H == d0:
            HSI_fs = HSI.copy()
        else:
            raise BadHsiShapeException("Shape of cube ",H,W," does not correspond to 0th order size ",d0,
                                       ". If you want to resize the cube via an objective lens and a field stop,"
                                       "please use --field_stop mode.")

    # Order 0
    verboseprint("Order 0")

    # We check if the geometry is possible, ie if no projection os overlapping
    if check_no_projection_overlap:
        projection_overlap = False
        test = np.zeros((L_CTIS,L_CTIS,nb_proj+1))

    # Actually creating the 0th order : for each slice of the cube, simply add it to the 0th order position
    for l in range(L):
        top_o0 = mid_H - half_d0
        left_o0 = mid_W - half_d0
        # The final pixel value is the value in the cube times the spectral sensistivy to that wavelength, the
        # relative energy in this order, and the 0th order attenuator effect.
        CTIS[top_o0 : top_o0 + d0, left_o0 : left_o0 + d0] += \
        HSI_fs[:, :, l] * spectral_sensitivity[l] * diffraction_brightnesses[0] * o0_attenuator
        if check_no_projection_overlap:
            test[top_o0 : top_o0 + d0, left_o0 : left_o0 + d0, 0] = 1

    # Add the 0th order's contribution to the system matrix
    if compute_system_matrix:
        for l in range(L):
            top_o0 = mid_H - half_d0
            left_o0 = mid_W - half_d0
            for y_csm in range(d0):
                for x_csm in range(d0):
                    S[get_position_of_point_in_CTIS_vector([y_csm + top_o0, x_csm + left_o0], [H, W]),
                      get_position_of_point_in_HSI_vector([y_csm, x_csm, l], [d0, d0, L])] = \
                        spectral_sensitivity[l] * diffraction_brightnesses[0] * o0_attenuator

    ### Other orders
    for e,p in enumerate(Projections):
        p.reproject_pts_to_ocv_space()
        start_pt = p.start_pt
        end_pt = p.end_pt
        order = p.order
        ys = np.linspace(start_pt[0], end_pt[0], num = L, dtype = np.int)
        xs = np.linspace(start_pt[1], end_pt[1], num = L, dtype = np.int)
        for l in range(L):
            top_current = ys[l] - half_d0
            left_current = xs[l] - half_d0
            CTIS[top_current : top_current+d0, left_current : left_current +d0] +=\
                HSI_fs[:,:,l] * spectral_sensitivity[l] * diffraction_brightnesses[order]
            if check_no_projection_overlap:
                test[top_current : top_current+d0, left_current : left_current +d0,e+1] = 1

            if compute_system_matrix:
                for x_csm in range(d0):
                    for y_csm in range(d0):
                        S[get_position_of_point_in_CTIS_vector([y_csm + top_current, x_csm + left_current], [H,W]),
                          get_position_of_point_in_HSI_vector([y_csm, x_csm, l], [d0, d0, L])] = \
                            spectral_sensitivity[l] * diffraction_brightnesses[order]

    # Check if no projection
    if check_no_projection_overlap:
        test = np.sum(test, axis = -1)

        if np.max(test) >= 2:
            test = test.astype(np.uint8)
            test = cv2.cvtColor(test,cv2.COLOR_GRAY2BGR)
            mask = cv2.inRange(test, (2,2,2),(255,255,255))
            test[mask == 255] = [0, 0, 255]
            test[np.where((test == [1,1,1]).all(axis=2))] = [0, 255, 0]

            cv2.imwrite(os.path.join(path_save,"debug_overlapping_projections.png"),test)
            raise OverlappingProjectionsException("Projections are overlapping. Change geometry, fpa_length, "
                                                  "0th order size and/or wavelength range. "
                                                  "See debug_overlapping_projections.png")

    # Multuply the final image by the gain
    CTIS *= gain
    # Reflect that change in the system matrix
    if compute_system_matrix:
        S = S.multiply(gain)

    return CTIS, S, HSI_fs
