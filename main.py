from helpers.helpers_misc import remove_if_existent, destroy_and_create, read_3D_im, imread, create_if_nonexistent,\
convert_and_write_3D_im
from helpers.helpers_args import valid_geometry, get_geom_params
from helpers.cube2ctis import create_CTIS_from_HSI
from scipy import sparse
import argparse
import os
import cv2
import numpy as np

######################
#   Main file
#   Converts a hyperspectral cube into a CTIS image
######################

if __name__ == "__main__":

    # getting arguments from command line
    parser = argparse.ArgumentParser(description='Hyperparameters')
    parser.add_argument('-I', '--input_dir', type=str, required = True,
                        help='input directory of hyperspectral cubes, in tiff format')
    parser.add_argument('-O', '--output_dir', type=str, required = True,
                        help='output directory for CTIS images')

    parser.add_argument('-o0', '--o0_attenuator', type = float, default = 0.1,
                        help = "attenuation coefficient of the 0th order")
    parser.add_argument('-g', '--gain', type = float, default = 50,
                        help = "gain of the CCD")
    parser.add_argument('-w_min', '--min_wavelength', type = int, default = 400,
                        help = "minimum wavelength of the hyperspectral cube")
    parser.add_argument('-w_max', '--max_wavelength', type = int, default = 1000,
                        help = "maximum wavelength of the hyperspectral cube")
    parser.add_argument('-s', '--spectral_sensitivity', type = str, default = None,
                        help = "spectral sensitivity of the sensor csv file")
    parser.add_argument('-gr', '--grating', type = valid_geometry, default = "1R60")

    parser.add_argument('-l', '--fpa_length', type=int, default=512,
                        help="length of the CTIS image")

    parser.add_argument('-m', '--compute_system_matrix', action = "store_true",
                        help = "do you want to compute the system matrix representing the CTIS action")
    parser.add_argument('-fs', '--field_stop', action = "store_true",
                        help = "are you giving in input full size hyperspectral cubes, and do you wish to apply an"
                               "objective lens and a field stop to reduce the cube ?")

    parser.add_argument('-v', '--verbose', action = "store_true")



    args = parser.parse_args()

    # setting paths for loading and saving.
    # This will create a folder <args.output_dir> if non-existent, and create sub-folders inside it to save the
    # outputs
    path_load_cubes = args.input_dir
    path_save_run = args.output_dir
    create_if_nonexistent(path_save_run)
    path_save_CTIS = os.path.join(path_save_run, "CTIS")
    path_save_CTIS_norm = os.path.join(path_save_run, "CTIS_norm")
    destroy_and_create(path_save_CTIS)
    destroy_and_create(path_save_CTIS_norm)
    simulate_field_stop_effect = args.field_stop
    if simulate_field_stop_effect:
        path_save_HSI_fs = os.path.join(path_save_run, "HSI_at_field_stop")
        destroy_and_create(path_save_HSI_fs)

    # optical parameters
    # Length of the FPA : size of the output CTIS image
    fpa_length = args.fpa_length
    # Optical gain : value by which all pixels in the CTIS image are multiplied
    gain = args.gain
    # Optical attenuator to balance the dynamic of 0 and 1 order : value by which all 0th order pixels are multiplied
    o0_attenuator = args.o0_attenuator
    # Should the system matrix be computed ? This can be very long.
    # The system matrix can be very big : <nb_voxels in cube> * <nb_pixels in CTIS image> elements
    compute_system_matrix = args.compute_system_matrix

    # wl range is important because dispersion is linear with wl as d=a*wl, so starting at 400 and stretching to
    # 1000 nec. creates a 40% empty space at the beginning of the spread. See Hagen06
    wl_min = args.min_wavelength
    wl_max = args.max_wavelength
    # Spectral sensitivity of the FPA : values by which each wavelength of the cube will be multiplied (optionnal)
    spectral_sensitivity = args.spectral_sensitivity
    # Diffraction grating parameters : size of 0th order, projections geometry.
    CTIS_geom_params = get_geom_params(args.grating)

    # Print stuff only if in verbose mode
    verboseprint = print if args.verbose else lambda *a, **k: None

    # Going through each cube in <args.input_dir>
    first_cube = True
    for name_HSI in os.listdir(path_load_cubes):

        verboseprint(" ")
        verboseprint("Starting to process ", os.path.join(path_load_cubes, name_HSI))

        # Open the cube
        name = name_HSI.split(".")[0]
        path_HSI = os.path.join(path_load_cubes, name_HSI)
        HSI = read_3D_im(path_HSI)
        H, W, L = HSI.shape

        # Create CTIS image
        CTIS, S, HSI_fs = create_CTIS_from_HSI(HSI = HSI, CTIS_geom_params = CTIS_geom_params, o0_attenuator = o0_attenuator,
                                               gain = gain, compute_system_matrix  = compute_system_matrix and first_cube,
                                               simulate_field_stop_effect = simulate_field_stop_effect,
                                               wl_range = [wl_min, wl_max],
                                               path_spectral_sensitivity = spectral_sensitivity,
                                               verbose = args.verbose,
                                               L_CTIS = fpa_length, check_no_projection_overlap = first_cube,
                                               path_save = path_save_run)


        # Save results
        verboseprint("Saving CTIS at ", os.path.join(path_save_run, name + ".png"))
        cv2.imwrite(os.path.join(path_save_CTIS, name + ".png"), CTIS)
        # For visualization, a normalized version of the CTIS
        CTIS_norm_for_visu = (CTIS-np.min(CTIS))/(np.max(CTIS)-np.min(CTIS))*255
        cv2.imwrite(os.path.join(path_save_CTIS_norm, name+".png"), CTIS_norm_for_visu.astype(np.uint8))

        # If the field stop effect was on, the resulting cube at field stop is saved
        if simulate_field_stop_effect:
            verboseprint("Saving HSI_fs at ", os.path.join(path_save_HSI_fs, name + ".tiff"))
            convert_and_write_3D_im(os.path.join(path_save_HSI_fs, name + ".tiff"), HSI_fs)
        # If the system matrix is computed, save it. It can be very big, so we save it in a sparse format.
        # Use scipy.sparse.load_npz to load it afterwards.
        if S is not None:
            verboseprint("Sparse matrix has nb elements ",S.nnz," which will take up ", S.nnz*2/(10**6),
                  "MB of memory.")
            path_save_S = os.path.join(path_save_run, "system_matrix.npz")
            remove_if_existent(path_save_S)
            sparse.save_npz(path_save_S, sparse.csc_matrix(S))

        first_cube = False

    print("Done.")