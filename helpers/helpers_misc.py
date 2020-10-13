import os
import shutil
import cv2
import imageio
import numpy as np

######################
# Misc helper functions
######################

# File management ==========================================
def destroy_and_create(path_):
	if os.path.exists(path_):
		shutil.rmtree(path_)
	os.makedirs(path_)

def create_if_nonexistent(path_):
	if not os.path.exists(path_):
		os.makedirs(path_)

def remove_if_existent(path_):
	if os.path.exists(path_):
		os.remove(path_)

# Image processing =========================================

def get_name_and_ext(im_name):
    im_name_tab = im_name.split(".")
    try:
        return im_name_tab[0], im_name_tab[1]
    except IndexError:
        print("file_im ", im_name)
        raise Exception


# Double security of reading in -1 (original channels) and checking for None
def imread(path_im):
    _, ext = get_name_and_ext(os.path.basename(path_im))

    if ext in ["png", "jpg"]:
        im = cv2.imread(path_im, -1)
    elif ext in ["tif", "tiff"]:
        im = load_3D_im(path_im)
    else:
        print("unknown extension ",ext)
        raise Exception
    if im is None:
        print(path_im, "does not exist")
        raise Exception
    return im

def read_3D_im(path_3Dim):
    if not os.path.exists(path_3Dim):
        print("No file at ",path_3Dim)
        raise Exception
    ext_HSI = path_3Dim.split(".")[-1]
    if ext_HSI == "npy":
        return np.load(path_3Dim)
    elif ext_HSI == "tiff":
        mim = imageio.mimread(path_3Dim)
        L = len(mim)
        [H, W] = mim[0].shape
        cube = np.zeros((H, W, L))
        for i in range(L):
            cube[:, :, i] = mim[i]
    else:
        print("Unknown extension")
        raise Exception
    return cube


def convert_and_write_3D_im(path_im_3D, im_3D):
    # Z ! This converts to uint8 for reading by IJ : may not be the intended behavior
    if path_im_3D.split(".")[-1] not in ["tif", "tiff"]:
        print("Error in saving 3D image : must be saved to tiff file")
        raise Exception
    imageio.mimwrite(path_im_3D, np.transpose(im_3D.astype(np.uint8), (2, 0, 1)))
