######################
# helper functions to compute the system matrix : conversion between cube/ctis space and vector space
######################

def get_position_of_point_in_HSI_vector(position_in_HSI, HSI_shape):
    # Decomposition in vector will be the following : for a given wl starting from the first one, line by line starting
    # from the top will be added to the vector.
    H, W, L = HSI_shape
    y, x, l = position_in_HSI
    return l * H * W + y * W + x

def get_position_of_point_in_CTIS_vector(position_in_FPA, FPA_shape):
    H, W = FPA_shape
    y, x = position_in_FPA
    return y * W + x
