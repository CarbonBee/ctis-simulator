import numpy as np
import math

######################
# A projection is one of the connected components of the CTIS image.
# The 0th order has one, the 1st order has several.
# This file regroups all methods linked to these projections, in particular methods which define the position of the
# projections as a function of the desired geometry
######################

class UnknownGeometryException(Exception):
    pass

# The Projection object essentially contains the start point and end point of each projection.
class Projection():

    def __init__(self, start_pt, azimuth, length, order, name = None):
        # s_p is the pixel [y,x] on the CTIS image where the *center* of the first cube spectral slice (css) should be
        # placed.
        # First css <> smallest wl <> closest to the CTIS center
        # e_p is the last point.
        self.name = name
        self.start_pt = start_pt

        self.azimuth = azimuth
        self.proj_angle = proj_angle

        self.length = length
        self.order = order
        self.end_pt = self.compute_end_point_from_start_point()

    def reproject_pts_to_ocv_space(self):
        self.start_pt = [int(mid_H - self.start_pt[0]), int(self.start_pt[1] + mid_W)]
        self.end_pt = [int(mid_H - self.end_pt[0]), int(self.end_pt[1] + mid_W)]

    def compute_end_point_from_start_point(self):
        return [self.start_pt[0] + self.length * np.sin(self.azimuth*np.pi/180),
                self.start_pt[1] + self.length * np.cos(self.azimuth * np.pi / 180)]

    def debug_print(self):
        print("Projection ",self.name," has start pt (y,x) at ",self.start_pt," and end point at ",self.end_pt)

# method to compute the Projections objects depending on the desired geometry
def compute_projections_given_geometry(geom, nb_orders, measures,fpa_length):
    # Coordinates : 0 is the center of the image. X and Y axis in traditional cartesian definitions (*!=* OCV)
    # Degrees rotating following trigonometric circle
    # ie counterclockwise and 0 corresponding to east

    # Getting a bunch of geometric information
    global full_ratio; global mid_H; global half_d0; global d0; global mid_W; global fpa_len
    fpa_len = fpa_length
    [H, d0, wl_range, m] = measures
    full_ratio = 1-wl_range[0]/wl_range[1]
    mid_H = int(H / 2)
    mid_W = mid_H
    half_d0 = d0/2
    global proj_angle
    proj_angle = 45

    # Pre-computing a bunch of distances on the FPA
    full_straight_range = mid_H - half_d0 - 2*m
    ext_straight_range =  (full_straight_range) * full_ratio
    mid_straight_range = (full_straight_range - ext_straight_range - m) * full_ratio
    ext_diagonal_range =  ext_straight_range * np.sqrt(2)
    half_diagonal_angle =  45/2
    ext_halfdiagonal_range =  ext_straight_range / np.cos(half_diagonal_angle*np.pi/180)
    mid_diagonal_range =  mid_straight_range * np.sqrt(2)
    ext_60degrees_range =  ext_straight_range
    mid_60degrees_range =  mid_straight_range * np.sqrt(3)/2
    start_position_ext = mid_H - m - ext_straight_range - half_d0
    start_position_mid = start_position_ext - half_d0 - m - mid_straight_range - half_d0

    # Depending on the geometry, compute the projections
    if geom == "rectangular":
        if nb_orders ==1:
            measures = [ext_straight_range, start_position_ext]
            projections = compute_projections_for_rectangular_o1(measures)
        elif nb_orders==2:
            measures = [start_position_ext, ext_straight_range, ext_diagonal_range, ext_halfdiagonal_range,
                        start_position_mid,mid_straight_range, mid_diagonal_range]
            projections = compute_projections_for_rectangular_o2(measures)
        elif nb_orders>=2:
            raise UnknownGeometryException("Rectangular >=3 currently not correctly implemented")
    elif geom == "hexagonal":
        if nb_orders ==1:
            measures = [start_position_ext, ext_straight_range, ext_60degrees_range]
            projections = compute_projections_for_hexagonal_o1(measures)
        elif nb_orders ==2:
            measures = [start_position_ext, ext_straight_range, ext_60degrees_range, start_position_mid, mid_straight_range,
     mid_60degrees_range]
            projections = compute_projections_for_hexagonal_o2(measures)
        else:
            raise UnknownGeometryException("Unknown geometry ",geom," ",nb_orders)
    elif geom == "radial":
        if nb_orders ==1:
            measures = [full_straight_range, ext_straight_range]
            projections = compute_projections_for_radial_o1(measures)
        else:
            raise UnknownGeometryException("Unknown geometry ",geom," ",nb_orders)
    else:
        raise UnknownGeometryException("Unknown geometry ", geom, " ", nb_orders)
    return projections

def compute_projections_for_radial_o1(measures):
    Projections = []
    [full_straight_range_bordabord, ext_straight_range] = measures

    complete_interior_length = full_straight_range_bordabord - ext_straight_range
    length_of_start_circle = 2 * np.pi * complete_interior_length
    nb_of_possible_rays = math.floor(length_of_start_circle / (d0*1.4))
    angles = [360 / nb_of_possible_rays * i for i in range(nb_of_possible_rays)]
    for angle in angles:
        p = Projection([np.sin(angle * np.pi / 180) * complete_interior_length,
                        np.cos(angle * np.pi / 180) * complete_interior_length],
                       angle, ext_straight_range, 1,
                       "R" + str(int(angle)) + "_1")

        Projections.append(p)
    return Projections

def compute_projections_for_hexagonal_o1(measures):
    Projections = []
    [start_position_ext, ext_straight_range, ext_60degrees_range] = measures

    o1_start_position = start_position_ext
    o1_straight_range = ext_straight_range
    o1_60degrees_range = ext_60degrees_range
    s_e1 = Projection([0, o1_start_position], 0, o1_straight_range, 1, "E1")
    s_w1 = Projection([s_e1.start_pt[0], -s_e1.start_pt[1]], 180, o1_straight_range, 1, "W1")
    Projections.extend([s_e1, s_w1])

    o1_start_position_sqrt32 = np.sqrt(3) / 2 * o1_start_position
    o1_start_position_half = o1_start_position / 2
    s_ne1 = Projection([o1_start_position_sqrt32, o1_start_position_half], 60, o1_60degrees_range, 1, "NE1")
    s_nw1 = Projection([s_ne1.start_pt[0], -s_ne1.start_pt[1]], 120, o1_60degrees_range, 1, "NW1")
    s_sw1 = Projection([-s_ne1.start_pt[0], -s_ne1.start_pt[1]], -120, o1_60degrees_range, 1, "SW1")
    s_se1 = Projection([-s_ne1.start_pt[0], s_ne1.start_pt[1]], -60, o1_60degrees_range, 1, "SE1")
    Projections.extend([s_ne1, s_nw1, s_sw1, s_se1])
    return Projections

def compute_projections_for_hexagonal_o2(measures):
    Projections = []
    [start_position_ext, ext_straight_range, ext_60degrees_range, start_position_mid, mid_straight_range,
     mid_60degrees_range] = measures

    # Order 1
    o1_start_position = start_position_mid
    o1_straight_range = mid_straight_range
    o1_60degrees_range = mid_60degrees_range
    s_e1 = Projection([0, o1_start_position], 0, o1_straight_range, 1, "E1")
    s_w1 = Projection([s_e1.start_pt[0], -s_e1.start_pt[1]], 180, o1_straight_range, 1, "W1")
    Projections.extend([s_e1, s_w1])
    o1_start_position_sqrt32 = np.sqrt(3) / 2 * o1_start_position
    o1_start_position_half = o1_start_position / 2
    s_n3e1 = Projection([o1_start_position_sqrt32, o1_start_position_half], 60, o1_60degrees_range, 1, "N3E1")
    s_n3w1 = Projection([s_n3e1.start_pt[0], -s_n3e1.start_pt[1]], 120, o1_60degrees_range, 1, "N3W1")
    s_s3w1 = Projection([-s_n3e1.start_pt[0], -s_n3e1.start_pt[1]], -120, o1_60degrees_range, 1, "S3W1")
    s_s3e1 = Projection([-s_n3e1.start_pt[0], s_n3e1.start_pt[1]], -60, o1_60degrees_range, 1, "S3E1")
    Projections.extend([s_n3e1, s_n3w1, s_s3w1, s_s3e1])

    # Order 2
    o2_start_position = start_position_ext
    o2_straight_range = ext_straight_range
    o2_60degrees_range = ext_60degrees_range
    o2_start_position_sqrt32 = np.sqrt(3) / 2 * o2_start_position
    o2_start_position_half = o2_start_position / 2
    s_n3e2 = Projection([o2_start_position_sqrt32, o2_start_position_half], 60, o2_60degrees_range, 2, "N3E2")
    s_n3w2 = Projection([s_n3e2.start_pt[0], -s_n3e2.start_pt[1]], 120, o2_60degrees_range, 2, "N3W2")
    s_s3w2 = Projection([-s_n3e2.start_pt[0], -s_n3e2.start_pt[1]], -120, o2_60degrees_range, 2, "S3W2")
    s_s3e2 = Projection([-s_n3e2.start_pt[0], s_n3e2.start_pt[1]], -60, o2_60degrees_range, 2, "S3E2")
    Projections.extend([s_n3e2, s_n3w2, s_s3w2, s_s3e2])

    s_e2 = Projection([0, o2_start_position], 0, o2_straight_range, 2, "E2")
    s_n2 = Projection([o2_start_position, 0], 90, o2_straight_range, 2, "N2")
    s_w2 = Projection([0, -o2_start_position], 180, o2_straight_range, 2, "W2")
    s_s2 = Projection([-o2_start_position, 0], -90, o2_straight_range, 2, "S2")
    Projections.extend([s_n2, s_e2, s_w2, s_s2])
    s_e3n2 = Projection([o2_start_position_half, o2_start_position_sqrt32], 30, o2_60degrees_range, 2, "E3N2")
    s_e3s2 = Projection([-s_e3n2.start_pt[0], s_e3n2.start_pt[1]], -30, o2_60degrees_range, 2, "E3S2")
    s_w3n2 = Projection([s_e3n2.start_pt[0], -s_e3n2.start_pt[1]], 150, o2_60degrees_range, 2, "W3N2")
    s_w3s2 = Projection([-s_e3n2.start_pt[0], -s_e3n2.start_pt[1]], -150, o2_60degrees_range, 2, "W3S2")
    Projections.extend([s_e3n2, s_e3s2, s_w3n2, s_w3s2])

    return Projections

def compute_projections_for_rectangular_o1(measures):
    Projections = []

    [ext_straight_range, start_position_ext] = measures

    o1_straight_range = ext_straight_range
    o1_diagonal_range = ext_straight_range
    o1_start_position = start_position_ext
    s_e1 = Projection([0, o1_start_position], 0, o1_straight_range,1, "E1")
    s_n1 = Projection([o1_start_position,0], 90, o1_straight_range,1, "N1")
    s_w1 = Projection([s_e1.start_pt[0], -s_e1.start_pt[1]], 180, o1_straight_range,1, "W1")
    s_s1 = Projection([-s_n1.start_pt[0], s_n1.start_pt[1]], -90, o1_straight_range,1, "S1")
    Projections.extend([s_e1, s_n1, s_w1, s_s1])
    s_ne1 = Projection([o1_start_position, o1_start_position], 45, o1_diagonal_range,1, "NE1")
    s_nw1 = Projection([s_ne1.start_pt[0],-s_ne1.start_pt[1]], 135,o1_diagonal_range,1, "NW1")
    s_sw1 = Projection([-s_ne1.start_pt[0], -s_ne1.start_pt[1]], -135,o1_diagonal_range,1, "SW1")
    s_se1 = Projection([-s_ne1.start_pt[0], s_ne1.start_pt[1]], -45,o1_diagonal_range,1, "SE1")
    Projections.extend([s_ne1, s_nw1, s_sw1, s_se1])

    return Projections

def compute_projections_for_rectangular_o2(measures):
    Projections = []

    [start_position_ext, ext_straight_range, ext_diagonal_range, ext_halfdiagonal_range, start_position_mid,
     mid_straight_range, mid_diagonal_range] = measures

    #1st order
    o1_start_position = start_position_mid
    o1_straight_range = mid_straight_range
    o1_diagonal_range = mid_diagonal_range
    s_e1 = Projection([0, o1_start_position], 0, o1_straight_range, 1, "E1")
    s_n1 = Projection([o1_start_position, 0], 90, o1_straight_range, 1, "N1")
    s_w1 = Projection([s_e1.start_pt[0], -s_e1.start_pt[1]], 180, o1_straight_range, 1, "W1")
    s_s1 = Projection([-s_n1.start_pt[0], s_n1.start_pt[1]], -90, o1_straight_range, 1, "S1")
    Projections.extend([s_e1, s_n1, s_w1, s_s1])
    s_ne1 = Projection([o1_start_position, o1_start_position], 45, o1_diagonal_range, 1, "NE1")
    s_nw1 = Projection([s_ne1.start_pt[0], -s_ne1.start_pt[1]], 135, o1_diagonal_range, 1, "NW1")
    s_sw1 = Projection([-s_ne1.start_pt[0], -s_ne1.start_pt[1]], -135, o1_diagonal_range, 1, "SW1")
    s_se1 = Projection([-s_ne1.start_pt[0], s_ne1.start_pt[1]], -45, o1_diagonal_range, 1, "SE1")
    Projections.extend([s_ne1, s_nw1, s_sw1, s_se1])

    # 2nd order
    o2_start_position = start_position_ext
    o2_straight_range = ext_straight_range
    o2_diagonal_range = ext_diagonal_range
    o2_halfdiagonal_range = ext_halfdiagonal_range
    s_e2 = Projection([0, o2_start_position], 0, o2_straight_range, 2, "E2")
    s_n2 = Projection([o2_start_position, 0], 90, o2_straight_range, 2, "N2")
    s_w2 = Projection([s_e2.start_pt[0], -s_e2.start_pt[1]], 180, o2_straight_range, 2, "W2")
    s_s2 = Projection([-s_n2.start_pt[0], s_n2.start_pt[1]], -90, o2_straight_range, 2, "S2")
    Projections.extend([s_e2, s_n2, s_w2, s_s2])
    s_ne2 = Projection([o2_start_position, o2_start_position], 45, o2_diagonal_range, 2, "NE2")
    s_nw2 = Projection([s_ne2.start_pt[0], -s_ne2.start_pt[1]], 135, o2_diagonal_range, 2, "NW2")
    s_sw2 = Projection([-s_ne2.start_pt[0], -s_ne2.start_pt[1]], -135, o2_diagonal_range, 2, "SW2")
    s_se2 = Projection([-s_ne2.start_pt[0], s_ne2.start_pt[1]], -45, o2_diagonal_range, 2, "SE2")
    Projections.extend([s_ne2, s_nw2, s_sw2, s_se2])
    # notation based on image of https://fr.wiktionary.org/wiki/nord-nord-est
    o2_start_position_middle = int(o2_start_position/2)
    s_ene2 = Projection([o2_start_position_middle, o2_start_position], 45/2, o2_halfdiagonal_range, 2, "ENE2")
    s_nne2 = Projection([o2_start_position, o2_start_position_middle], 45*3/2, o2_halfdiagonal_range, 2, "NNE2")
    s_nnw2 = Projection([o2_start_position, -o2_start_position_middle], 90+45/2, o2_halfdiagonal_range, 2, "NNW2")
    s_wnw2 = Projection([o2_start_position_middle, - o2_start_position], 90+45*3/2, o2_halfdiagonal_range, 2, "WNW2")
    s_wsw2 = Projection([- s_ene2.start_pt[0], -s_ene2.start_pt[1]], -(90+45*3/2), o2_halfdiagonal_range, 2, "WSW2")
    s_ssw2 = Projection([- s_nne2.start_pt[0], -s_nne2.start_pt[1]], -(90+45/2), o2_halfdiagonal_range, 2, "SSW2")
    s_sse2 = Projection([- s_nnw2.start_pt[0], -s_nnw2.start_pt[1]], -45*3/2, o2_halfdiagonal_range, 2, "SSE2")
    s_ese2 = Projection([- s_wnw2.start_pt[0], -s_wnw2.start_pt[1]], -45/2, o2_halfdiagonal_range, 2, "ESE2")
    Projections.extend([s_ene2, s_nne2, s_nnw2, s_wnw2, s_wsw2, s_ssw2, s_sse2, s_ese2])

    return Projections