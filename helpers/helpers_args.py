import argparse

######################
# helper functions to check and process args from argparse
######################

def valid_geometry(geometry):
    invalid_geometry = False
    if len(geometry)<=2:
        invalid_geometry = True
    else:
        if geometry[0] not in ["1","2"]:
            invalid_geometry = True
        if geometry[1] not in ["R","H","C"]:
            invalid_geometry = True
        try:
            s = int(geometry[2:])
            if s >= 116 or s == 0 or s%2 != 0:
                invalid_geometry = True
        except ValueError:
            invalid_geometry = True

    if invalid_geometry:
        msg = geometry + " is not a valid geometry. This argument must be of the form " \
                         "<1|2><R|C|H><even positive int less than 116>."
        raise argparse.ArgumentTypeError(msg)
    else:
        return geometry

def get_geom_params(param_string):
    codes = {"H" : "hexagonal","R" : "rectangular","C":"radial"}
    return [ codes[param_string[1]], int(param_string[0]), int(param_string[2:])]